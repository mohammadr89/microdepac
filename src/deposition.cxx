/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * INPUTS AND OUTPUTS
 * 
 * === INPUTS (*.ini file) ===
 * [deposition]
 * swdeposition          = boolean : Enable/disable deposition module (default: false)
 * use_depac             = boolean : Use DEPAC model vs. original model (default: true)
 * 
 * # DEPAC-specific parameters (only used when use_depac=true):
 * iratns                = int     : NH3/SO2 ratio regime (required)
 * hlaw                  = float   : Henry's law constant (required)
 * react                 = float   : Reactivity factor (required)
 * c_ave_prev_nh3        = float   : Previous NH3 concentration [ug m-3] (required)
 * nwet_veg              = int     : Vegetation wetness indicator (required)
 * nwet_soil             = int     : Soil wetness indicator (required)
 * nwet_wet              = int     : Wet surface wetness indicator (required)
 * sw_override_ccomp     = boolean : Override compensation point calculation (default: false)
 * ccomp_override_value  = float   : Fixed compensation point value [ug m-3] (default: 0.0)
 * 
 * # Original model parameters (only used when use_depac=false):
 * deposition_var        = float   : Deposition variance (default: 1e5)
 * henry_so2             = float   : SO2 Henry's law constant (default: 1e5)
 * rsoil_so2             = float   : SO2 soil resistance [s m-1] (default: 250.0)
 * rwat_so2              = float   : SO2 water resistance [s m-1] (default: 1.0)
 * rws_so2               = float   : SO2 wet skin resistance [s m-1] (default: 100.0)
 * 
 * [thermo]
 * pbot                  = float   : Surface pressure [Pa] (required for DEPAC)
 * 
 * === OUTPUTS (Cross-sections) ===
 * # Deposition velocities (per tile and grid-mean):
 * vdnh3_veg, vdnh3_soil, vdnh3_wet     : NH3 deposition velocity [m s-1]
 * 
 * # Resistances (DEPAC only):
 * ra, rb                               : Aerodynamic and quasi-laminar resistances [s m-1]
 * ra_veg, ra_soil, ra_wet              : Per-tile aerodynamic resistance [s m-1]
 * rb_veg, rb_soil, rb_wet              : Per-tile quasi-laminar resistance [s m-1]
 * rc_tot, rc_tot_veg, rc_tot_soil, rc_tot_wet : Total canopy resistance [s m-1]
 * rc_eff, rc_eff_veg, rc_eff_soil, rc_eff_wet : Effective canopy resistance [s m-1]
 * cw, cstom, csoil_eff                 : External leaf, stomatal, soil resistances [s m-1]
 * cw_veg, cstom_veg, csoil_eff_veg     : Per-tile resistances [s m-1]
 * 
 * # Compensation points (DEPAC only):
 * ccomp_tot, ccomp_tot_veg, ccomp_tot_soil, ccomp_tot_wet : Total compensation point [ug m-3]
 * cw_out, cstom_out, csoil_out         : Component compensation points [ug m-3]
 * cw_out_veg, cstom_out_veg, etc.      : Per-tile compensation points [ug m-3]
 * 
 * # Surface meteorology (DEPAC only):
 * T_surface, T_surface_veg, T_surface_soil, T_surface_wet     : Surface temperature [K]
 * rh_surface, rh_surface_veg, rh_surface_soil, rh_surface_wet : Surface RH [%]
 * obuk_veg, obuk_soil, obuk_wet        : Obukhov length [m]
 * ustar_veg, ustar_soil, ustar_wet     : Friction velocity [m s-1]
 */

#include "boundary.h"
#include "boundary_surface_lsm.h"
#include "chemistry.h"
#include "constants.h"
#include "cross.h"
#include "deposition.h"
#include "fields.h"
#include "grid.h"
#include "master.h"
#include "netcdf_interface.h"
#include "radiation.h"
#include "stats.h"
#include "thermo.h"
#include "timeloop.h"
#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <utility>
#include "radiation_rrtmgp_functions.h"
#include "radiation_prescribed.h"


// Added: C linkage for DEPAC Fortran wrapper
// The bridge between C++ and Fortran
extern "C" { //This function is written in another language (C/Fortran), but let me use it here in C++
    void depac_wrapper(
            const char* compnam,
            int day_of_year,
            float lat,
            float t,
            float ust, //how strongly the wind transfers momentum to the surface, setting the scale for turbulence and mixing near theground
            float glrad,
            float sinphi,
            float rh,
            float lai,
            float sai,
            int nwet,
            int lu,
            int iratns,
            float *rc_tot,
            float *ccomp_tot,
            float hlaw,
            float react,
            int *status,
            float c_ave_prev_nh3,
            float ra,
            float rb,
            float catm,
            float *rc_eff,
            float *gsoil_eff_out,
            float *rsoil_eff_out,
            float p,  
            float *gw_out,
            float *gstom_out,
            float *cw_out,
            float *cstom_out,
            float *csoil_out,
            bool use_input_ccomp  // flag to use input ccomp_tot
    );
}

namespace {
    template<typename TF>
    void calc_tiled_mean(
        TF* const restrict fld, //This is where the result is stored (an output array).
        const TF* const restrict f_veg, //These are the fractions (between 0 and 1) of each land type at each point.
        const TF* const restrict f_soil,
        const TF* const restrict f_wet,
        const TF* const restrict fld_veg, //These are the values (e.g. NH₃ flux, temperature, etc.) for each tile type.
        const TF* const restrict fld_soil,
        const TF* const restrict fld_wet,
        const TF fac, //A constant multiplier for all final values
        const int istart, const int iend,
        const int jstart, const int jend,
        const int icells)
        {
            for (int j=jstart; j<jend; ++j)
            #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*icells;
                    fld[ij] = (
                            f_veg [ij] * fld_veg [ij] +
                            f_soil[ij] * fld_soil[ij] +
                            f_wet [ij] * fld_wet [ij] ) * fac;
                }
        }

    template<typename TF>
    void calc_vd_water(
        TF* const restrict fld,
        const TF* const restrict ra,
        const TF* const restrict ustar,
        const int* const restrict water_mask,
        const TF diff_scl,
        const TF rwat,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int icells)
        {
            const TF ckarman = 0.4;

            for (int j=jstart; j<jend; ++j)
            #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*icells;

                    if (water_mask[ij] == 1)
                    {
                        const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl;
                        fld[ij] = (TF)1.0 / (ra[ij] + rb + rwat);
                    }
                }
        }

    template<typename TF>
    void calc_spatial_avg_deposition(
        TF* const restrict fld,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int icells)
        {
            //// Calculate sum and count
            //TF n_dep = (TF)0.0;
            //TF sum_dep = (TF)0.0;

            //for (int j=jstart; j<jend; ++j)
            //    #pragma ivdep
            //    for (int i=istart; i<iend; ++i)
            //    {
            //        const int ij = i + j*icells;
            //        sum_dep += fld[ij];
            //        n_dep += 1.0;
            //    }

            //// Calculate and apply average
            //TF avg_dep = sum_dep / n_dep;  // Calculate average across entire domain

            //for (int j=jstart; j<jend; ++j)
            //    #pragma ivdep
            //    for (int i=istart; i<iend; ++i)
            //    {
            //        const int ij = i + j*icells;
            //        fld[ij] = avg_dep; //Replace ALL Values with That Average
            //    }
        }

    template<typename TF>
    void calc_deposition_per_tile_orig(
        const std::basic_string<char> lu_type,
        TF* restrict vdnh3,
        const TF* const restrict lai,
        const TF* const restrict c_veg,
        const TF* const restrict rs,
        const TF* const restrict rs_veg,
        const TF* const restrict ra,
        const TF* const restrict ustar,
        const TF* const restrict fraction,
        const TF* const restrict rmes,
        const TF* const restrict rsoil,
        const TF* const restrict rcut,
        const TF* const restrict rws,
        const TF* const restrict rwat,
        const TF* const restrict diff_scl,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int jj)
    {
        const int ntrac_vd = 1; // One tracer (NH3)
        const TF ckarman = (TF)0.4;
        const TF hc = (TF)10.0; // Canopy height - constant for now...

        // lu_type comes from the boundary/land surface model via
        // boundary.get_tiles() 
        //     ↓ 
        // returns map with keys: "veg", "soil", "wet"
        //     ↓
        // for (auto& tile : tiles)
        //     ↓  
        // tile.first = "veg" (first iteration)
        // tile.first = "soil" (second iteration)  
        // tile.first = "wet" (third iteration)
        //     ↓
        // calc_deposition_per_tile(master, tile.first, ...)
        //     ↓
        // lu_type parameter = "veg", "soil", or "wet"

        // Regarding tile.first:
        // auto& tiles = boundary.get_tiles();
        // // tiles = std::map<std::string, Surface_tile<TF>>
        // //                      ↑              ↑
        // //                   tile.first    tile.second
        // 
        // for (auto& tile : tiles) {
        //     calc_deposition_per_tile<TF>(
        //         master,
        //         tile.first,    // ← This is std::string: "veg", "soil", or "wet"
        //         // ...
        //     );
        // }

        if (lu_type == "veg")
        {
            // Note: I think memory-wise it's more efficient to first loop over ij and then over species,
            // because otherwise rb and rc vectors must be allocated for the entire grid instead of for
            // the number of tracers. Also, it avoids the use of if statements (e.g. "if (t==0) vdnh3[ij] = ...")
            std::vector<TF> rmes_local = {rmes[0]};
            std::vector<TF> rb(ntrac_vd, (TF)0.0); //the vector rb starts filled with zero(s).
            std::vector<TF> rc(ntrac_vd, (TF)0.0);

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i) {
                    const int ij = i + j*jj;

                    //Do not proceed in loop if tile fraction is small
                    if (fraction[ij] < (TF)1e-12)
                        continue;

                    //rmes for NO and NO2 requires multiplication with rs, according to Ganzeveld et al. (1995)
                    const TF ra_inc = (TF)14. * hc * lai[ij] / ustar[ij];

                    for (int t=0; t<ntrac_vd; ++t)
                    {
                        rb[t] = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[t];
                        rc[t] = TF(1.0) / ((TF)1.0 / (diff_scl[t] + rs[ij] + rmes_local[t]) + (TF)1.0 / rcut[t] + (TF)1.0 / (ra_inc + rsoil[t]));
                    }

                    vdnh3[ij]   = (TF)1.0 / (ra[ij] + rb[0] + rc[0]);
                }

        }
        else if (lu_type == "soil")
        {
            std::vector<TF> rb(ntrac_vd, (TF)0.0);

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i) {
                    const int ij = i + j*jj;

                    //Do not proceed in loop if tile fraction is small
                    if (fraction[ij] < (TF)1e-12) continue;

                    for (int t=0; t<ntrac_vd; ++t)
                    {
                        rb[t] = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[t];
                    }
                    vdnh3[ij]   = (TF)1.0 / (ra[ij] + rb[0] + rsoil[0]);
                }
        }
        else if (lu_type == "wet")
        // Mixes contributions from soil and vegetation:
        // Computes two parallel deposition paths:
        //   1. Through the wet leaf surface
        //   2. Through the wet soil surface
        // Blends both paths using c_veg[ij] (canopy cover fraction)
        {
            std::vector<TF> rb_veg(ntrac_vd, (TF)0.0);
            std::vector<TF> rb_soil(ntrac_vd, (TF)0.0);
            std::vector<TF> rc(ntrac_vd, (TF)0.0);
            std::vector<TF> rmes_local = {rmes[0]};

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;

                    // Do not proceed in loop if tile fraction is small
                    if (fraction[ij] < (TF)1e-12) continue;
                    const TF ra_inc = (TF)14. * hc * lai[ij] / ustar[ij];

                    //Note that in rc calculation, rcut is replaced by rws for calculating wet skin uptake
                    for (int t=0; t<ntrac_vd; ++t)
                    {
                        rb_veg[t] = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[t];
                        rb_soil[t] = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[t];
                        rc[t] = TF(1.0) / ((TF)1.0 / (diff_scl[t] + rs_veg[ij] + rmes_local[t]) + (TF)1.0 / rws[t] + (TF)1.0 / (ra_inc + rsoil[t]));
                    }

                    // Calculate vd for wet skin tile as the weighted average of vd to wet soil and to wet vegetation
                    vdnh3[ij]   = c_veg[ij] / (ra[ij] + rb_veg[0] + rc[0]) + ((TF)1.0 - c_veg[ij]) / (ra[ij] + rb_soil[0] + rsoil[0]);
                }
        }
    }

    template<typename TF>
    void calc_deposition_per_tile_depac(
                Master& master,           // Add Master reference as parameter
                const std::string lu_type,
                TF* restrict vdnh3,              // Output: NH3 deposition velocity
                const TF* const restrict lai,     // Leaf Area Index
                const TF* const restrict c_veg,   // Vegetation fraction
                const TF* const restrict rs,      // Surface resistance
                const TF* const restrict rs_veg,  // Vegetation surface resistance
                const TF* const restrict ra,      // Aerodynamic resistance (from MicroHH)
                const TF* const restrict ustar,   // Friction velocity
                const TF* const restrict fraction, // Tile fraction
                const TF* const restrict nh3_concentration,
                const TF* const restrict diff_scl, // Diffusion scaling factors
                const TF* const restrict rho, 
                const TF glrad,          // Global radiation
                const TF sinphi,         // Solar elevation angle
                const TF temperature,    // Temperature [K]
                const TF rh,            // Relative humidity [%]
                const TF* const restrict T_a,
                const TF* const restrict RH_a,
                const TF sai,           // Stem Area Index
                const TF lat,           // Latitude [degrees]
                const int day_of_year,  // Day of year
                const int nwet,         // Surface wetness
    		    const int nwet_veg,      // Add these three parameters
    		    const int nwet_soil,
    		    const int nwet_wet,
                const int lu,           // Land use type
                const int iratns,       // NH3 compensation point option
                const TF hlaw,          // Henry's law constant
                const TF react,         // Reactivity factor
                const TF c_ave_prev_nh3, // Previous NH3 concentration
                const TF catm,          // Atmospheric NH3 concentration
        	    const TF c_ug,          // Concentration conversion factor
                const TF pressure,      // Added pressure parameter
                const bool sw_override_ccomp,        // 
                const TF ccomp_override_value,       // 
                Deposition_tile_map<TF>& deposition_tiles,  // Stores and access DEPAC-specific output data
                                                            // This parameter enables comprehensive DEPAC diagnostics by:
                                                            // Storing all DEPAC intermediate outputs
                                                            // Enabling spatial analysis of resistance components
                                                            // Supporting cross-section outputs for debugging
                                                            // Facilitating model validation and research  
                const Chemistry<TF>& chemistry,        // Chemistry reference
                const int istart, const int iend,
                const int jstart, const int jend,
                const int jj,
                const int kstart,
                const int ijcells)
                {
                    const TF ckarman = 0.4;  // von Karman constant
                    const int STATUS_OK = 0;  // Status code for successful DEPAC calls
                    const TF xmnh3 = 17.031;  // Molar mass of NH3 [g/mol]
                    const TF xmair = 28.9647;       // Molar mass of dry air  [kg kmol-1]
                    const TF xmair_i = TF(1) / xmair;
                    //const TF c_ug = TF(1.0e9) * rhoref[kstart] * xmnh3 * xmair_i;   // mol/mol to ug/m3
                    //const TF c_ug = TF(1.0e9) * rho[kstart] * xmnh3 * xmair_i;   // mol/mol to ug/m3

                    // Define component name outside the loops (doesn't change)
                    char compnam[4] = "NH3";

                    if (lu_type == "veg") {
                        for (int j=jstart; j<jend; ++j)
                            for (int i=istart; i<iend; ++i) {
                                const int ij = i + j*jj;
                                const int ijk = i + j*jj + kstart*ijcells;

                                //if (i == istart && j == jstart) {
                                //    master.print_message("DEPAC call - day_of_year: %d (passing to Fortran)\n", day_of_year);
                                //}

                                //std::cout << "rho= " << rho[kstart] << std::endl;
                                //std::cout <<" c_ave_prev_nh3 * c_ug= " <<c_ave_prev_nh3 * c_ug << std::endl;

                                //std::cout << "VEG tile: i=" << i << ", j=" << j << ", ijk=" << ijk << std::endl;
                                //std::cout << "  NH3 conc = " << nh3_concentration[ijk]
                                //    << ", glrad = " << glrad
                                //    << ", rh = " << rh
                                //    << ", sinphi = " << sinphi << std::endl;

                                if (fraction[ij] < (TF)1e-12)
                                    continue;
                                int local_lu;
                                TF local_sai;

                                if (lai[ij] <= 3.5) {
                                    local_lu = 1;  // grass
                                    local_sai = lai[ij];  // For grass, SAI = LAI
                                } else {
                                    local_lu = 5;  // deciduous forest
                                    local_sai = lai[ij] + 1.0;  // For forest, add stem area
                                }

                                // Keep IFS Ra and use vegetation Rb scaling
                                const TF rb = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[0];

				                // Added this line to store rb
				                deposition_tiles.at(lu_type).rb.data()[ij] = rb;

                                //const TF nh3_ugm3 = nh3_concentration[ijk] * c_ug; // mol/mol to ug/m3 conversion

                                // Use optimal target height concentration from Chemistry module
                                const TF nh3_conc_value = chemistry.get_c_target()[ij];
                                const TF nh3_ugm3 = nh3_conc_value * c_ug; 

                                // debug print for temperature, RH and NH3 concentration passed to depac
                                //std::cout << "Grid points: i=" << i << ", j=" << j
                                //    << ", kstart=" << kstart
                                //    << ", ijk=" << ijk
                                //    << ", NH3=" << nh3_ugm3
                                //    << ", T=" << temperature
                                //    << ", RH=" << rh << std::endl;
                                //std::cout << ", kstart=" << kstart
                                //    << ", NH3=" << nh3_ugm3
                                //    << ", T=" << temperature-273.15
                                //    << ", RH=" << rh << std::endl;

                                // Conductance/resistance variables
                                float rc_tot;           // total canopy resistance Rc (s/m)
                                float ccomp_tot = 0.0;  // total compensation point (ug/m3)
                                float rc_eff;           // effective total canopy resistance (s/m)
                                float gsoil_eff_out;    // effective soil conductance (m/s)
                                float rsoil_eff_out;    // effective soil resistance (s/m) - not directly in DEPAC interface
                                float gw_out;           // external leaf surface conductance (m/s)
                                float gstom_out;        // stomatal conductance (m/s)
                                
                                // Compensation point variables
                                float cw_out;           // external leaf surface compensation point (ug/m3)
                                // Note: cw vs. cw_out! "cw" is inverse of gw_out and is "external leaf surface resistance" 
                                float cstom_out;        // stomatal compensation point (ug/m3)
                                float csoil_out;        // soil compensation point (ug/m3)
                                
                                // Status indicator
                                int status;             // error status (0 = success, >0 = error)

                                // Initialize ccomp_tot with the override value or 0 if no override
                                bool use_input_ccomp = false;

                                // If override is enabled, set the flag and the compensation point value
                                if (sw_override_ccomp) {
                                    ccomp_tot = ccomp_override_value;
                                    use_input_ccomp = true;
                                }

                                depac_wrapper(
                                        compnam,
                                        day_of_year,
                                        lat,
                                        T_a[ij] - 273.15,
                                        ustar[ij],
                                        glrad,
                                        sinphi,
                                        RH_a[ij],
                                        lai[ij],
                                        //sai,
                                        local_sai,      // CHANGED: Use calculated SAI
                                        nwet_veg,       // nwet = 0 for dry vegetation
                                        //lu,
                                        local_lu,       //CHANGED: Use LAI-determined land use type
                                        iratns,
                                        &rc_tot,
                                        &ccomp_tot,
                                        hlaw,
                                        react,
                                        &status,
                                        //c_ave_prev_nh3 * c_ug,
                                        c_ave_prev_nh3,
                                        ra[ij],
                                        rb,
                                        nh3_ugm3,
                                        //catm,
                                        &rc_eff,
                                        &gsoil_eff_out,
                                        &rsoil_eff_out,
                                        pressure,
                                        &gw_out,        
                                        &gstom_out,    
                                        &cw_out,      
                                        &cstom_out,  
                                        &csoil_out,
                                        use_input_ccomp  // Pass the flag
                                );
                                // ccomp_tot = 5.0;
                                // gw_out = 5.0;
                                // gstom_out = 5.0;
                                // gsoil_eff_out = 5.0;
                                // gsoil_eff_out = 5.0;
                                // cw_out; = 5.0;
                                // cstom_out = 5.0;
                                // csoil_out = 5.0;
                                // rc_tot = 80.0;
                                // rc_eff = 80.0;

                                if (status == STATUS_OK) {
                                    deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot; // Store total canopy resistance Rc (s/m)
                                    deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot; // Store ccomp_tot value
                                    deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff; // Store effective total canopy resistance (s/m)
                                    deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? (TF)1.0 / gsoil_eff_out : (TF)9999.0;  // Store effective soil resistance (inverting conductance)
                                    deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? (TF)1.0 / gw_out : (TF)9999.0; // Store external leaf resistance (inverting conductances)
                                    deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? (TF)1.0 / gstom_out : (TF)9999.0; // Store stomatal resistance (inverting conductance)
                                    deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out; //Store external leaf surface compensation point (ug/m3)
                                    deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out; //Store stomatal compensation point (ug/m3)
                                    deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out; //Store stomatal conductance (m/s)

                                    // Calculate deposition velocity using resistance analogy
                                    vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rc_eff);
                                }
                            }
                    }
                    else if (lu_type == "soil") {
                        // Bare soil tile handling  
                        for (int j=jstart; j<jend; ++j)
                            for (int i=istart; i<iend; ++i) {
                                const int ij = i + j*jj;
                                const int ijk = i + j*jj + kstart*ijcells;  // Added this for surface level

                                //std::cout << "VEG tile: i=" << i << ", j=" << j << ", ijk=" << ijk << std::endl;
                                //std::cout << "  NH3 conc = " << nh3_concentration[ijk]
                                //    << ", glrad = " << glrad
                                //    << ", rh = " << rh
                                //    << ", sinphi = " << sinphi << std::endl;

                                if (fraction[ij] < (TF)1e-12)
                                    continue;

                                // Use soil Rb scaling
                                const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[0];

				                // Added this line to store rb
				                deposition_tiles.at(lu_type).rb.data()[ij] = rb;
                                //const TF nh3_ugm3 = nh3_concentration[ijk] * c_ug; // mol/mol to ug/m3 conversion

                                // Use optimal target height concentration from Chemistry module
                                const TF nh3_conc_value = chemistry.get_c_target()[ij];
                                const TF nh3_ugm3 = nh3_conc_value * c_ug;

                                float rc_tot;           // total canopy resistance Rc (s/m)
                                float ccomp_tot = 0.0;  // total compensation point (ug/m3)
                                float rc_eff;           // effective total canopy resistance (s/m)
                                float gsoil_eff_out;    // effective soil conductance (m/s)
                                float rsoil_eff_out;    // effective soil resistance (s/m) - not directly in DEPAC interface
                                float gw_out;           // external leaf surface conductance (m/s)
                                float gstom_out;        // stomatal conductance (m/s)
                                
                                // Compensation point variables
                                float cw_out;           // external leaf surface compensation point (ug/m3)
                                // Note: cw vs. cw_out! "cw" is inverse of gw_out and is "external leaf surface resistance" 
                                float cstom_out;        // stomatal compensation point (ug/m3)
                                float csoil_out;        // soil compensation point (ug/m3)
                                
                                // Status indicator
                                int status;             // error status (0 = success, >0 = error)

                                // Initialize ccomp_tot with the override value or 0 if no override
                                bool use_input_ccomp = false;

                                // If override is enabled, set the flag and the compensation point value
                                if (sw_override_ccomp) {
                                    ccomp_tot = ccomp_override_value;
                                    use_input_ccomp = true;
                                }

                                depac_wrapper(
                                        compnam,
                                        day_of_year,
                                        lat,
                                        T_a[ij] - 273.15,
                                        ustar[ij],
                                        glrad,
                                        sinphi,
                                        RH_a[ij],
                                        lai[ij],
                                        sai,
                                        nwet_soil,  // nwet = 0 for dry soil
                                        lu,
                                        iratns,
                                        &rc_tot,
                                        &ccomp_tot,
                                        hlaw,
                                        react,
                                        &status,
                                        // c_ave_prev_nh3 * c_ug,
                                        c_ave_prev_nh3,
                                        ra[ij],
                                        rb,
                                        nh3_ugm3,
                                        //catm,
                                        &rc_eff,
                                        &gsoil_eff_out,
                                        &rsoil_eff_out,
                                        pressure,
                                        &gw_out,        
                                        &gstom_out,    
                                        &cw_out,      
                                        &cstom_out,  
                                        &csoil_out,
                                        use_input_ccomp  // Pass the flag
                                );

                                if (status == STATUS_OK) {
                                    deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot; // Store total canopy resistance Rc (s/m)
                                    deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot; // Store ccomp_tot value
                                    deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff; // Store effective total canopy resistance (s/m)
                                    deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? (TF)1.0 / gsoil_eff_out : (TF)9999.0;  // Store effective soil resistance (inverting conductance)
                                    deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? (TF)1.0 / gw_out : (TF)9999.0; // Store external leaf resistance (inverting conductances)
                                    deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? (TF)1.0 / gstom_out : (TF)9999.0; // Store stomatal resistance (inverting conductance)
                                    deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out; //Store external leaf surface compensation point (ug/m3)
                                    deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out; //Store stomatal compensation point (ug/m3)
                                    deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out; //Store stomatal conductance (m/s)

                                    vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rsoil_eff_out);
                                }
                            }
                    }
                    else if (lu_type == "wet") {
                        // Wet surfaces handling (both vegetation and soil)
                        for (int j=jstart; j<jend; ++j)
                            for (int i=istart; i<iend; ++i) {
                                const int ij = i + j*jj;
                                const int ijk = i + j*jj + kstart*ijcells;  // Added this for surface level
                                //const TF nh3_ugm3 = nh3_concentration[ijk] * c_ug; // mol/mol to ug/m3 conversion

                                const TF nh3_conc_value = chemistry.get_c_target()[ij];
                                const TF nh3_ugm3 = nh3_conc_value * c_ug;

                                // std::cout << "VEG tile: i=" << i << ", j=" << j << ", ijk=" << ijk << std::endl;
                                // std::cout << "  NH3 conc = " << nh3_concentration[ijk]
                                //     << ", glrad = " << glrad
                                //     << ", rh = " << rh
                                //     << ", sinphi = " << sinphi << std::endl;

                                if (fraction[ij] < (TF)1e-12)
                                    continue;

                                float rc_tot;           // total canopy resistance Rc (s/m)
                                float ccomp_tot = 0.0;  // total compensation point (ug/m3)
                                float rc_eff;           // effective total canopy resistance (s/m)
                                float gsoil_eff_out;    // effective soil conductance (m/s)
                                float rsoil_eff_out;    // effective soil resistance (s/m) - not directly in DEPAC interface
                                float gw_out;           // external leaf surface conductance (m/s)
                                float gstom_out;        // stomatal conductance (m/s)
                                
                                // Compensation point variables
                                float cw_out;           // external leaf surface compensation point (ug/m3)
                                // Note: cw vs. cw_out! "cw" is inverse of gw_out and is "external leaf surface resistance" 
                                float cstom_out;        // stomatal compensation point (ug/m3)
                                float csoil_out;        // soil compensation point (ug/m3)
                                
                                // Status indicator
                                int status;             // error status (0 = success, >0 = error)

                                // Initialize ccomp_tot with the override value or 0 if no override
                                bool use_input_ccomp = false;

                                // If override is enabled, set the flag and the compensation point value
                                if (sw_override_ccomp) {
                                    ccomp_tot = ccomp_override_value;
                                    use_input_ccomp = true;
                                }

                                if (c_veg[ij] > 0) {
                                    // NEW: Added same LAI-based determination for wet vegetation
                                    int local_lu;
                                    TF local_sai;
                                    if (lai[ij] <= 3.5) {
                                        local_lu = 1;  // grass
                                        local_sai = lai[ij];  // For grass, SAI = LAI
                                    } else {
                                        local_lu = 5;  // deciduous forest
                                        local_sai = lai[ij] + 1.0;  // For forest, add stem area
                                    }
                                    // Wet vegetation case

                                    const TF rb = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[0];

				                    // Added this line to store rb
				                    deposition_tiles.at(lu_type).rb.data()[ij] = rb;

                                    depac_wrapper(
                                            compnam,
                                            day_of_year,
                                            lat,
                                            T_a[ij] - 273.15,
                                            ustar[ij],
                                            glrad,
                                            sinphi,
                                            RH_a[ij],
                                            lai[ij],
                                            //sai,
                                            local_sai,      // CHANGED: Use calculated SAI
                                            nwet_wet,       // nwet = 1 for wet conditions
                                            //lu,
                                            local_lu,       // CHANGED: Use LAI-determined land use type
                                            iratns,
                                            &rc_tot,
                                            &ccomp_tot,
                                            hlaw,
                                            react,
                                            &status,
                                            //c_ave_prev_nh3 * c_ug,
                                            c_ave_prev_nh3,
                                            ra[ij],
                                            rb,
                                            nh3_ugm3,
                                            //catm,
                                            &rc_eff,
                                            &gsoil_eff_out,
                                            &rsoil_eff_out,
                                            pressure,
                                            &gw_out,        
                                            &gstom_out,    
                                            &cw_out,      
                                            &cstom_out,  
                                            &csoil_out,
                                            use_input_ccomp  // Pass the flag
                                    );

                                            if (status == STATUS_OK) {
                                                deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot; // Store total canopy resistance Rc (s/m)
                                                deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot; // Store ccomp_tot value
                                                deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff; // Store effective total canopy resistance (s/m)
                                                deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? (TF)1.0 / gsoil_eff_out : (TF)9999.0;  // Store effective soil resistance (inverting conductance)
                                                deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? (TF)1.0 / gw_out : (TF)9999.0; // Store external leaf resistance (inverting conductances)
                                                deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? (TF)1.0 / gstom_out : (TF)9999.0; // Store stomatal resistance (inverting conductance)
                                                deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out; //Store external leaf surface compensation point (ug/m3)
                                                deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out; //Store stomatal compensation point (ug/m3)
                                                deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out; //Store stomatal conductance (m/s)

                                                vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rc_eff);
                                            }
                                }
                                else {
                                    // Wet soil case
                                    const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[0];
				                    // Added this line to store rb
				                    deposition_tiles.at(lu_type).rb.data()[ij] = rb;

                                    depac_wrapper(
                                            compnam,
                                            day_of_year,
                                            lat,
                                            T_a[ij] - 273.15,
                                            ustar[ij],
                                            glrad,
                                            sinphi,
                                            RH_a[ij],
                                            lai[ij],
                                            sai,
                                            nwet_wet,  // nwet = 1 for wet conditions
                                            lu,
                                            iratns,
                                            &rc_tot,
                                            &ccomp_tot,
                                            hlaw,
                                            react,
                                            &status,
                                            //c_ave_prev_nh3 * c_ug,
                                            c_ave_prev_nh3,
                                            ra[ij],
                                            rb,
                                            nh3_ugm3,
                                            //catm,
                                            &rc_eff,
                                            &gsoil_eff_out,
                                            &rsoil_eff_out,
                                            pressure,
                                            &gw_out,          
                                            &gstom_out,      
                                            &cw_out,        
                                            &cstom_out,    
                                            &csoil_out,
                                            use_input_ccomp  // Pass the flag
                                    );

                                    if (status == STATUS_OK) {
                                        deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot; // Store total canopy resistance Rc (s/m)
                                        deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot; // Store ccomp_tot value
                                        deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff; // Store effective total canopy resistance (s/m)
                                        deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? (TF)1.0 / gsoil_eff_out : (TF)9999.0;  // Store effective soil resistance (inverting conductance)
                                        deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? (TF)1.0 / gw_out : (TF)9999.0; // Store external leaf resistance (inverting conductances)
                                        deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? (TF)1.0 / gstom_out : (TF)9999.0; // Store stomatal resistance (inverting conductance)
                                        deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out; //Store external leaf surface compensation point (ug/m3)
                                        deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out; //Store stomatal compensation point (ug/m3)
                                        deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out; //Store stomatal conductance (m/s)

                                        vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rsoil_eff_out);
                                    }
                                }
                            }
                    }
                }

    template<typename TF>
    void calc_deposition_per_tile(
            Master& master,
            const std::string lu_type,
            TF* restrict vdnh3,
            const TF* const restrict lai,
            const TF* const restrict c_veg,
            const TF* const restrict rs,
            const TF* const restrict rs_veg,
            const TF* const restrict ra,
            const TF* const restrict ustar,
            const TF* const restrict fraction,
            const TF* const restrict nh3_concentration,
            const TF* const restrict rmes,
            const TF* const restrict rsoil,
            const TF* const restrict rcut,
            const TF* const restrict rws,
            const TF* const restrict rwat,
            const TF* const restrict diff_scl,
            const TF* const restrict rho,
            const bool use_depac,  // Switch parameter
            // DEPAC-specific parameters below (only used when use_depac is true)
            const TF glrad,
            const TF sinphi,
            const TF temperature,
            const TF rh,
            const TF* const restrict T_a,
            const TF* const restrict RH_a,
            const TF sai,
            const TF lat,
            const int day_of_year,
            const int nwet,
    	    const int nwet_veg,      
    	    const int nwet_soil,
    	    const int nwet_wet,
            const int lu,
            const int iratns,
            const TF hlaw,
            const TF react,
            const TF c_ave_prev_nh3,
            const TF catm,
            const TF c_ug,          // Concentration conversion factor
            const TF pressure,
            const bool sw_override_ccomp,
            const TF ccomp_override_value,
            Deposition_tile_map<TF>& deposition_tiles,
            const Chemistry<TF>& chemistry,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jj,
            const int kstart,
            const int ijcells)
    {
        if (use_depac) {
            // Use DEPAC model
            calc_deposition_per_tile_depac<TF>(
                    master,
                    lu_type,
                    vdnh3,
                    lai,
                    c_veg,
                    rs,
                    rs_veg,
                    ra,
                    ustar,
                    fraction,
                    nh3_concentration,
                    diff_scl,
                    rho,
                    glrad,
                    sinphi,
                    temperature,
                    rh,
                    T_a,
                    RH_a,
                    sai,
                    lat,
                    day_of_year,
                    nwet,
		            nwet_veg, 
    		        nwet_soil,
    		        nwet_wet,
                    lu,
                    iratns,
                    hlaw,
                    react,
                    c_ave_prev_nh3,
                    catm,
		            c_ug,
                    pressure,
                    sw_override_ccomp,
                    ccomp_override_value,
                    deposition_tiles,
                    chemistry,
                    istart, iend,
                    jstart, jend,
                    jj,
                    kstart,
                    ijcells);
        } else {
            // Use original model
            calc_deposition_per_tile_orig(
                    lu_type,
                    vdnh3,
                    lai,
                    c_veg,
                    rs,
                    rs_veg,
                    ra,
                    ustar,
                    fraction,
                    rmes,
                    rsoil,
                    rcut,
                    rws,
                    rwat,
                    diff_scl,
                    istart, iend,
                    jstart, jend,
                    jj);
        }
    }
}

/*
 * The following is a constructor for the Deposition class
 *
 * Constructor - Initializes deposition model and parameters
 * 
 * What the Constructor Does:
 * 1. Reads configuration settings from input file:
 *    - sw_deposition: Boolean flag to enable/disable deposition
 *    - use_depac: Boolean flag to choose between DEPAC and original model
 *
 * 2. Logs which deposition model is being used (DEPAC vs original)
 *
 * 3. Initializes DEPAC-specific parameters by reading them from the input file:
 *    - iratns: NH3/SO2 ratio regime
 *    - hlaw: Henry's law constant  
 *    - react: Reactivity factor
 *    - c_ave_prev_nh3: Previous NH3 concentration
 *    - Various wetness parameters (nwet_veg, nwet_soil, etc.)
 *
 * 4. Sets up override options for compensation points
 */
template<typename TF>
Deposition<TF>::Deposition(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, 
        Radiation<TF>& radiationin, Chemistry<TF>& chemistryin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), radiation(radiationin), chemistry(chemistryin)
{
    sw_deposition = inputin.get_item<bool>("deposition", "swdeposition", "", false);
    use_depac = inputin.get_item<bool>("deposition", "use_depac", "", true);  // Default to DEPAC

    // Log which mode is being used
    if (sw_deposition) {
        if (use_depac) {
            master.print_message("Deposition: Using DEPAC model for NH3 deposition\n");
        } else {
            master.print_message("Deposition: Using original model for NH3 deposition\n");
        }
    }

    // Added: Initialize DEPAC parameters for NH3 deposition
    iratns = inputin.get_item<int>("deposition", "iratns", ""); 
    hlaw = inputin.get_item<TF>("deposition", "hlaw", "");     
    react = inputin.get_item<TF>("deposition", "react", "");                 // Reactivity factor
    //c_ave_prev_nh3 = inputin.get_item<TF>("deposition", "c_ave_prev_nh3", ""); // Previous NH3 concentration (mol/mol, then it converts to ug/m3)
    c_ave_prev_nh3 = inputin.get_item<TF>("deposition", "c_ave_prev_nh3", ""); // Previous NH3 concentration (ug/m3)
    pressure = inputin.get_item<TF>("thermo", "pbot", "");  // Get pressure from thermo settings
    sw_override_ccomp = inputin.get_item<bool>("deposition", "sw_override_ccomp", "", false);
    ccomp_override_value = inputin.get_item<TF>("deposition", "ccomp_override_value", "", TF(0.0));
    nwet_veg = inputin.get_item<int>("deposition", "nwet_veg", "");  // Vegetation wetness
    nwet_soil = inputin.get_item<int>("deposition", "nwet_soil", ""); // Soil wetness
    nwet_wet = inputin.get_item<int>("deposition", "nwet_wet", "");  // Wet surface wetness

    // Read sinphi calculation method (default: prescribed from NetCDF)
    sw_sinphi_prescr = inputin.get_item<bool>("deposition", "sw_sinphi_prescr", "", true);

    if (sw_sinphi_prescr) {
        // Read sunrise/sunset from input NetCDF file
        Netcdf_file input_nc(master, "plume_chem_input.nc", Netcdf_mode::Read);
        
        t_sunrise = input_nc.get_variable<TF>("t_sunrise");
        t_sunset = input_nc.get_variable<TF>("t_sunset");
        
        if (sw_deposition && use_depac) {
            master.print_message("DEPAC: Using prescribed sinphi from NetCDF (sunrise=%.1fh, sunset=%.1fh)\n", 
                                t_sunrise, t_sunset);
        }
    } else {
        if (sw_deposition && use_depac) {
            master.print_message("DEPAC: Using astronomical sinphi\n");
        }
    }

    //// Debug print
    //if (sw_override_ccomp) {
    //    master.print_message("DEPAC: Compensation point override ENABLED. Using value: %f\n", ccomp_override_value);
    //} else {
    //    master.print_message("DEPAC: Compensation point override DISABLED. Using calculated values.\n");
    //}

}

// The destructor for the Deposition class
template <typename TF>
Deposition<TF>::~Deposition()
{
}

template <typename TF>
void Deposition<TF>::init(Input& inputin)
{
    // Read the default deposition velocities. They are needed by 
    // chemistry, even if deposition is disabled.

    // vd_o3   = inputin.get_item<TF>("deposition", "vdo3", "", (TF)0.005);
    // vd_no   = inputin.get_item<TF>("deposition", "vdno", "", (TF)0.002);
    // vd_no2  = inputin.get_item<TF>("deposition", "vdno2", "", (TF)0.005);
    // vd_hno3 = inputin.get_item<TF>("deposition", "vdhno3", "", (TF)0.040);
    // vd_h2o2 = inputin.get_item<TF>("deposition", "vdh2o2", "", (TF)0.018);
    // vd_rooh = inputin.get_item<TF>("deposition", "vdrooh", "", (TF)0.008);
    // vd_hcho = inputin.get_item<TF>("deposition", "vdhcho", "", (TF)0.0033);
    vd_nh3  = inputin.get_item<TF>("deposition", "vdnh3", "", (TF)0.0);  // Added NH3

    if (!sw_deposition)
        return;

    // Get Grid Data Reference
    auto& gd = grid.get_grid_data();

    // Create surface tiles for different surface types for deposition:
    // emplace() creates new Deposition_tile<TF> object for each surface type
    for (auto& name : deposition_tile_names)
        deposition_tiles.emplace(name, Deposition_tile<TF>{});

    // Allocate Memory for Each Tile
    for (auto& tile : deposition_tiles)
    {
        // "second" refers to the value part of a key-value pair in the map/dictionary.

        // tile.second.vdo3.resize(gd.ijcells);
        // tile.second.vdno.resize(gd.ijcells);
        // tile.second.vdno2.resize(gd.ijcells);
        // tile.second.vdhno3.resize(gd.ijcells);
        // tile.second.vdh2o2.resize(gd.ijcells);
        // tile.second.vdrooh.resize(gd.ijcells);
        // tile.second.vdhcho.resize(gd.ijcells);
        tile.second.vdnh3.resize(gd.ijcells);  // allocate memory for arrays
        tile.second.ra.resize(gd.ijcells);
        tile.second.rb.resize(gd.ijcells);
        tile.second.obuk.resize(gd.ijcells);
        tile.second.ustar.resize(gd.ijcells);
        tile.second.T_surface.resize(gd.ijcells);
        tile.second.rh_surface.resize(gd.ijcells);
        // to clarify the line "tile.second.ra.resize(gd.ijcells)"
        // tile.second.ra: 
        //      Type: std::vector<TF>& (reference to a vector)
        //      What it is: Vector/Array - container holding aerodynamic resistance values
        //      Purpose: Stores one resistance value per grid cell
        //      Member of: The Deposition_tile<TF> object

	    if (use_depac){
            tile.second.ccomp_tot.resize(gd.ijcells);
            tile.second.cw.resize(gd.ijcells);
            tile.second.cstom.resize(gd.ijcells);
            tile.second.csoil_eff.resize(gd.ijcells);
            tile.second.cw_out.resize(gd.ijcells);
            tile.second.cstom_out.resize(gd.ijcells);
            tile.second.csoil_out.resize(gd.ijcells);
            tile.second.rc_tot.resize(gd.ijcells);
            tile.second.rc_eff.resize(gd.ijcells);
	    }
    }
    // Initialize grid-mean arrays (Create arrays for tile-averaged values)
    ra_mean.resize(gd.ijcells);
    rb_mean.resize(gd.ijcells);
    std::fill(rb_mean.begin(), rb_mean.end(), TF(0.0));
    T_surface_mean.resize(gd.ijcells);
    rh_surface_mean.resize(gd.ijcells);

    // DEPAC Grid-Mean Arrays
    if (use_depac){
        ccomp_mean.resize(gd.ijcells);
        cw_mean.resize(gd.ijcells);
        cstom_mean.resize(gd.ijcells);
        csoil_eff_mean.resize(gd.ijcells);
        cw_out_mean.resize(gd.ijcells);
        cstom_out_mean.resize(gd.ijcells);
        csoil_out_mean.resize(gd.ijcells);
        rc_tot_mean.resize(gd.ijcells);
        rc_eff_mean.resize(gd.ijcells);
    }

    // Set Tile Descriptions (Human-readable names for output/debugging)
    deposition_tiles.at("veg" ).long_name = "vegetation";
    deposition_tiles.at("soil").long_name = "bare soil";
    deposition_tiles.at("wet" ).long_name = "wet skin";

    // Read Additional Parameters
    deposition_var = inputin.get_item<TF>("deposition", "deposition_var","", (TF)1e5);
    henry_so2 = inputin.get_item<TF>("deposition", "henry_so2", "", (TF)1e5);
    rsoil_so2 = inputin.get_item<TF>("deposition", "rsoil_so2", "", (TF)250.0);
    rwat_so2 = inputin.get_item<TF>("deposition", "rwat_so2", "", (TF)1.0);
    rws_so2 = inputin.get_item<TF>("deposition", "rws_so2", "", (TF)100.0);

    // Note: rmes for NO and NO2 (indices 1 and 2) will still be scaled with rs
    // rmes     = {(TF)1.0, (TF)5.0, (TF)0.5, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    // rsoil    = {(TF)400.0, (TF)1e5, (TF)600.0, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    // rcut     = {(TF)1e5, (TF)1e5, (TF)1e5, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    // rws      = {(TF)2000.0, (TF)1e5, (TF)1e5, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    // rwat     = {(TF)2000.0, (TF)1e5, (TF)1e5, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    // diff     = {(TF)0.13, (TF)0.16, (TF)0.13, (TF)0.11, (TF)0.15, (TF)0.13, (TF)0.16};
    // diff_scl = {(TF)1.6, (TF)1.3, (TF)1.6, (TF)1.9, (TF)1.4, (TF)1.6, (TF)1.3};
    // henry    = {(TF)0.01, (TF)2e-3, (TF)0.01, (TF)1e14, (TF)1e5, (TF)240., (TF)6e3};
    // f0       = {(TF)1.0, (TF)0.0, (TF)0.1, (TF)0.0, (TF)1.0, (TF)0.1, (TF)0.0};

    // Initialize Resistance Arrays for NH3 
    rmes     = {(TF)0.0};  
    rsoil    = {(TF)100.0};
    rcut     = {(TF)1e5};
    rws      = {(TF)10.0};
    rwat     = {(TF)10.0};
    diff     = {(TF)0.13};
    diff_scl = {(TF)1.0};
    // henry    = {(TF)2e4};
    // f0       = {(TF)1.0};

    // Define uninitialized resistance values by scaling with O3 and SO2 resistances (Wesely 1989)
    // for (int i=3; i<7; i++)
    // {
    //     rmes[i]  = (TF)1.0 / (henry[i] / (TF)3000.0 + (TF)100.0 * f0[i]);
    //     rsoil[i] = (TF)1.0 / (henry[i] / (henry_so2 + rsoil_so2) + f0[i] / rsoil[0]);
    //     rcut[i]  = (TF)1.0 / (henry[i] / henry_so2  + f0[i]) * rcut[0];
    //     rws[i]   = (TF)1.0 / (TF(1.0) / ((TF)3.0*rws_so2) + (TF)1e-7 * henry[i] + f0[i] / rws[0]);
    //     rwat[i]  = (TF)1.0 / (henry[i] / (henry_so2 + rwat_so2)  + f0[i] / rwat[0]);
    // }

    // Change diff_scl to diff_scl^(2/3) for use in rb calculation
    for (int i=0; i<1; i++) diff_scl[i] = pow(diff_scl[i], (TF)2.0/(TF)3.0);  // Modified for NH3 only

    // Initialize Tile Deposition Velocities
    for (auto& tile : deposition_tiles)
    {
        // std::fill(tile.second.vdo3.begin(),tile.second.vdo3.end(), vd_o3);
        // std::fill(tile.second.vdno.begin(),tile.second.vdno.end(), vd_no);
        // std::fill(tile.second.vdno2.begin(),tile.second.vdno2.end(), vd_no2);
        // std::fill(tile.second.vdhno3.begin(),tile.second.vdhno3.end(), vd_hno3);
        // std::fill(tile.second.vdh2o2.begin(),tile.second.vdh2o2.end(), vd_h2o2);
        // std::fill(tile.second.vdrooh.begin(),tile.second.vdrooh.end(), vd_rooh);
        // std::fill(tile.second.vdhcho.begin(),tile.second.vdhcho.end(), vd_hcho);
        std::fill(tile.second.vdnh3.begin(),tile.second.vdnh3.end(), vd_nh3);  // Added NH3
    }
}

template <typename TF>
void Deposition<TF>::create(Stats<TF>& stats, Cross<TF>& cross)
{
    if (!sw_deposition)
        return;

    // add cross-sections
    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = {
            // "vdo3_soil", "vdno_soil", "vdno2_soil", "vdhno3_soil", "vdh2o2_soil", "vdrooh_soil", "vdhcho_soil",
            // "vdo3_wet", "vdno_wet", "vdno2_wet", "vdhno3_wet", "vdh2o2_wet", "vdrooh_wet", "vdhcho_wet",
            // "vdo3_veg", "vdno_veg", "vdno2_veg", "vdhno3_veg", "vdh2o2_veg", "vdrooh_veg", "vdhcho_veg"
            "vdnh3_soil", "vdnh3_wet", "vdnh3_veg",
            "ra_soil", "ra_wet", "ra_veg",
            "rb", "rb_soil", "rb_wet", "rb_veg",
            "obuk_soil", "obuk_wet", "obuk_veg",
            "ustar_soil", "ustar_wet", "ustar_veg",
            "ccomp_tot_soil", "ccomp_tot_wet", "ccomp_tot_veg",
            "ra", "ccomp_tot",
            "cw", "cstom", "csoil_eff",
            "cw_soil", "cw_wet", "cw_veg",
            "cstom_soil", "cstom_wet", "cstom_veg",
            "csoil_eff_soil", "csoil_eff_wet", "csoil_eff_veg",

            "cw_out", "cstom_out", "csoil_out", 
            "cw_out_soil", "cw_out_wet", "cw_out_veg",
            "cstom_out_soil", "cstom_out_wet", "cstom_out_veg",
            "csoil_out_soil", "csoil_out_wet", "csoil_out_veg",
            "rc_tot", "rc_tot_veg", "rc_tot_soil", "rc_tot_wet",
            "rc_eff", "rc_eff_veg", "rc_eff_soil", "rc_eff_wet",
            "T_surface",        // grid-mean surface temperature
            "T_surface_veg",    // vegetation tile surface temperature  
            "T_surface_soil",   // soil tile surface temperature
            "T_surface_wet",    // wet tile surface temperature
            "rh_surface",       // grid-mean surface RH
            "rh_surface_veg",   // vegetation tile surface RH
            "rh_surface_soil",  // soil tile surface RH
            "rh_surface_wet"    // wet tile surface RH
        };

        cross_list = cross.get_enabled_variables(allowed_crossvars);
    }
}

template <typename TF>
void Deposition<TF>::update_time_dependent(
        Timeloop<TF>& timeloop,
        Boundary<TF>& boundary,
        Thermo<TF>& thermo,
        // TF* restrict vdo3,
        // TF* restrict vdno,
        // TF* restrict vdno2,
        // TF* restrict vdhno3,
        // TF* restrict vdh2o2,
        // TF* restrict vdrooh,
        // TF* restrict vdhcho,
        TF* restrict vdnh3  // Added NH3
        )
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    const std::vector<TF>& rho = thermo.get_basestate_vector("rho");

    std::vector<TF> T_a(gd.ijcells);   // Temperature array
    std::vector<TF> RH_a(gd.ijcells); // Relative humidity array

    // Only retrieve DEPAC-specific values if using DEPAC
    if (use_depac)
    {
        // Initialize for DEPAC
        day_of_year = int(timeloop.calc_day_of_year());
        lat = gd.lat;
        
        // Calculate sinphi
        if (sw_sinphi_prescr)
        {
            // Prescribed: matches radiation schedule from NetCDF
            const TF hour = timeloop.calc_hour_of_day();
            const TF pi = 3.14159265358979323846;
            const TF dlen = t_sunset - t_sunrise;
            sinphi = std::sin(pi * (hour - t_sunrise) / dlen);
            
        }
        else
        {
            // Astronomical: calculate from date/location
            const int year = timeloop.get_year();
            const TF secs = TF(timeloop.calc_hour_of_day() * 3600);
            TF azimuth;
            std::tie(sinphi, azimuth) = Radiation_rrtmgp_functions::calc_cos_zenith_angle(
                    lat, gd.lon, day_of_year, secs, year);
        }

        // DEBUG: Print hour and sinphi every 100 iterations
        // if (timeloop.get_iteration() % 100 == 0) {
        //     const TF hour = timeloop.calc_hour_of_day();
        //     master.print_message("Hour=%.3f, Day=%d, sinphi=%.6f\n", hour, day_of_year, sinphi);
        // }

        const Radiation_prescribed<TF>& radiation_prescribed = static_cast<const Radiation_prescribed<TF>&>(radiation);

        //master.print_message("DEBUG: About to access radiation, hour = %f\n", timeloop.calc_hour_of_day());

        // Access the first element of the array
        glrad = radiation_prescribed.get_sw_flux_dn()[0];  // Use index 0 to get the correct value
                                                           //master.print_message("DEBUG: Using sw_flux_dn from radiation module as glrad: %f W/m2\n", glrad);

                                                           ///  // Get RH from thermo and convert to %
                                                           ///  auto tmp1 = fields.get_tmp();
                                                           ///  thermo.get_thermo_field(*tmp1, "rh", true, false);
                                                           ///  rh = tmp1->fld.data()[0] * 100.0;
                                                           ///  fields.release_tmp(tmp1);

                                                           ///  // Get temperature from thermo and convert to Celsius
                                                           ///  auto tmp2 = fields.get_tmp();
                                                           ///  thermo.get_thermo_field(*tmp2, "T", true, false);
                                                           ///  temperature = tmp2->fld.data()[0];
                                                           ///  fields.release_tmp(tmp2);

                                                           ///  //// debug prints
                                                           ///  //std::cout << "Temperature from MicroHH (K): " << temperature << std::endl;
                                                           ///  //std::cout << "Temperature passed to DEPAC (C): " << temperature << std::endl;

    // Ask the fields object for a temporary field (a kind of temporary storage or workspace) 
    // and store it in a variable called tmp2.
    // 'fields' is an object that manages field data (like velocity, pressure, etc.)

    // DEBUG: Print every 100 iterations
    if (timeloop.get_iteration() % 100 == 0)
    {
        master.print_message("Iter=%d, Hour=%.3f, sinphi=%.6f, glrad=%.1f W/m2\n", 
                            timeloop.get_iteration(),
                            timeloop.calc_hour_of_day(), 
                            sinphi,
                            glrad);
    }

    auto tmp2 = fields.get_tmp();

    // Send the temporary field tmp2 to the thermo object.
    // Ask thermo to fill that field with the temperature field, "T".

    // 'thermo' is an object that calculates thermodynamic properties (like temperature)
    thermo.get_thermo_field(*tmp2, "T", true, false);
    
    // Extract temperature at kstart level into 2D array
    for (int j = gd.jstart; j < gd.jend; ++j)
        for (int i = gd.istart; i < gd.iend; ++i)
        {
            const int ij = i + j * gd.icells;
            const int ijk = i + j * gd.icells + gd.kstart * gd.ijcells;
            T_a[ij] = tmp2->fld[ijk];
        }
    
    fields.release_tmp(tmp2);
    
    // Extract relative humidity field
    auto tmp1 = fields.get_tmp();
    thermo.get_thermo_field(*tmp1, "rh", true, false);
    
    // Extract RH at kstart level into 2D array
    for (int j = gd.jstart; j < gd.jend; ++j)
        for (int i = gd.istart; i < gd.iend; ++i)
        {
            const int ij = i + j * gd.icells;
            const int ijk = i + j * gd.icells + gd.kstart * gd.ijcells;
            RH_a[ij] = tmp1->fld[ijk] * 100.0;  // Convert to percentage
        }
    
    fields.release_tmp(tmp1);
    
    // Set representative values for synchronization (use first valid grid point)
    temperature = T_a[gd.istart + gd.jstart * gd.icells];
    rh = RH_a[gd.istart + gd.jstart * gd.icells];

    // Store surface temperature and RH in tiles
    for (const auto& tile_name : deposition_tile_names)
    {
        if (deposition_tiles.count(tile_name) > 0)
        {
            auto& dep_tile = deposition_tiles.at(tile_name);
            
            // Copy surface temperature and RH to each tile
            std::copy(T_a.begin(), T_a.end(), dep_tile.T_surface.begin());
            std::copy(RH_a.begin(), RH_a.end(), dep_tile.rh_surface.begin());
        }
    }
    
    // Calculate grid-mean surface temperature and RH
    auto& tiles = boundary.get_tiles();
    get_tiled_mean(T_surface_mean.data(), "T_surface", (TF)1.0, 
                   tiles.at("veg").fraction.data(), 
                   tiles.at("soil").fraction.data(), 
                   tiles.at("wet").fraction.data());
    
    get_tiled_mean(rh_surface_mean.data(), "rh_surface", (TF)1.0,
                   tiles.at("veg").fraction.data(), 
                   tiles.at("soil").fraction.data(), 
                   tiles.at("wet").fraction.data());
    }

    // get information from lsm:
    auto& tiles = boundary.get_tiles();
    auto& lai = boundary.get_lai();
    auto& water_mask = boundary.get_water_mask();
    auto& c_veg = boundary.get_c_veg();
    //auto& sw_flux_dn = radiation.get_c_veg(); //how to call a variable from a remote scheme (maarten)


    // Calculate conversion factor for NH3
    TF xmnh3 = 17.031;  // Molar mass of NH3 [g/mol]
    TF xmair = 28.9647; // Molar mass of dry air [kg kmol-1]
    TF xmair_i = TF(1) / xmair;
    TF c_ug = TF(1.0e9) * fields.rhoref[gd.kstart] * xmnh3 * xmair_i;
    
    // Copy values from boundary tiles to deposition tiles
    for (const auto& tile_name : deposition_tile_names)
    {
        if (tiles.count(tile_name) > 0 && deposition_tiles.count(tile_name) > 0)
        {
            const auto& boundary_tile = tiles.at(tile_name);
            auto& dep_tile = deposition_tiles.at(tile_name);

            std::copy(boundary_tile.ra.begin(), boundary_tile.ra.end(), dep_tile.ra.begin());
            std::copy(boundary_tile.obuk.begin(), boundary_tile.obuk.end(), dep_tile.obuk.begin());
            std::copy(boundary_tile.ustar.begin(), boundary_tile.ustar.end(), dep_tile.ustar.begin());
        }
    }

    // calculate deposition per tile:
    for (auto& tile : tiles)
    {
        calc_deposition_per_tile<TF>(
                master,
                tile.first,
                // deposition_tiles.at(tile.first).vdo3.data(),
                // deposition_tiles.at(tile.first).vdno.data(),
                // deposition_tiles.at(tile.first).vdno2.data(),
                // deposition_tiles.at(tile.first).vdhno3.data(),
                // deposition_tiles.at(tile.first).vdh2o2.data(),
                // deposition_tiles.at(tile.first).vdrooh.data(),
                // deposition_tiles.at(tile.first).vdhcho.data(),

                //// Code explanation... it gets lu_type and what "first" and "second" are!
                    // Function Definition:
                    // template<typename TF>
                    // void calc_deposition_per_tile(
                    //         Master& master,
                    //         const std::string lu_type,  // ← HERE: lu_type is a parameter
                    //         TF* restrict vdnh3,
                    //         // ... other parameters
                    //
                    // 2. Where It Gets Called:
                    // In the update_time_dependent function:
                    // // calculate deposition per tile:
                    // for (auto& tile : tiles)  // ← tiles is a map/container
                    // {
                    //     calc_deposition_per_tile<TF>(
                    //             master,
                    //             tile.first,        // ← THIS becomes lu_type!
                    //             // ... other arguments
                    //
                    // 3. What tile.first Contains:
                    // The tiles container has key-value pairs where:
                    // •	tile.first = the key (string) = "veg", "soil", or "wet"
                    // •	tile.second = the value (tile data)
                    //
                    // 4. Where tiles Comes From:
                    // auto& tiles = boundary.get_tiles();
                    // This gets the tiles from the boundary/land surface model.
                    //
                    // Summary:
                    // tile = pair object
                    // tile.second = tile object
                    // tile.second.ra = vector
                    // resize() = vector method
                    // gd.ijcells = integer value

                deposition_tiles.at(tile.first).vdnh3.data(),  // Added NH3
                lai.data(),
                c_veg.data(),
                tile.second.rs.data(),
                tiles.at("veg").rs.data(),
                tile.second.ra.data(),
                tile.second.ustar.data(),
                tile.second.fraction.data(),
                fields.sp.at("nh3")->fld.data(),  // Pass NH3 concentration directly from Fields
                rmes.data(), 
                rsoil.data(), 
                rcut.data(),
                rws.data(), 
                rwat.data(),
                diff_scl.data(),   
                rho.data(),
                use_depac,  // Add the switch parameter
                // Added: DEPAC parameters
                glrad,          // Now using calculated time-dependent radiation
                sinphi,         // Sine of solar elevation
                temperature,
                rh,
                T_a.data(),    // Air temperature
                RH_a.data(),            // Relative humidity
                sai,           // Stem area index
                lat,           // Latitude
                day_of_year,   // Day of year
                nwet,          // Surface wetness
                nwet_veg,      // Pass the class member variables
                nwet_soil,
                nwet_wet,
                lu,            // Land use type
                iratns,        // NH3 compensation point option
                hlaw,          // Henry's law constant
                react,         // Reactivity factor
                c_ave_prev_nh3, // Previous NH3 concentration
                catm,          // Atmospheric NH3 concentration
		        c_ug,
                pressure,
                sw_override_ccomp,              // NEW argument
                ccomp_override_value,           // NEW argument
                deposition_tiles,               // NEW argument
                chemistry,        // Chemistry reference
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells,
                gd.kstart,           // Added this
                gd.ijcells);         // added this
    }

    // Calculate tile-mean deposition for chemistry
    // get_tiled_mean(vdo3,"o3",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    // get_tiled_mean(vdno,"no",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    // get_tiled_mean(vdno2,"no2",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    // get_tiled_mean(vdhno3,"hno3",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    // get_tiled_mean(vdh2o2,"h2o2",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    // get_tiled_mean(vdrooh,"rooh",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    // get_tiled_mean(vdhcho,"hcho",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdnh3,"nh3",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());  // Added NH3

    // Only calculate DEPAC-specific means if using DEPAC
    if (use_depac) {
        get_tiled_mean(ra_mean.data(), "ra", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(rb_mean.data(), "rb", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(ccomp_mean.data(), "ccomp_tot", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());

        get_tiled_mean(cw_mean.data(), "cw", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(cstom_mean.data(), "cstom", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(csoil_eff_mean.data(), "csoil_eff", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());

        // Added: Calculate grid-mean values for compensation points
        get_tiled_mean(cw_out_mean.data(), "cw_out", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(cstom_out_mean.data(), "cstom_out", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(csoil_out_mean.data(), "csoil_out", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(rc_tot_mean.data(), "rc_tot", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(rc_eff_mean.data(), "rc_eff", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    }

    // cmk: we use the wet-tile info for u* and ra, since these are calculated in lsm with f_wet = 100%
    // update_vd_water(vdo3,"o3",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    // update_vd_water(vdno,"no",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    // update_vd_water(vdno2,"no2",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    // update_vd_water(vdhno3,"hno3",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    // update_vd_water(vdh2o2,"h2o2",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    // update_vd_water(vdrooh,"rooh",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    // update_vd_water(vdhcho,"hcho",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    update_vd_water(vdnh3,"nh3",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());  // Added NH3

    // spatial_avg_vd(vdo3);
    // spatial_avg_vd(vdno);
    // spatial_avg_vd(vdno2);
    // spatial_avg_vd(vdhno3);
    // spatial_avg_vd(vdh2o2);
    // spatial_avg_vd(vdrooh);
    // spatial_avg_vd(vdhcho);
    // spatial_avg_vd(vdnh3);  // Added NH3
}

template<typename TF>
void Deposition<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    const TF no_offset = TF(0);

    for (auto& name : cross_list)
    {
        // if (name == "vdo3_veg")
        //     cross.cross_plane(deposition_tiles.at("veg").vdo3.data(), name, iotime);
        // else if (name == "vdno_veg")
        //     cross.cross_plane(deposition_tiles.at("veg").vdno.data(), name, iotime);
        // else if (name == "vdno2_veg")
        //     cross.cross_plane(deposition_tiles.at("veg").vdno2.data(), name, iotime);
        // else if (name == "vdhno3_veg")
        //     cross.cross_plane(deposition_tiles.at("veg").vdhno3.data(), name, iotime);
        // else if (name == "vdh2o2_veg")
        //     cross.cross_plane(deposition_tiles.at("veg").vdh2o2.data(), name, iotime);
        // else if (name == "vdrooh_veg")
        //     cross.cross_plane(deposition_tiles.at("veg").vdrooh.data(), name, iotime);
        // else if (name == "vdhcho_veg")
        //     cross.cross_plane(deposition_tiles.at("veg").vdhcho.data(), name, iotime);
        // else if (name == "vdo3_soil")
        //     cross.cross_plane(deposition_tiles.at("soil").vdo3.data(), name, iotime);
        // else if (name == "vdno_soil")
        //     cross.cross_plane(deposition_tiles.at("soil").vdno.data(), name, iotime);
        // else if (name == "vdno2_soil")
        //     cross.cross_plane(deposition_tiles.at("soil").vdno2.data(), name, iotime);
        // else if (name == "vdhno3_soil")
        //     cross.cross_plane(deposition_tiles.at("soil").vdhno3.data(), name, iotime);
        // else if (name == "vdh2o2_soil")
        //     cross.cross_plane(deposition_tiles.at("soil").vdh2o2.data(), name, iotime);
        // else if (name == "vdrooh_soil")
        //     cross.cross_plane(deposition_tiles.at("soil").vdrooh.data(), name, iotime);
        // else if (name == "vdhcho_soil")
        //     cross.cross_plane(deposition_tiles.at("soil").vdhcho.data(), name, iotime);
        // else if (name == "vdo3_wet")
        //     cross.cross_plane(deposition_tiles.at("wet").vdo3.data(), name, iotime);
        // else if (name == "vdno_wet")
        //     cross.cross_plane(deposition_tiles.at("wet").vdno.data(), name, iotime);
        // else if (name == "vdno2_wet")
        //     cross.cross_plane(deposition_tiles.at("wet").vdno2.data(), name, iotime);
        // else if (name == "vdhno3_wet")
        //     cross.cross_plane(deposition_tiles.at("wet").vdhno3.data(), name, iotime);
        // else if (name == "vdh2o2_wet")
        //     cross.cross_plane(deposition_tiles.at("wet").vdh2o2.data(), name, iotime);
        // else if (name == "vdrooh_wet")
        //     cross.cross_plane(deposition_tiles.at("wet").vdrooh.data(), name, iotime);
        // else if (name == "vdhcho_wet")
        //     cross.cross_plane(deposition_tiles.at("wet").vdhcho.data(), name, iotime);
        if (name == "vdnh3_veg")
            cross.cross_plane(deposition_tiles.at("veg").vdnh3.data(), no_offset, name, iotime);
        else if (name == "vdnh3_soil")
            cross.cross_plane(deposition_tiles.at("soil").vdnh3.data(), no_offset, name, iotime);
        else if (name == "vdnh3_wet")
            cross.cross_plane(deposition_tiles.at("wet").vdnh3.data(), no_offset, name, iotime);
        //else if (name == "ra_veg")
        //    cross.cross_plane(deposition_tiles.at("veg").ra.data(), no_offset, name, iotime);
        //else if (name == "ra_soil")
        //    cross.cross_plane(deposition_tiles.at("soil").ra.data(), no_offset, name, iotime);
        //else if (name == "ra_wet")
        //    cross.cross_plane(deposition_tiles.at("wet").ra.data(), no_offset, name, iotime);
        else if (name == "rb")
            cross.cross_plane(rb_mean.data(), no_offset, name, iotime);
        else if (name == "rb_veg")
            cross.cross_plane(deposition_tiles.at("veg").rb.data(), no_offset, name, iotime);
        else if (name == "rb_soil")
            cross.cross_plane(deposition_tiles.at("soil").rb.data(), no_offset, name, iotime);
        else if (name == "rb_wet")
            cross.cross_plane(deposition_tiles.at("wet").rb.data(), no_offset, name, iotime);
        //else if (name == "obuk_veg")
        //    cross.cross_plane(deposition_tiles.at("veg").obuk.data(), no_offset, name, iotime);
        //else if (name == "obuk_soil")
        //    cross.cross_plane(deposition_tiles.at("soil").obuk.data(), no_offset, name, iotime);
        //else if (name == "obuk_wet")
        //    cross.cross_plane(deposition_tiles.at("wet").obuk.data(), no_offset, name, iotime);
        //else if (name == "ustar_veg")
        //    cross.cross_plane(deposition_tiles.at("veg").ustar.data(), no_offset, name, iotime);
        //else if (name == "ustar_soil")
        //    cross.cross_plane(deposition_tiles.at("soil").ustar.data(), no_offset, name, iotime);
        //else if (name == "ustar_wet")
        //    cross.cross_plane(deposition_tiles.at("wet").ustar.data(), no_offset, name, iotime);
        else if (name == "ra")
            cross.cross_plane(ra_mean.data(), no_offset, name, iotime);
       else if (use_depac) {
            if (name == "ccomp_tot")
                cross.cross_plane(ccomp_mean.data(), no_offset, name, iotime);
            else if (name == "ccomp_tot_veg")
                cross.cross_plane(deposition_tiles.at("veg").ccomp_tot.data(), no_offset, name, iotime);
            else if (name == "ccomp_tot_soil")
                cross.cross_plane(deposition_tiles.at("soil").ccomp_tot.data(), no_offset, name, iotime);
            else if (name == "ccomp_tot_wet")
                cross.cross_plane(deposition_tiles.at("wet").ccomp_tot.data(), no_offset, name, iotime);
            else if (name == "cw")
                cross.cross_plane(cw_mean.data(), no_offset, name, iotime);
            else if (name == "cstom")
                cross.cross_plane(cstom_mean.data(), no_offset, name, iotime);
            else if (name == "csoil_eff")
                cross.cross_plane(csoil_eff_mean.data(), no_offset, name, iotime);
            else if (name == "cw_out")
                cross.cross_plane(cw_out_mean.data(), no_offset, name, iotime);
            else if (name == "cstom_out")
                cross.cross_plane(cstom_out_mean.data(), no_offset, name, iotime);
            else if (name == "csoil_out")
                cross.cross_plane(csoil_out_mean.data(), no_offset, name, iotime);
            else if (name == "cw_veg")
                cross.cross_plane(deposition_tiles.at("veg").cw.data(), no_offset, name, iotime);
            else if (name == "rc_tot")
                cross.cross_plane(rc_tot_mean.data(), no_offset, name, iotime);
            else if (name == "rc_tot_veg")
                cross.cross_plane(deposition_tiles.at("veg").rc_tot.data(), no_offset, name, iotime);
            else if (name == "rc_tot_soil")
                cross.cross_plane(deposition_tiles.at("soil").rc_tot.data(), no_offset, name, iotime);
            else if (name == "rc_tot_wet")
                cross.cross_plane(deposition_tiles.at("wet").rc_tot.data(), no_offset, name, iotime);
            else if (name == "rc_eff")
                cross.cross_plane(rc_eff_mean.data(), no_offset, name, iotime);
            else if (name == "rc_eff_veg")
                cross.cross_plane(deposition_tiles.at("veg").rc_eff.data(), no_offset, name, iotime);
            else if (name == "rc_eff_soil")
                cross.cross_plane(deposition_tiles.at("soil").rc_eff.data(), no_offset, name, iotime);
            else if (name == "rc_eff_wet")
                cross.cross_plane(deposition_tiles.at("wet").rc_eff.data(), no_offset, name, iotime);
            else if (name == "T_surface")
                cross.cross_plane(T_surface_mean.data(), no_offset, name, iotime);
            else if (name == "T_surface_veg")
                cross.cross_plane(deposition_tiles.at("veg").T_surface.data(), no_offset, name, iotime);
            else if (name == "T_surface_soil")
                cross.cross_plane(deposition_tiles.at("soil").T_surface.data(), no_offset, name, iotime);
            else if (name == "T_surface_wet")
                cross.cross_plane(deposition_tiles.at("wet").T_surface.data(), no_offset, name, iotime);
            else if (name == "rh_surface")
                cross.cross_plane(rh_surface_mean.data(), no_offset, name, iotime);
            else if (name == "rh_surface_veg")
                cross.cross_plane(deposition_tiles.at("veg").rh_surface.data(), no_offset, name, iotime);
            else if (name == "rh_surface_soil")
                cross.cross_plane(deposition_tiles.at("soil").rh_surface.data(), no_offset, name, iotime);
            else if (name == "rh_surface_wet")
                cross.cross_plane(deposition_tiles.at("wet").rh_surface.data(), no_offset, name, iotime);
       }
    }
}


template<typename TF>
const TF Deposition<TF>::get_vd(const std::string& name) const
{
    // if (name == "o3")
    //     return vd_o3;
    // else if (name == "no")
    //     return vd_no;
    // else if (name == "no2")
    //     return vd_no2;
    // else if (name == "hno3")
    //     return vd_hno3;
    // else if (name == "h2o2")
    //     return vd_h2o2;
    // else if (name == "rooh")
    //     return vd_rooh;
    // else if (name == "hcho")
    //     return vd_hcho;
    if (name == "nh3")  // Added NH3
        return vd_nh3;
    else
    {
        std::string error = "Deposition::get_vd() can't return \"" + name + "\"";
        throw std::runtime_error(error);
    }
}

template<typename TF>
void Deposition<TF>::get_tiled_mean(
        TF* restrict fld_out, std::string name, const TF fac,
        const TF* const restrict fveg,
        const TF* const restrict fsoil,
        const TF* const restrict fwet)
//TF* restrict fld_out,           // OUTPUT: Where results go
                                  // - Pre-allocated array for calculated averages
                                  // - 'restrict' = no memory overlap (compiler optimization)
                                  // - Example: vdnh3 array gets filled with NH3 deposition velocities
        
// Example call:
// get_tiled_mean(vdnh3, "nh3", 1.0, veg_fraction, soil_fraction, wet_fraction);
//                  ↑      ↑     ↑         ↑            ↑            ↑
//               fld_out  name  fac      fveg        fsoil        fwet
//
{
    auto& gd = grid.get_grid_data();
    TF* fld_veg;
    TF* fld_soil;
    TF* fld_wet;

    // Yikes..
    // if (name == "o3")
    // {
    //     fld_veg  = deposition_tiles.at("veg").vdo3.data();
    //     fld_soil = deposition_tiles.at("soil").vdo3.data();
    //     fld_wet  = deposition_tiles.at("wet").vdo3.data();
    // }
    // else if (name == "no")
    // {
    //     fld_veg  = deposition_tiles.at("veg").vdno.data();
    //     fld_soil = deposition_tiles.at("soil").vdno.data();
    //     fld_wet  = deposition_tiles.at("wet").vdno.data();
    // }
    // else if (name == "no2")
    // {
    //     fld_veg  = deposition_tiles.at("veg").vdno2.data();
    //     fld_soil = deposition_tiles.at("soil").vdno2.data();
    //     fld_wet  = deposition_tiles.at("wet").vdno2.data();
    // }
    // else if (name == "hno3")
    // {
    //     fld_veg  = deposition_tiles.at("veg").vdhno3.data();
    //     fld_soil = deposition_tiles.at("soil").vdhno3.data();
    //     fld_wet  = deposition_tiles.at("wet").vdhno3.data();
    // }
    // else if (name == "h2o2")
    // {
    //     fld_veg  = deposition_tiles.at("veg").vdh2o2.data();
    //     fld_soil = deposition_tiles.at("soil").vdh2o2.data();
    //     fld_wet  = deposition_tiles.at("wet").vdh2o2.data();
    // }
    // else if (name == "rooh")
    // {
    //     fld_veg  = deposition_tiles.at("veg").vdrooh.data();
    //     fld_soil = deposition_tiles.at("soil").vdrooh.data();
    //     fld_wet  = deposition_tiles.at("wet").vdrooh.data();
    // }
    // else if (name == "hcho")
    // {
    //     fld_veg  = deposition_tiles.at("veg").vdhcho.data();
    //     fld_soil = deposition_tiles.at("soil").vdhcho.data();
    //     fld_wet  = deposition_tiles.at("wet").vdhcho.data();
    // }
    if (name == "nh3")  // Added NH3
    {
        fld_veg  = deposition_tiles.at("veg").vdnh3.data();
        fld_soil = deposition_tiles.at("soil").vdnh3.data();
        fld_wet  = deposition_tiles.at("wet").vdnh3.data();
    }
    else if (name == "ra") {
        fld_veg  = deposition_tiles.at("veg").ra.data();
        fld_soil = deposition_tiles.at("soil").ra.data();
        fld_wet  = deposition_tiles.at("wet").ra.data();
    }
    else if (name == "rb") {
        fld_veg  = deposition_tiles.at("veg").rb.data();
        fld_soil = deposition_tiles.at("soil").rb.data();
        fld_wet  = deposition_tiles.at("wet").rb.data();
    }
    else if (name == "ccomp_tot") {
        fld_veg  = deposition_tiles.at("veg").ccomp_tot.data();
        fld_soil = deposition_tiles.at("soil").ccomp_tot.data();
        fld_wet  = deposition_tiles.at("wet").ccomp_tot.data();
    }
    else if (name == "cw") {
        fld_veg  = deposition_tiles.at("veg").cw.data();
        fld_soil = deposition_tiles.at("soil").cw.data();
        fld_wet  = deposition_tiles.at("wet").cw.data();
    }
    else if (name == "cstom") {
        fld_veg  = deposition_tiles.at("veg").cstom.data();
        fld_soil = deposition_tiles.at("soil").cstom.data();
        fld_wet  = deposition_tiles.at("wet").cstom.data();
    }
    else if (name == "csoil_eff") {
        fld_veg  = deposition_tiles.at("veg").csoil_eff.data();
        fld_soil = deposition_tiles.at("soil").csoil_eff.data();
        fld_wet  = deposition_tiles.at("wet").csoil_eff.data();
    }
    else if (name == "cw_out") {
        fld_veg  = deposition_tiles.at("veg").cw_out.data();
        fld_soil = deposition_tiles.at("soil").cw_out.data();
        fld_wet  = deposition_tiles.at("wet").cw_out.data();
    }
    else if (name == "cstom_out") {
        fld_veg  = deposition_tiles.at("veg").cstom_out.data();
        fld_soil = deposition_tiles.at("soil").cstom_out.data();
        fld_wet  = deposition_tiles.at("wet").cstom_out.data();
    }
    else if (name == "csoil_out") {
        fld_veg  = deposition_tiles.at("veg").csoil_out.data();
        fld_soil = deposition_tiles.at("soil").csoil_out.data();
        fld_wet  = deposition_tiles.at("wet").csoil_out.data();
    }
    else if (name == "rc_tot") {
        fld_veg  = deposition_tiles.at("veg").rc_tot.data();
        fld_soil = deposition_tiles.at("soil").rc_tot.data();
        fld_wet  = deposition_tiles.at("wet").rc_tot.data();
    }
    else if (name == "rc_eff") {
        fld_veg  = deposition_tiles.at("veg").rc_eff.data();
        fld_soil = deposition_tiles.at("soil").rc_eff.data();
        fld_wet  = deposition_tiles.at("wet").rc_eff.data();
    }
    else if (name == "T_surface") { 
        fld_veg  = deposition_tiles.at("veg").T_surface.data();
        fld_soil = deposition_tiles.at("soil").T_surface.data();
        fld_wet  = deposition_tiles.at("wet").T_surface.data();
    }
    else if (name == "rh_surface") {
        fld_veg  = deposition_tiles.at("veg").rh_surface.data();
        fld_soil = deposition_tiles.at("soil").rh_surface.data();
        fld_wet  = deposition_tiles.at("wet").rh_surface.data();
    }
    else
        throw std::runtime_error("Cannot calculate tiled mean for variable \"" + name + "\"\\n");

    calc_tiled_mean(
            fld_out,
            fveg,
            fsoil,
            fwet,
            fld_veg,
            fld_soil,
            fld_wet,
            fac,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
}

template<typename TF>
void Deposition<TF>::update_vd_water(
        TF* restrict fld_out, std::string name,
        const TF* const restrict ra,
        const TF* const restrict ustar,
        const int* const restrict water_mask,
        const TF* const restrict diff_scl,
        const TF* const restrict rwat)
{
    auto& gd = grid.get_grid_data();

    // TF* fld;
    TF diff_scl_val;
    TF rwat_val;

    // Yikes...
    // if (name == "o3")
    // {
    //     // fld = vd_o3.data();
    //     diff_scl_val = diff_scl[0];
    //     rwat_val = rwat[0];
    // }
    // else if (name == "no")
    // {
    //     // fld = vd_no.data();
    //     diff_scl_val = diff_scl[1];
    //     rwat_val = rwat[1];
    // }
    // else if (name == "no2")
    // {
    //     // fld = vd_no2.data();
    //     diff_scl_val = diff_scl[2];
    //     rwat_val = rwat[2];
    // }
    // else if (name == "hno3")
    // {
    //     // fld = vd_hno3.data();
    //     diff_scl_val = diff_scl[3];
    //     rwat_val = rwat[3];
    // }
    // else if (name == "h2o2")
    // {
    //     // fld = vd_h2o2.data();
    //     diff_scl_val = diff_scl[4];
    //     rwat_val = rwat[4];
    // }
    // else  if (name == "rooh")
    // {
    //     // fld = vd_rooh.data();
    //     diff_scl_val = diff_scl[5];
    //     rwat_val = rwat[5];
    // }
    // else if (name == "hcho")
    // {
    //     // fld = vd_hcho.data();
    //     diff_scl_val = diff_scl[6];
    //     rwat_val = rwat[6];
    // }
    if (name == "nh3")  // Added NH3
    {
        diff_scl_val = diff_scl[0];  // Using first index for NH3
        rwat_val = rwat[0];  // Using first index for NH3
    }
    else
        throw std::runtime_error("Cannot update vd to water for variable \"" + name + "\"\\n");

    calc_vd_water(
            fld_out,
            ra,
            ustar,
            water_mask,
            diff_scl_val,
            rwat_val,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
}

template<typename TF>
void Deposition<TF>::spatial_avg_vd(
        TF* restrict fld_out)
{
    auto& gd = grid.get_grid_data();

    calc_spatial_avg_deposition(
            fld_out,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
}

template class Deposition<double>;
//:template class Chemistry<float>;

