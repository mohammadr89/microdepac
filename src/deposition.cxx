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
#include <cmath>
#include "radiation.h"
#include "radiation_prescribed.h"
#include "radiation_rrtmgp_functions.h"
#include "stats.h"
#include "thermo.h"
#include "timeloop.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <utility>


// Added: C linkage for DEPAC Fortran wrapper
// The bridge between C++ and Fortran
extern "C"
{
    void depac_wrapper
        (
        const char* compnam,
        int day_of_year,
        float lat,
        float t,
        float ust, 
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
        bool use_input_ccomp 
        );
}

namespace
{
    template<typename TF>
    void calc_tiled_mean
        (
        TF* const restrict fld, 
        const TF* const restrict f_veg, 
        const TF* const restrict f_soil,
        const TF* const restrict f_wet,
        const TF* const restrict fld_veg, 
        const TF* const restrict fld_soil,
        const TF* const restrict fld_wet,
        const TF fac, 
        const int istart, const int iend,
        const int jstart, const int jend,
        const int icells
        )
    {
        for (int j=jstart; j<jend; ++j)
        #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                fld[ij] =
                (f_veg [ij] * fld_veg [ij] +
                 f_soil[ij] * fld_soil[ij] +
                 f_wet [ij] * fld_wet [ij] ) * fac;
            }
    }

    template<typename TF>
    void calc_vd_water
        (
        TF* const restrict fld,
        const TF* const restrict ra,
        const TF* const restrict ustar,
        const int* const restrict water_mask,
        const TF diff_scl,
        const TF rwat,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int icells
        )
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
    void calc_deposition_per_tile_orig
        (
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
        const int jj
        )
    {
        const int ntrac_vd = 1; 
        const TF ckarman = (TF)0.4;
        const TF hc = (TF)10.0; 
        if (lu_type == "veg")
        {
            std::vector<TF> rmes_local = {rmes[0]};
            std::vector<TF> rb(ntrac_vd, (TF)0.0); 
            std::vector<TF> rc(ntrac_vd, (TF)0.0);

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
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
                for (int i=istart; i<iend; ++i)
                {
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
                    if (fraction[ij] < (TF)1e-12) continue;
                    const TF ra_inc = (TF)14. * hc * lai[ij] / ustar[ij];
                    for (int t=0; t<ntrac_vd; ++t)
                    {
                        rb_veg[t] = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[t];
                        rb_soil[t] = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[t];
                        rc[t] = TF(1.0) / ((TF)1.0 / (diff_scl[t] + rs_veg[ij] + rmes_local[t]) + (TF)1.0 / rws[t] + (TF)1.0 / (ra_inc + rsoil[t]));
                    }
                    vdnh3[ij]   = c_veg[ij] / (ra[ij] + rb_veg[0] + rc[0]) + ((TF)1.0 - c_veg[ij]) / (ra[ij] + rb_soil[0] + rsoil[0]);
                }
        }
    }

    template<typename TF>
    void calc_deposition_per_tile_depac
        (
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
        const TF* const restrict diff_scl,
        const TF* const restrict rho,
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
        const TF c_ug,
        const TF pressure,
        const bool sw_override_ccomp,
        const TF ccomp_override_value,
        Deposition_tile_map<TF>& deposition_tiles,
        const Chemistry<TF>& chemistry,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int jj,
        const int* const restrict lu_map,
        const int kstart,
        const int ijcells
        )
    {
        const TF ckarman = 0.4;
        const int STATUS_OK = 0;
        const TF xmnh3 = 17.031;
        const TF xmair = 28.9647;
        const TF xmair_i = TF(1) / xmair;
        char compnam[4] = "NH3";

        if (lu_type == "veg") 
        {
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + kstart*ijcells;
                    if (fraction[ij] < (TF)1e-12)
                        continue;

                    const int local_lu = lu_map[ij];
                    // Forest types in DEPAC: 4=coniferous, 5=deciduous
                    const bool is_forest = (local_lu == 4 || local_lu == 5 || local_lu == 12);
                    const TF local_sai   = is_forest ? lai[ij] + TF(1.0) : lai[ij];

                    // Keep IFS Ra and use vegetation Rb scaling
                    const TF rb = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[0];
                    deposition_tiles.at(lu_type).rb.data()[ij] = rb;

                    // Use optimal target height concentration from Chemistry module
                    const TF nh3_conc_value = chemistry.get_c_target()[ij];
                    const TF nh3_ugm3 = nh3_conc_value * c_ug; 
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
                    int status;             // error status (0 = success, >0 = error)
                    // Initialize ccomp_tot with the override value or 0 if no override
                    bool use_input_ccomp = false;
                    // If override is enabled, set the flag and the compensation point value
                    if (sw_override_ccomp)
                    {
                        ccomp_tot = ccomp_override_value;
                        use_input_ccomp = true;
                    }

                    depac_wrapper
                        (
                        compnam,
                        day_of_year,
                        lat,
                        T_a[ij] - 273.15,
                        ustar[ij],
                        glrad,
                        sinphi,
                        RH_a[ij],
                        lai[ij],
                        local_sai,
                        nwet_veg,
                        local_lu,
                        iratns,
                        &rc_tot,
                        &ccomp_tot,
                        hlaw,
                        react,
                        &status,
                        c_ave_prev_nh3,
                        ra[ij],
                        rb,
                        nh3_ugm3,
                        &rc_eff,
                        &gsoil_eff_out,
                        &rsoil_eff_out,
                        pressure,
                        &gw_out,        
                        &gstom_out,    
                        &cw_out,      
                        &cstom_out,  
                        &csoil_out,
                        use_input_ccomp
                        );

                    if (status == STATUS_OK)
                    {
                        deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot;
                        deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot;
                        deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff;
                        deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? (TF)1.0 / gsoil_eff_out : (TF)9999.0;
                        deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? (TF)1.0 / gw_out : (TF)9999.0;
                        deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? (TF)1.0 / gstom_out : (TF)9999.0;
                        deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out;
                        deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out;
                        deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out;

                        const TF total_resistance = ra[ij] + rb + rc_eff;
                        if (std::abs(total_resistance) > (TF)1e-6)
                        {
                            vdnh3[ij] = (TF)1.0 / total_resistance;
                            // Note: vdnh3 can be negative (emission) if total_resistance is negative
                        }
                        else
                        {
                            vdnh3[ij] = (TF)0.0;  // Only set to zero if resistance is very close to zero
                        }
                    }
                }
        }
        else if (lu_type == "soil")
        {
            // Bare soil tile handling  
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + kstart*ijcells;  

                    if (fraction[ij] < (TF)1e-12)
                        continue;

                    const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[0];
                    deposition_tiles.at(lu_type).rb.data()[ij] = rb;
                    const TF nh3_conc_value = chemistry.get_c_target()[ij];
                    const TF nh3_ugm3 = nh3_conc_value * c_ug;
                    float rc_tot;
                    float ccomp_tot = 0.0;
                    float rc_eff;
                    float gsoil_eff_out;
                    float rsoil_eff_out;
                    float gw_out;
                    float gstom_out;
                    float cw_out;
                    float cstom_out;
                    float csoil_out;
                    int status;
                    bool use_input_ccomp = false;

                    if (sw_override_ccomp)
                    {
                        ccomp_tot = ccomp_override_value;
                        use_input_ccomp = true;
                    }

                    depac_wrapper
                        (
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
                        nwet_soil,
                        lu,
                        iratns,
                        &rc_tot,
                        &ccomp_tot,
                        hlaw,
                        react,
                        &status,
                        c_ave_prev_nh3,
                        ra[ij],
                        rb,
                        nh3_ugm3,
                        &rc_eff,
                        &gsoil_eff_out,
                        &rsoil_eff_out,
                        pressure,
                        &gw_out,        
                        &gstom_out,    
                        &cw_out,      
                        &cstom_out,  
                        &csoil_out,
                        use_input_ccomp
                        );

                    if (status == STATUS_OK)
                    {
                        deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot;
                        deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot;
                        deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff;
                        deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? (TF)1.0 / gsoil_eff_out : (TF)9999.0;
                        deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? (TF)1.0 / gw_out : (TF)9999.0;
                        deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? (TF)1.0 / gstom_out : (TF)9999.0;
                        deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out;
                        deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out;
                        deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out;

                        const TF total_resistance = ra[ij] + rb + rsoil_eff_out;
                        if (std::abs(total_resistance) > (TF)1e-6)
                        {
                            vdnh3[ij] = (TF)1.0 / total_resistance;
                        }
                        else
                        {
                            vdnh3[ij] = (TF)0.0;
                        }
                    }
                }
        }
        else if (lu_type == "wet")
        {
            // Wet surfaces handling (both vegetation and soil)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + kstart*ijcells; 
                    const TF nh3_conc_value = chemistry.get_c_target()[ij];
                    const TF nh3_ugm3 = nh3_conc_value * c_ug;
                    if (fraction[ij] < (TF)1e-12)
                        continue;
                    float rc_tot;
                    float ccomp_tot = 0.0;
                    float rc_eff;
                    float gsoil_eff_out;
                    float rsoil_eff_out;
                    float gw_out;
                    float gstom_out;
                    float cw_out;
                    float cstom_out;
                    float csoil_out;
                    int status;
                    bool use_input_ccomp = false;
                    if (sw_override_ccomp)
                    {
                        ccomp_tot = ccomp_override_value;
                        use_input_ccomp = true;
                    }
                    if (c_veg[ij] > 0)
                    {
                        const int local_lu = lu_map[ij];
                        const bool is_forest = (local_lu == 4 || local_lu == 5 || local_lu == 12);
                        const TF local_sai   = is_forest ? lai[ij] + TF(1.0) : lai[ij];
                        const TF rb = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[0];
                        deposition_tiles.at(lu_type).rb.data()[ij] = rb;

                        depac_wrapper
                            (
                            compnam,
                            day_of_year,
                            lat,
                            T_a[ij] - 273.15,
                            ustar[ij],
                            glrad,
                            sinphi,
                            RH_a[ij],
                            lai[ij],
                            local_sai,
                            nwet_wet,
                            local_lu,
                            iratns,
                            &rc_tot,
                            &ccomp_tot,
                            hlaw,
                            react,
                            &status,
                            c_ave_prev_nh3,
                            ra[ij],
                            rb,
                            nh3_ugm3,
                            &rc_eff,
                            &gsoil_eff_out,
                            &rsoil_eff_out,
                            pressure,
                            &gw_out,        
                            &gstom_out,    
                            &cw_out,      
                            &cstom_out,  
                            &csoil_out,
                            use_input_ccomp
                            );

                        if (status == STATUS_OK)
                        {
                            deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot;
                            deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot;
                            deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff;
                            deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? (TF)1.0 / gsoil_eff_out : (TF)9999.0;
                            deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? (TF)1.0 / gw_out : (TF)9999.0;
                            deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? (TF)1.0 / gstom_out : (TF)9999.0;
                            deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out;
                            deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out;
                            deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out;
                            const TF total_resistance = ra[ij] + rb + rc_eff;
                            if (std::abs(total_resistance) > (TF)1e-6)
                            {
                                vdnh3[ij] = (TF)1.0 / total_resistance;
                            }
                            else
                            {
                                vdnh3[ij] = (TF)0.0;
                            }
                        }
                    }
                    else
                    {
                        const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[0];
                        deposition_tiles.at(lu_type).rb.data()[ij] = rb;

                        depac_wrapper
                            (
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
                            nwet_wet,
                            lu,
                            iratns,
                            &rc_tot,
                            &ccomp_tot,
                            hlaw,
                            react,
                            &status,
                            c_ave_prev_nh3,
                            ra[ij],
                            rb,
                            nh3_ugm3,
                            &rc_eff,
                            &gsoil_eff_out,
                            &rsoil_eff_out,
                            pressure,
                            &gw_out,          
                            &gstom_out,      
                            &cw_out,        
                            &cstom_out,    
                            &csoil_out,
                            use_input_ccomp
                            );

                        if (status == STATUS_OK)
                        {
                            deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot;
                            deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot;
                            deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff;
                            deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? (TF)1.0 / gsoil_eff_out : (TF)9999.0;
                            deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? (TF)1.0 / gw_out : (TF)9999.0;
                            deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? (TF)1.0 / gstom_out : (TF)9999.0;
                            deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out;
                            deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out;
                            deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out;
	                        const TF total_resistance = ra[ij] + rb + rsoil_eff_out;
	                        if (std::abs(total_resistance) > (TF)1e-6)
	                        {
	                            vdnh3[ij] = (TF)1.0 / total_resistance;
	                        }
                            else
                            {
                                vdnh3[ij] = (TF)0.0;
	                        }
                        }
                    }
                }
        }
    }

    template<typename TF>
    void calc_deposition_per_tile
        (
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
        const bool use_depac,
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
        const TF c_ug,
        const TF pressure,
        const bool sw_override_ccomp,
        const TF ccomp_override_value,
        Deposition_tile_map<TF>& deposition_tiles,
        const Chemistry<TF>& chemistry,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int jj,
        const int* const restrict lu_map,
        const int kstart,
        const int ijcells
        )
    {
        if (use_depac)
        {
            calc_deposition_per_tile_depac<TF>
                (
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
                lu_map,
                kstart,
                ijcells
                );
        }
        else
        {
            calc_deposition_per_tile_orig
                (
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
                jj
                );
        }
    }
}

template<typename TF>
Deposition<TF>::Deposition(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Radiation<TF>& radiationin, Chemistry<TF>& chemistryin, Input& inputin) :
master(masterin), grid(gridin), fields(fieldsin), radiation(radiationin), chemistry(chemistryin)
{
    sw_deposition = inputin.get_item<bool>("deposition", "swdeposition", "", false);
    use_depac = inputin.get_item<bool>("deposition", "use_depac", "", true);
    iratns = inputin.get_item<int>("deposition", "iratns", "");
    hlaw = inputin.get_item<TF>("deposition", "hlaw", "");
    react = inputin.get_item<TF>("deposition", "react", "");
    c_ave_prev_nh3 = inputin.get_item<TF>("deposition", "c_ave_prev_nh3", "");
    pressure = inputin.get_item<TF>("thermo", "pbot", "");
    sw_override_ccomp = inputin.get_item<bool>("deposition", "sw_override_ccomp", "", false);
    ccomp_override_value = inputin.get_item<TF>("deposition", "ccomp_override_value", "", TF(0.0));
    nwet_veg = inputin.get_item<int>("deposition", "nwet_veg", "");
    nwet_soil = inputin.get_item<int>("deposition", "nwet_soil", "");
    nwet_wet = inputin.get_item<int>("deposition", "nwet_wet", "");
    sw_sinphi_prescr = inputin.get_item<bool>("deposition", "sw_sinphi_prescr", "", true);
    if (sw_sinphi_prescr)
    {
        Netcdf_file input_nc(master, "plume_chem_input.nc", Netcdf_mode::Read);
        t_sunrise = input_nc.get_variable<TF>("t_sunrise");
        t_sunset = input_nc.get_variable<TF>("t_sunset");
    } 
}

template <typename TF>
Deposition<TF>::~Deposition()
{
}

template <typename TF>
void Deposition<TF>::init(Input& inputin)
{
    // Read the default deposition velocities. They are needed by 
    // chemistry, even if deposition is disabled.
    vd_nh3  = inputin.get_item<TF>("deposition", "vdnh3", "", (TF)0.0); 

    if (!sw_deposition)
        return;

    if (use_depac)
    {
        // Read per-cell DEPAC land use map from binary file written by plume_chem_input.py.
        // lu_map.0000000 contains one integer per cell: 1=grass, 4=coniferous, 5=deciduous etc.
        auto& gd = grid.get_grid_data();
        lu_map.resize(gd.ijcells);
        std::string lu_map_file = "lu_map.0000000";
        std::ifstream file(lu_map_file, std::ios::binary);
        if (!file)
            throw std::runtime_error("Cannot open " + lu_map_file);
        std::vector<double> buf(gd.imax * gd.jmax);
        file.read(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(double));
        file.close();
        for (int j = gd.jstart; j < gd.jend; ++j)
            for (int i = gd.istart; i < gd.iend; ++i)
            {
                const int ij     = i + j * gd.icells;
                const int ij_buf = (i - gd.istart) + (j - gd.jstart) * gd.imax;
                lu_map[ij] = static_cast<int>(buf[ij_buf]);
            }
    }

    auto& gd = grid.get_grid_data();

    for (auto& name : deposition_tile_names)
        deposition_tiles.emplace(name, Deposition_tile<TF>{});

    for (auto& tile : deposition_tiles)
    {
        tile.second.vdnh3.resize(gd.ijcells);
        tile.second.ra.resize(gd.ijcells);
        tile.second.rb.resize(gd.ijcells);
        tile.second.obuk.resize(gd.ijcells);
        tile.second.ustar.resize(gd.ijcells);
        tile.second.T_surface.resize(gd.ijcells);
        tile.second.rh_surface.resize(gd.ijcells);

	    if (use_depac)
        {
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
    ra_mean.resize(gd.ijcells);
    rb_mean.resize(gd.ijcells);
    std::fill(rb_mean.begin(), rb_mean.end(), TF(0.0));
    T_surface_mean.resize(gd.ijcells);
    rh_surface_mean.resize(gd.ijcells);

    if (use_depac)
    {
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

    deposition_tiles.at("veg" ).long_name = "vegetation";
    deposition_tiles.at("soil").long_name = "bare soil";
    deposition_tiles.at("wet" ).long_name = "wet skin";
    deposition_var = inputin.get_item<TF>("deposition", "deposition_var","", (TF)1e5);
    henry_so2 = inputin.get_item<TF>("deposition", "henry_so2", "", (TF)1e5);
    rsoil_so2 = inputin.get_item<TF>("deposition", "rsoil_so2", "", (TF)250.0);
    rwat_so2 = inputin.get_item<TF>("deposition", "rwat_so2", "", (TF)1.0);
    rws_so2 = inputin.get_item<TF>("deposition", "rws_so2", "", (TF)100.0);
    rmes     = {(TF)0.0};
    rsoil    = {(TF)100.0};
    rcut     = {(TF)1e5};
    rws      = {(TF)10.0};
    rwat     = {(TF)10.0};
    diff     = {(TF)0.13};
    diff_scl = {(TF)1.0};

    // Change diff_scl to diff_scl^(2/3) for use in rb calculation
    for (int i=0; i<1; i++) diff_scl[i] = pow(diff_scl[i], (TF)2.0/(TF)3.0);  // Modified for NH3 only

    for (auto& tile : deposition_tiles)
    {
        std::fill(tile.second.vdnh3.begin(),tile.second.vdnh3.end(), vd_nh3); 
    }
}

template <typename TF>
void Deposition<TF>::create(Stats<TF>& stats, Cross<TF>& cross)
{
    if (!sw_deposition)
        return;
    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = 
        {
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
            "T_surface",
            "T_surface_veg",
            "T_surface_soil",
            "T_surface_wet",
            "rh_surface",
            "rh_surface_veg",
            "rh_surface_soil",
            "rh_surface_wet"
        };
        cross_list = cross.get_enabled_variables(allowed_crossvars);
    }
}

template <typename TF>
void Deposition<TF>::update_time_dependent
    (
    Timeloop<TF>& timeloop,
    Boundary<TF>& boundary,
    Thermo<TF>& thermo,
    TF* restrict vdnh3 
    )
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();
    const std::vector<TF>& rho = thermo.get_basestate_vector("rho");

    std::vector<TF> T_a(gd.ijcells);
    std::vector<TF> RH_a(gd.ijcells);

    // Only retrieve DEPAC-specific values if using DEPAC
    if (use_depac)
    {
        day_of_year = int(timeloop.calc_day_of_year());
        lat = gd.lat;
        if (sw_sinphi_prescr)
        {
            const TF hour = timeloop.calc_hour_of_day();
            const TF pi = 3.14159265358979323846;
            const TF dlen = t_sunset - t_sunrise;
            sinphi = std::sin(pi * (hour - t_sunrise) / dlen);
        }
        else
        {
            const int year = timeloop.get_year();
            const TF secs = TF(timeloop.calc_hour_of_day() * 3600);
            TF azimuth;
            std::tie(sinphi, azimuth) = Radiation_rrtmgp_functions::calc_cos_zenith_angle(lat, gd.lon, day_of_year, secs, year);
        }

        const Radiation_prescribed<TF>& radiation_prescribed = static_cast<const Radiation_prescribed<TF>&>(radiation);
        glrad = radiation_prescribed.get_sw_flux_dn()[0];
        auto tmp2 = fields.get_tmp();
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
                RH_a[ij] = tmp1->fld[ijk] * 100.0;
            }
        fields.release_tmp(tmp1);
        temperature = T_a[gd.istart + gd.jstart * gd.icells];
        rh = RH_a[gd.istart + gd.jstart * gd.icells];
        
        for (const auto& tile_name : deposition_tile_names)
        {
            if (deposition_tiles.count(tile_name) > 0)
            {
                auto& dep_tile = deposition_tiles.at(tile_name);
                std::copy(T_a.begin(), T_a.end(), dep_tile.T_surface.begin());
                std::copy(RH_a.begin(), RH_a.end(), dep_tile.rh_surface.begin());
            }
        }
        
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
    
    auto& tiles = boundary.get_tiles();
    auto& lai = boundary.get_lai();
    auto& water_mask = boundary.get_water_mask();
    auto& c_veg = boundary.get_c_veg();
    TF xmnh3 = 17.031;
    TF xmair = 28.9647;
    TF xmair_i = TF(1) / xmair;
    TF c_ug = TF(1.0e9) * fields.rhoref[gd.kstart] * xmnh3 * xmair_i;
    
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
        calc_deposition_per_tile<TF>
            (
            master,
            tile.first,
            deposition_tiles.at(tile.first).vdnh3.data(),
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
            use_depac,
            glrad,
            sinphi,
            temperature,
            rh,
            T_a.data(),
            RH_a.data(),
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
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells,
            lu_map.data(),
            gd.kstart,
            gd.ijcells
            );
    }
    get_tiled_mean(vdnh3,"nh3",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());

    // Only calculate DEPAC-specific means if using DEPAC
    if (use_depac)
    {
        get_tiled_mean(ra_mean.data(), "ra", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(rb_mean.data(), "rb", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(ccomp_mean.data(), "ccomp_tot", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(cw_mean.data(), "cw", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(cstom_mean.data(), "cstom", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(csoil_eff_mean.data(), "csoil_eff", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(cw_out_mean.data(), "cw_out", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(cstom_out_mean.data(), "cstom_out", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(csoil_out_mean.data(), "csoil_out", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(rc_tot_mean.data(), "rc_tot", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
        get_tiled_mean(rc_eff_mean.data(), "rc_eff", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    }
    
    // cmk: we use the wet-tile info for u* and ra, since these are calculated in lsm with f_wet = 100%
    update_vd_water(vdnh3,"nh3",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
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
        else if (use_depac)
        {
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
    if (name == "nh3")
        return vd_nh3;
    else
    {
        std::string error = "Deposition::get_vd() can't return \"" + name + "\"";
        throw std::runtime_error(error);
    }
}

template<typename TF>
void Deposition<TF>::get_tiled_mean
    (
    TF* restrict fld_out, std::string name, const TF fac,
    const TF* const restrict fveg,
    const TF* const restrict fsoil,
    const TF* const restrict fwet
    )
{
    auto& gd = grid.get_grid_data();
    TF* fld_veg;
    TF* fld_soil;
    TF* fld_wet;
    if (name == "nh3") 
    {
        fld_veg  = deposition_tiles.at("veg").vdnh3.data();
        fld_soil = deposition_tiles.at("soil").vdnh3.data();
        fld_wet  = deposition_tiles.at("wet").vdnh3.data();
    }
    else if (name == "ra")
    {
        fld_veg  = deposition_tiles.at("veg").ra.data();
        fld_soil = deposition_tiles.at("soil").ra.data();
        fld_wet  = deposition_tiles.at("wet").ra.data();
    }
    else if (name == "rb") 
    {
        fld_veg  = deposition_tiles.at("veg").rb.data();
        fld_soil = deposition_tiles.at("soil").rb.data();
        fld_wet  = deposition_tiles.at("wet").rb.data();
    }
    else if (name == "ccomp_tot")
    {
        fld_veg  = deposition_tiles.at("veg").ccomp_tot.data();
        fld_soil = deposition_tiles.at("soil").ccomp_tot.data();
        fld_wet  = deposition_tiles.at("wet").ccomp_tot.data();
    }
    else if (name == "cw")
    {
        fld_veg  = deposition_tiles.at("veg").cw.data();
        fld_soil = deposition_tiles.at("soil").cw.data();
        fld_wet  = deposition_tiles.at("wet").cw.data();
    }
    else if (name == "cstom")
    {
        fld_veg  = deposition_tiles.at("veg").cstom.data();
        fld_soil = deposition_tiles.at("soil").cstom.data();
        fld_wet  = deposition_tiles.at("wet").cstom.data();
    }
    else if (name == "csoil_eff")
    {
        fld_veg  = deposition_tiles.at("veg").csoil_eff.data();
        fld_soil = deposition_tiles.at("soil").csoil_eff.data();
        fld_wet  = deposition_tiles.at("wet").csoil_eff.data();
    }
    else if (name == "cw_out")
    {
        fld_veg  = deposition_tiles.at("veg").cw_out.data();
        fld_soil = deposition_tiles.at("soil").cw_out.data();
        fld_wet  = deposition_tiles.at("wet").cw_out.data();
    }
    else if (name == "cstom_out")
    {
        fld_veg  = deposition_tiles.at("veg").cstom_out.data();
        fld_soil = deposition_tiles.at("soil").cstom_out.data();
        fld_wet  = deposition_tiles.at("wet").cstom_out.data();
    }
    else if (name == "csoil_out")
    {
        fld_veg  = deposition_tiles.at("veg").csoil_out.data();
        fld_soil = deposition_tiles.at("soil").csoil_out.data();
        fld_wet  = deposition_tiles.at("wet").csoil_out.data();
    }
    else if (name == "rc_tot")
    {
        fld_veg  = deposition_tiles.at("veg").rc_tot.data();
        fld_soil = deposition_tiles.at("soil").rc_tot.data();
        fld_wet  = deposition_tiles.at("wet").rc_tot.data();
    }
    else if (name == "rc_eff")
    {
        fld_veg  = deposition_tiles.at("veg").rc_eff.data();
        fld_soil = deposition_tiles.at("soil").rc_eff.data();
        fld_wet  = deposition_tiles.at("wet").rc_eff.data();
    }
    else if (name == "T_surface")
    { 
        fld_veg  = deposition_tiles.at("veg").T_surface.data();
        fld_soil = deposition_tiles.at("soil").T_surface.data();
        fld_wet  = deposition_tiles.at("wet").T_surface.data();
    }
    else if (name == "rh_surface")
    {
        fld_veg  = deposition_tiles.at("veg").rh_surface.data();
        fld_soil = deposition_tiles.at("soil").rh_surface.data();
        fld_wet  = deposition_tiles.at("wet").rh_surface.data();
    }
    else
        throw std::runtime_error("Cannot calculate tiled mean for variable \"" + name + "\"\\n");

    calc_tiled_mean
        (
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
        gd.icells
        );
}

template<typename TF>
void Deposition<TF>::update_vd_water
    (
    TF* restrict fld_out, std::string name,
    const TF* const restrict ra,
    const TF* const restrict ustar,
    const int* const restrict water_mask,
    const TF* const restrict diff_scl,
    const TF* const restrict rwat
    )
{
    auto& gd = grid.get_grid_data();
    TF diff_scl_val;
    TF rwat_val;
    if (name == "nh3") 
    {
        diff_scl_val = diff_scl[0];  // Using first index for NH3
        rwat_val = rwat[0];  // Using first index for NH3
    }
    else
        throw std::runtime_error("Cannot update vd to water for variable \"" + name + "\"\\n");

    calc_vd_water
        (
        fld_out,
        ra,
        ustar,
        water_mask,
        diff_scl_val,
        rwat_val,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells
        );
}

template class Deposition<double>;
