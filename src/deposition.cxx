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
 * c_ave_prev_so2        = float   : Long-term mean SO2 concentration [ug m-3] (default: 0). 0 disables the NH3 co-deposition factor   // ADDED   // CHANGED: optional
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
 * vdnh3_veg, vdnh3_soil, vdnh3_wet     : NH3 exchange velocity v_e [m s-1] (equals deposition velocity when compensation points are absent)
 * 
 * # Resistances (DEPAC only):
 * ra, rb                               : Aerodynamic and quasi-laminar resistances [s m-1]
 * ra_veg, ra_soil, ra_wet              : Per-tile aerodynamic resistance [s m-1]
 * rb_veg, rb_soil, rb_wet              : Per-tile quasi-laminar resistance [s m-1]
 * rc_tot, rc_tot_veg, rc_tot_soil, rc_tot_wet : Total canopy resistance [s m-1]
 * rc_eff, rc_eff_veg, rc_eff_soil, rc_eff_wet : Effective canopy resistance [s m-1]
 * cw, cstom, csoil_eff                 : External leaf, stomatal, effective soil pathway resistances R_w, R_stom, R_soil,eff [s m-1]
 * cw_out, cstom_out, csoil_out         : External leaf, stomatal, soil compensation points chi_w, chi_stom, chi_soil,eff [ug m-3]
 * ccomp_tot, ccomp_tot_veg, ccomp_tot_soil, ccomp_tot_wet : Total compensation point [ug m-3]
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
#include "field3d_io.h"   // ADDED: MPI-aware 2D binary loader (same one the LSM uses for c_veg, lai, ...)
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
#include <cstdlib>   // DEBUG: getenv
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdexcept>      // ADDED: std::runtime_error is thrown below
#include <utility>
#include "monin_obukhov.h"


// C linkage for DEPAC Fortran wrapper (iso_c_binding, Sect. 3)
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
        int use_input_ccomp,   // CHANGED: was bool (1 byte). The Fortran wrapper reads integer(c_int) (4 bytes), so the upper bytes could be garbage
        float tsea,            // ADDED: water surface temperature [K]. DEPAC reads it for water land-use classes (le_drydepos_gas_depac.f90 l. 1755)   // CHANGED: comma added
        float c_ave_prev_so2   // ADDED: long-term mean SO2 [ug m-3]. Last argument, same position as in depac_wrapper_c (wrapper_depac.f90)
        );
}

namespace
{
    // =====================================================================================
    // ADDED: land-use handling. This file contains NO land-use class numbers any more.
    //
    // plume_chem_input.py is the single place where land use is defined (table ifs_to_depac and table depac_is_forest).
    // It writes two 2D binaries in the model precision, like the other LSM binaries:
    //   lu_map.0000000     DEPAC land-use class 1..18 of every column (depac_lu.inc)
    //   sai_add.0000000    SAI - LAI of every column: 1 for forest, 0 otherwise
    // The kernels use lu_map as the DEPAC class and SAI = LAI + sai_add.
    // REDUNDANT (removed): the hard-coded test (local_lu == 4 || local_lu == 5 || local_lu == 12).
    // =====================================================================================
    constexpr int n_depac_lu = 18;   // number of DEPAC land-use classes in depac_lu.inc (1-indexed)

    // Convert one value read from lu_map to a valid DEPAC class. Throws if it is not an integer in 1..n_depac_lu.
    template<typename TF>
    int checked_depac_lu(const TF value)
    {
        const int lu = static_cast<int>(std::lround(value));
        if (std::abs(value - TF(lu)) > TF(1e-6) || lu < 1 || lu > n_depac_lu)
            throw std::runtime_error("Deposition: lu_map holds " + std::to_string(value) + ", which is not a DEPAC land-use class 1.." + std::to_string(n_depac_lu) + " (wrong file, wrong precision, or the IFS to DEPAC table was not applied in plume_chem_input.py)");
        return lu;
    }

    // Convert one value read from sai_add to 0 or 1. Throws otherwise.
    template<typename TF>
    TF checked_sai_add(const TF value)
    {
        if (std::abs(value) > TF(1e-6) && std::abs(value - TF(1.0)) > TF(1e-6))
            throw std::runtime_error("Deposition: sai_add holds " + std::to_string(value) + ", expected 0 or 1");
        return (value > TF(0.5)) ? TF(1.0) : TF(0.0);
    }

    // DEBUG (temporary, remove after the crash is found): with the environment variable DEPAC_DEBUG set, every call of
    // depac_wrapper prints its inputs to stderr BEFORE the call, so the last line of a rank before a segfault is the crashing call.
    template<typename TF>
    void depac_debug_print(const char* tile, const int i, const int j, const int lu, const TF lai, const TF sai, const TF c_veg, const TF fraction,
                           const int nwet, const TF ustar, const TF ra, const TF rb, const TF T_c, const TF rh, const TF nh3)
    {
        static const bool on = (std::getenv("DEPAC_DEBUG") != nullptr);
        if (!on) return;
        const char* r = std::getenv("PMIX_RANK");
        std::fprintf(stderr, "DEPAC_DEBUG rank=%s tile=%s i=%d j=%d lu=%d lai=%g sai=%g c_veg=%g frac=%g nwet=%d ustar=%g ra=%g rb=%g T=%g rh=%g nh3=%g\n",
                     r ? r : "?", tile, i, j, lu, double(lai), double(sai), double(c_veg), double(fraction), nwet, double(ustar), double(ra), double(rb), double(T_c), double(rh), double(nh3));
        std::fflush(stderr);
    }

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
        // REDUNDANT (removed): const TF sai,   // SAI is now lai + sai_add per cell
        const TF lat,
        const int day_of_year,
        // REDUNDANT (removed): const int nwet,   // was never assigned and never used in the kernels
        const int nwet_veg,
        const int nwet_soil,
        const int nwet_wet,
        // REDUNDANT (removed): const int lu,   // land-use class is now taken per cell from lu_map
        const int iratns,
        const TF hlaw,
        const TF react,
        const TF c_ave_prev_nh3,
        const TF c_ave_prev_so2,   // ADDED: long-term mean SO2 [ug m-3]
        const TF catm,
        const TF pressure,
        const bool sw_override_ccomp,
        const TF ccomp_override_value,
        Deposition_tile_map<TF>& deposition_tiles,
        const Chemistry<TF>& chemistry,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int jj,
        const int* const restrict lu_map,
        const TF* const restrict sai_add,   // ADDED: SAI - LAI per cell (0 or 1), from sai_add.0000000
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

                    // CHANGED: lu_map is the DEPAC class and sai_add is 0 or 1 (forest), both from plume_chem_input.py.
                    // REDUNDANT (removed): const bool is_forest = (local_lu == 4 || local_lu == 5 || local_lu == 12);
                    const int local_lu = lu_map[ij];
                    const TF local_sai = lai[ij] + sai_add[ij];

                    // Quasi-laminar resistance R_b for vegetated tile: beta=2
                    const TF rb = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[0];
                    deposition_tiles.at(lu_type).rb.data()[ij] = rb;

                    // Retrieve NH3 mole fraction c_target at reference height z_ref from Chemistry module
                    const TF nh3_conc_value = chemistry.get_c_target()[ij];
                    const TF nh3_ugm3 = nh3_conc_value * chemistry.get_rho_target()[ij] * TF(1e9) * xmnh3 * xmair_i;  // convert mol mol-1 to ug m-3 for depac (chi_a, table 2)
                    // Conductance/resistance variables

                    float rc_tot;           // total canopy resistance R_c,tot [s m-1]
                    float ccomp_tot = 0.0;  // total canopy compensation point chi_c [ug m-3]
                    float rc_eff;           // effective canopy resistance R_c,eff [s m-1]
                    float gsoil_eff_out;    // effective soil conductance G_soil,eff [m s-1]
                    float rsoil_eff_out;    // effective soil resistance R_soil,eff [s m-1]
                    float gw_out;           // external leaf surface conductance G_w [m s-1]
                    float gstom_out;        // stomatal conductance G_stom [m s-1]
                    
                    // Compensation points returned by DEPAC
                    // Note: tile fields cw/cstom/csoil_eff store the corresponding resistances (1/G); cw_out/cstom_out/csoil_out store the compensation points
                    float cw_out;           // external leaf surface compensation point chi_w [ug m-3]
                    float cstom_out;        // stomatal compensation point chi_stom [ug m-3]
                    float csoil_out;        // soil compensation point chi_soil,eff [ug m-3]
                    int status;             // DEPAC return status (0 = success)
                    // Initialize compensation point: use override value if sw_override_ccomp, otherwise zero
                    bool use_input_ccomp = false;
                    if (sw_override_ccomp)
                    {
                        ccomp_tot = ccomp_override_value;
                        use_input_ccomp = true;
                    }

                    depac_debug_print<TF>("veg", i, j, local_lu, lai[ij], local_sai, c_veg[ij], fraction[ij], nwet_veg, ustar[ij], ra[ij], rb, T_a[ij]-TF(273.15), RH_a[ij], nh3_ugm3);   // DEBUG
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
                        use_input_ccomp,
                        T_a[ij],  // ADDED: tsea [K]. CHECK: air temperature at the first model level is used as water surface temperature   // CHANGED: comma added
                        c_ave_prev_so2   // ADDED: long-term mean SO2 [ug m-3], last argument
                        );

                    // ADDED: a non-zero DEPAC status was silently ignored before, and the tile then kept its old v_e
                    if (status != STATUS_OK)
                        throw std::runtime_error("Deposition: DEPAC returned status " + std::to_string(status) + " for tile " + lu_type + " at i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", lu=" + std::to_string(local_lu));
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
                        // ADDED: rc_tot = -9999 means DEPAC found no deposition path (le_drydepos_gas_depac.f90 l. 165-167, 1631-1632).
                        // Then rc_eff = -9999 as well, and 1/(ra + rb - 9999) < 0 would be a spurious emission.
                        if (rc_tot < TF(-9998))
                        {
                            vdnh3[ij] = TF(0);
                        }
                        else if (std::abs(total_resistance) > (TF)1e-6)
                        {
                            vdnh3[ij] = (TF)1.0 / total_resistance;
                            // v_e < 0 when R_c,eff < 0, indicating net NH3 emission
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
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + kstart*ijcells;  

                    if (fraction[ij] < (TF)1e-12)
                        continue;

                    // CHANGED: the soil tile used the members "lu" and "sai", which were never assigned (undefined values).
                    // It now uses the land-use class and SAI of the cell itself.
                    // CHECK: modelling choice. The soil path of DEPAC depends on the class (rsoil, rinc via b(lu)*h(lu)*SAI/u*).
                    const int local_lu = lu_map[ij];
                    const TF local_sai = lai[ij] + sai_add[ij];

                    const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[0];
                    deposition_tiles.at(lu_type).rb.data()[ij] = rb;
                    const TF nh3_conc_value = chemistry.get_c_target()[ij];
                    const TF nh3_ugm3 = nh3_conc_value * chemistry.get_rho_target()[ij] * TF(1e9) * xmnh3 * xmair_i;
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

                    depac_debug_print<TF>("soil", i, j, local_lu, lai[ij], local_sai, c_veg[ij], fraction[ij], nwet_soil, ustar[ij], ra[ij], rb, T_a[ij]-TF(273.15), RH_a[ij], nh3_ugm3);   // DEBUG
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
                        local_sai,   // CHANGED: was the unassigned member "sai"
                        nwet_soil,
                        local_lu,    // CHANGED: was the unassigned member "lu"
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
                        use_input_ccomp,
                        T_a[ij],  // ADDED: tsea [K]. CHECK: air temperature at the first model level is used as water surface temperature   // CHANGED: comma added
                        c_ave_prev_so2   // ADDED: long-term mean SO2 [ug m-3], last argument
                        );

                    // ADDED: a non-zero DEPAC status was silently ignored before, and the tile then kept its old v_e
                    if (status != STATUS_OK)
                        throw std::runtime_error("Deposition: DEPAC returned status " + std::to_string(status) + " for tile " + lu_type + " at i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", lu=" + std::to_string(local_lu));
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
                        // ADDED: rsoil_eff_out = -9999 means no soil path (depac_wrapper sets it when gsoil_eff_out <= 0,
                        // e.g. grass and "other", where DEPAC has R_inc = -999, le_drydepos_gas_depac.f90 l. 1605-1606).
                        // Then v_e = 0, not 1/(ra + rb - 9999) < 0 (a spurious emission).
                        if (rsoil_eff_out < TF(0))
                        {
                            vdnh3[ij] = TF(0);
                        }
                        else if (std::abs(total_resistance) > (TF)1e-6)
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
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    const int ijk = i + j*jj + kstart*ijcells; 
                    const TF nh3_conc_value = chemistry.get_c_target()[ij];
                    const TF nh3_ugm3 = nh3_conc_value * chemistry.get_rho_target()[ij] * TF(1e9) * xmnh3 * xmair_i;
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
                        // CHANGED: same as the veg tile (REDUNDANT: hard-coded is_forest test removed).
                        const int local_lu = lu_map[ij];
                        const TF local_sai = lai[ij] + sai_add[ij];
                        const TF rb = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[0];
                        deposition_tiles.at(lu_type).rb.data()[ij] = rb;

                        depac_debug_print<TF>("wet_veg", i, j, local_lu, lai[ij], local_sai, c_veg[ij], fraction[ij], nwet_wet, ustar[ij], ra[ij], rb, T_a[ij]-TF(273.15), RH_a[ij], nh3_ugm3);   // DEBUG
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
                            use_input_ccomp,
                            T_a[ij],  // ADDED: tsea [K]. CHECK: air temperature at the first model level is used as water surface temperature   // CHANGED: comma added
                            c_ave_prev_so2   // ADDED: long-term mean SO2 [ug m-3], last argument
                            );

                        // ADDED: a non-zero DEPAC status was silently ignored before, and the tile then kept its old v_e
                        if (status != STATUS_OK)
                            throw std::runtime_error("Deposition: DEPAC returned status " + std::to_string(status) + " for tile " + lu_type + " at i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", lu=" + std::to_string(local_lu));
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
                            // ADDED: rc_tot = -9999 means DEPAC found no deposition path (le_drydepos_gas_depac.f90 l. 165-167, 1631-1632).
                            // Then rc_eff = -9999 as well, and 1/(ra + rb - 9999) < 0 would be a spurious emission.
                            if (rc_tot < TF(-9998))
                            {
                                vdnh3[ij] = TF(0);
                            }
                            else if (std::abs(total_resistance) > (TF)1e-6)
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
                        // CHANGED: this branch used the unassigned members "lu" and "sai". It now uses the cell's own class and SAI.
                        // CHECK: modelling choice, see the soil tile above.
                        const int local_lu = lu_map[ij];
                        const TF local_sai = lai[ij] + sai_add[ij];

                        const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[0];
                        deposition_tiles.at(lu_type).rb.data()[ij] = rb;

                        depac_debug_print<TF>("wet_noveg", i, j, local_lu, lai[ij], local_sai, c_veg[ij], fraction[ij], nwet_wet, ustar[ij], ra[ij], rb, T_a[ij]-TF(273.15), RH_a[ij], nh3_ugm3);   // DEBUG
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
                            local_sai,   // CHANGED: was the unassigned member "sai"
                            nwet_wet,
                            local_lu,    // CHANGED: was the unassigned member "lu"
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
                            use_input_ccomp,
                            T_a[ij],  // ADDED: tsea [K]. CHECK: air temperature at the first model level is used as water surface temperature   // CHANGED: comma added
                            c_ave_prev_so2   // ADDED: long-term mean SO2 [ug m-3], last argument
                            );

                        // ADDED: a non-zero DEPAC status was silently ignored before, and the tile then kept its old v_e
                        if (status != STATUS_OK)
                            throw std::runtime_error("Deposition: DEPAC returned status " + std::to_string(status) + " for tile " + lu_type + " at i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", lu=" + std::to_string(local_lu));
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
	                        // ADDED: rsoil_eff_out = -9999 means no soil path (depac_wrapper sets it when gsoil_eff_out <= 0,
	                        // e.g. grass and "other", where DEPAC has R_inc = -999, le_drydepos_gas_depac.f90 l. 1605-1606).
	                        // Then v_e = 0, not 1/(ra + rb - 9999) < 0 (a spurious emission).
	                        if (rsoil_eff_out < TF(0))
	                        {
	                            vdnh3[ij] = TF(0);
	                        }
	                        else if (std::abs(total_resistance) > (TF)1e-6)
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
        // REDUNDANT (removed): const TF sai,   // SAI is now lai + sai_add per cell
        const TF lat,
        const int day_of_year,
        // REDUNDANT (removed): const int nwet,   // was never assigned and never used in the kernels
        const int nwet_veg,      
        const int nwet_soil,
        const int nwet_wet,
        // REDUNDANT (removed): const int lu,   // land-use class is now taken per cell from lu_map
        const int iratns,
        const TF hlaw,
        const TF react,
        const TF c_ave_prev_nh3,
        const TF c_ave_prev_so2,   // ADDED: long-term mean SO2 [ug m-3]
        const TF catm,
        const TF pressure,
        const bool sw_override_ccomp,
        const TF ccomp_override_value,
        Deposition_tile_map<TF>& deposition_tiles,
        const Chemistry<TF>& chemistry,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int jj,
        const int* const restrict lu_map,
        const TF* const restrict sai_add,   // ADDED: SAI - LAI per cell (0 or 1), from sai_add.0000000
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
                // REDUNDANT (removed): sai,   // argument removed together with the parameter
                lat,
                day_of_year,
                // REDUNDANT (removed): nwet,   // argument removed together with the parameter
                nwet_veg, 
                nwet_soil,
                nwet_wet,
                // REDUNDANT (removed): lu,   // argument removed together with the parameter
                iratns,
                hlaw,
                react,
                c_ave_prev_nh3,
                c_ave_prev_so2,   // ADDED
                catm,
                pressure,
                sw_override_ccomp,
                ccomp_override_value,
                deposition_tiles,
                chemistry,
                istart, iend,
                jstart, jend,
                jj,
                lu_map,
                sai_add,   // ADDED
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
    c_ave_prev_so2 = inputin.get_item<TF>("deposition", "c_ave_prev_so2", "", TF(0.));   // CHANGED: optional, 0 = no NH3 co-deposition factor
    if (use_depac && c_ave_prev_so2 <= TF(0.))   // ADDED: make the "off" state visible in the log, since a missing key is silent otherwise
        master.print_message("Deposition: c_ave_prev_so2 <= 0, the NH3 co-deposition factor in DEPAC is not used\n");
    pressure = inputin.get_item<TF>("thermo", "pbot", "");
    sw_override_ccomp = inputin.get_item<bool>("deposition", "sw_override_ccomp", "", false);
    ccomp_override_value = inputin.get_item<TF>("deposition", "ccomp_override_value", "", TF(0.0));
    nwet_veg = inputin.get_item<int>("deposition", "nwet_veg", "");
    nwet_soil = inputin.get_item<int>("deposition", "nwet_soil", "");
    nwet_wet = inputin.get_item<int>("deposition", "nwet_wet", "");
    catm = TF(0.0);   // ADDED: was passed on uninitialised. Not used by the DEPAC kernels (nh3_ugm3 is computed from Chemistry).
    sw_sinphi_prescr = inputin.get_item<bool>("deposition", "sw_sinphi_prescr", "", true);
    // REDUNDANT (moved to Deposition::create): the read of t_sunrise and t_sunset. Reason: the file name "plume_chem_input.nc" was hard-coded.
    // if (sw_sinphi_prescr)
    // {
    //     Netcdf_file input_nc(master, "plume_chem_input.nc", Netcdf_mode::Read);
    //     t_sunrise = input_nc.get_variable<TF>("t_sunrise");
    //     t_sunset = input_nc.get_variable<TF>("t_sunset");
    // }
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

    // REDUNDANT (moved to Deposition::create): the lu_map reader that was here.
    // Reason: it read imax*jmax values from the start of the global file on every MPI rank, so with npx*npy > 1
    // every rank got the wrong sub-domain. The old code, for reference:
    //     std::vector<double> buf(gd.imax * gd.jmax);
    //     file.read(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(double));
    //     lu_map[ij] = static_cast<int>(buf[ij_buf]);

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
    for (int i=0; i<1; i++) diff_scl[i] = pow(diff_scl[i], (TF)2.0/(TF)3.0);  // (Sc/Pr)^(2/3) scaling for R_b; set to 1.0 for NH3

    for (auto& tile : deposition_tiles)
    {
        std::fill(tile.second.vdnh3.begin(),tile.second.vdnh3.end(), vd_nh3); 
    }
}

template <typename TF>
void Deposition<TF>::create(Stats<TF>& stats, Cross<TF>& cross, Netcdf_handle& input_nc)   // CHANGED: added input_nc (<sim_name>_input.nc)
{
    if (!sw_deposition)
        return;

    // ADDED: read the sunrise and sunset times from <sim_name>_input.nc (moved here from the constructor)
    if (use_depac && sw_sinphi_prescr)
    {
        t_sunrise = input_nc.get_variable<TF>("t_sunrise");
        t_sunset  = input_nc.get_variable<TF>("t_sunset");
    }

    // ADDED: load lu_map (DEPAC class) and sai_add (0/1) from plume_chem_input.py. Done here, because it needs the grid and the
    // fields, as the LSM does for its own binaries.
    if (use_depac)
    {
        auto& gd = grid.get_grid_data();

        // Both files are written in the model precision (TF) like the other LSM binaries.
        // Field3d_io::load_xy_slice hands every MPI rank its own sub-domain.
        // CHECK: signature copied from how MicroHH's LSM loads its 2D binaries
        //        (field3d_io.load_xy_slice(field, tmp->fld.data(), filename), returns non-zero on failure).
        //        I could not compile against your MicroHH tree. If the compiler rejects this call, compare with load_2d_field() in boundary_surface_lsm.cxx.
        Field3d_io<TF> field3d_io(master, grid);
        auto tmp = fields.get_tmp();

        auto load_2d = [&](const std::string& filename, std::vector<TF>& field)
        {
            field.assign(gd.ijcells, TF(0.0));
            master.print_message("Loading \"%s\" ... ", filename.c_str());
            if (field3d_io.load_xy_slice(field.data(), tmp->fld.data(), filename.c_str()))
            {
                master.print_message("FAILED\n");
                fields.release_tmp(tmp);
                throw std::runtime_error("Deposition: cannot load " + filename);
            }
            master.print_message("OK\n");
        };

        std::vector<TF> lu_in, sai_add_in;
        load_2d("lu_map.0000000", lu_in);
        load_2d("sai_add.0000000", sai_add_in);
        fields.release_tmp(tmp);

        // Check every interior cell. Only interior cells are used by the kernels.
        lu_map.assign(gd.ijcells, 0);
        sai_add.assign(gd.ijcells, TF(0.0));
        for (int j = gd.jstart; j < gd.jend; ++j)
            for (int i = gd.istart; i < gd.iend; ++i)
            {
                const int ij = i + j * gd.icells;
                lu_map[ij]  = checked_depac_lu<TF>(lu_in[ij]);
                sai_add[ij] = checked_sai_add<TF>(sai_add_in[ij]);
            }
    }

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

    // Adjust ra from zsl to z_target using stability function ratio
    if (use_depac && chemistry.get_sw_adapt_ref_height())
    {
	namespace most = Monin_obukhov;
        const TF z_sl = gd.z[gd.kstart];
        const std::vector<TF>& z_target = chemistry.get_z_target();
        
        for (const auto& tile_name : deposition_tile_names)
        {
            if (deposition_tiles.count(tile_name) > 0)
            {
                auto& dep_tile = deposition_tiles.at(tile_name);
                
                for (int j = gd.jstart; j < gd.jend; ++j)
                    for (int i = gd.istart; i < gd.iend; ++i)
                    {
                        const int ij = i + j * gd.icells;
                        
                        if (dep_tile.ra[ij] > (TF)0 && z_target[ij] > (TF)0)
                        {
                            const TF L = dep_tile.obuk[ij];
                            // const TF fh_zsl = most::fh(z_sl, boundary.get_tiles().at(tile_name).z0h[ij], L);
                            // const TF fh_ztarget = most::fh(z_target[ij], boundary.get_tiles().at(tile_name).z0h[ij], L);
			    const TF fh_zsl = most::fh(z_sl, boundary.get_z0h()[ij], L);
			    const TF fh_ztarget = most::fh(z_target[ij], boundary.get_z0h()[ij], L);
                            
                            if (std::abs(fh_ztarget) > (TF)1e-15)
                                dep_tile.ra[ij] *= (fh_zsl / fh_ztarget);
                        }
                    }
            }
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
            fields.sp.at("nh3")->fld.data(),  // nh3_concentration: unused in DEPAC path; c_target from Chemistry module is used instead
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
            // REDUNDANT (removed): sai,   // argument removed together with the parameter
            lat,
            day_of_year,
            // REDUNDANT (removed): nwet,   // argument removed together with the parameter
            nwet_veg,
            nwet_soil,
            nwet_wet,
            // REDUNDANT (removed): lu,   // argument removed together with the parameter
            iratns,
            hlaw,
            react,
            c_ave_prev_nh3,
            c_ave_prev_so2,   // ADDED
            catm,
            pressure,
            sw_override_ccomp,
            ccomp_override_value,
            deposition_tiles,
            chemistry,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells,
            lu_map.data(),
            sai_add.data(),   // ADDED
            gd.kstart,
            gd.ijcells
            );
    }
    get_tiled_mean(vdnh3,"nh3",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());

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
    
    // For open water cells, recalculate v_e using wet-tile R_a and u_* (f_wet = 100% in LSM)
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
