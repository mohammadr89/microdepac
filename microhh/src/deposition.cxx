/* 
 * This code has been modified to handle only NH3 deposition calculations.
 * All other chemical species (O3, NO, NO2, HNO3, H2O2, ROOH, HCHO) have been removed.
 * 
 * Original species order was:
 * [0] = O3
 * [1] = NO   
 * [2] = NO2
 * [3] = HNO3
 * [4] = H2O2
 * [5] = ROOH
 * [6] = HCHO
 *
 * Modified version now only contains:
 * [0] = NH3
 */

/*
 *Integration of DEPAC Resistance Calculations in MicroHH
 *
 *1. Problem:
 *   - MicroHH uses IFS land surface scheme with vegetation, soil, and wet surface tiles.
 *   - DEPAC resistance calculations need integration while keeping IFS's Ra and Rb.
 *   - Vegetation types (grass vs forest) were not properly distinguished by LAI.
 *
 *2. Key Changes:
 *   a) LAI-Based Land Use Determination:
 *      - Grass: LAI ≤ 3.5 (land use = 1, SAI = LAI)
 *      - Forest: LAI > 3.5 (land use = 4, SAI = LAI + 1.0)
 *
 *   b) Unit Conversion:
 *      - mol/mol → µg/m³: `conc_ugm3 = conc_molmol * pressure * M_NH3 / (R * T) * 1e6`
 *      - mol/mol → ppb: `conc_ppb = conc_molmol * 1e9`
 *
 *3. Implementation:
 *   - Vegetation (lu_type = "veg"): LAI-based land use and SAI applied in DEPAC.
 *   - Soil (lu_type = "soil"): No changes, fixed soil parameters.
 *   - Wet (lu_type = "wet"): Split into wet vegetation and wet soil.
 *      - Wet vegetation follows the same rules as dry vegetation.
 *
 *4. DEPAC Integration:
 *   - Vegetation presence: `LAI_present = (lai > 0.0)`
 *   - Stomatal conductance: `Gsto = lai * gsto`
 *   - Total resistance: `1/RC = 1/Rstom + 1/Rw + 1/Rsoil_eff`
 *
 *5. Benefits:
 *   - Accurate vegetation parameters based on LAI.
 *   - Differentiates grass vs forest with proper resistance scaling.
 *   - Handles both dry and wet conditions effectively.
 */


//#include <cstdio>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <utility>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "stats.h"
#include "netcdf_interface.h"
#include "constants.h"
#include "timeloop.h"
#include "deposition.h"
#include "boundary_surface_lsm.h"
#include "boundary.h"
#include "cross.h"


// Added: C linkage for DEPAC Fortran wrapper
extern "C" {
    void depac_wrapper(
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
            float p  // Added pressure parameter
                );
}

namespace {

    template<typename TF>
        void calc_tiled_mean(
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
            //TF avg_dep = sum_dep / n_dep;

            //for (int j=jstart; j<jend; ++j)
            //    #pragma ivdep
            //    for (int i=istart; i<iend; ++i)
            //    {
            //        const int ij = i + j*icells;
            //        fld[ij] = avg_dep;
            //    }
        }


    template<typename TF>
        void calc_deposition_per_tile(
                const std::string lu_type,
                TF* restrict vdnh3,              // Output: NH3 deposition velocity
                const TF* const restrict lai,     // Leaf Area Index
                const TF* const restrict c_veg,   // Vegetation fraction
                const TF* const restrict rs,      // Surface resistance
                const TF* const restrict rs_veg,  // Vegetation surface resistance  
                const TF* const restrict ra,      // Aerodynamic resistance (from MicroHH)
                const TF* const restrict ustar,   // Friction velocity
                const TF* const restrict fraction, // Tile fraction
                const TF* const restrict diff_scl, // Diffusion scaling factors
                const TF glrad,          // Global radiation
                const TF sinphi,         // Solar elevation angle
                const TF temperature,    // Temperature [K]
                const TF rh,            // Relative humidity [%]
                const TF sai,           // Stem Area Index
                const TF lat,           // Latitude [degrees]
                const int day_of_year,  // Day of year
                const int nwet,         // Surface wetness 
                const int lu,           // Land use type
                const int iratns,       // NH3 compensation point option
                const TF hlaw,          // Henry's law constant
                const TF react,         // Reactivity factor 
                const TF c_ave_prev_nh3, // Previous NH3 concentration
                const TF catm,          // Atmospheric NH3 concentration
                const TF pressure,      // Added pressure parameter
                const int istart, const int iend,
                const int jstart, const int jend,
                const int jj)
                {
                    const TF ckarman = 0.4;  // von Karman constant
                    const int STATUS_OK = 0;  // Status code for successful DEPAC calls

                    if (lu_type == "veg") {
                        // Vegetation tile handling
                        for (int j=jstart; j<jend; ++j)
                            for (int i=istart; i<iend; ++i) {
                                const int ij = i + j*jj;

                                if (fraction[ij] < (TF)1e-12)
                                    continue;
                                // NEW: Automatic determination of land use type and SAI based on LAI
                                // This allows DEPAC to use different parameters for grass vs forest
                                int local_lu;
                                TF local_sai;
                                if (lai[ij] <= 3.5) {
                                    local_lu = 1;  // grass
                                    local_sai = lai[ij];  // For grass, SAI = LAI
                                } else {
                                    local_lu = 4;  // coniferous forest
                                    local_sai = lai[ij] + 1.0;  // For forest, add stem area
                                }

                                // Keep IFS Ra and use vegetation Rb scaling
                                const TF rb = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[0];

                                // Call DEPAC wrapper for dry vegetation
                                char compnam[4] = "NH3";
                                float rc_tot, ccomp_tot, rc_eff;
                                float gsoil_eff_out, rsoil_eff_out;
                                int status;

                                depac_wrapper(
                                        compnam,
                                        day_of_year,
                                        lat,
                                        temperature,
                                        ustar[ij],
                                        glrad,
                                        sinphi,
                                        rh,
                                        lai[ij],
                                        //sai,
                                        local_sai,        // CHANGED: Use calculated SAI
                                        0,  // nwet = 0 for dry vegetation
                                            //lu,
                                        local_lu,         // CHANGED: Use LAI-determined land use type
                                        iratns,
                                        &rc_tot,
                                        &ccomp_tot,
                                        hlaw,
                                        react,
                                        &status,
                                        c_ave_prev_nh3,
                                        ra[ij],
                                        rb,
                                        catm,
                                        &rc_eff,
                                        &gsoil_eff_out,
                                        &rsoil_eff_out,
                                        pressure
                                            );

                                        // Calculate deposition velocity using resistance analogy
                                        if (status == STATUS_OK) {
                                            vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rc_eff);
                                        }
                            }
                    }
                    else if (lu_type == "soil") {
                        // Bare soil tile handling  
                        for (int j=jstart; j<jend; ++j)
                            for (int i=istart; i<iend; ++i) {
                                const int ij = i + j*jj;

                                if (fraction[ij] < (TF)1e-12)
                                    continue;

                                // Use soil Rb scaling
                                const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[0];

                                // Call DEPAC wrapper for dry soil
                                char compnam[4] = "NH3";
                                float rc_tot, ccomp_tot, rc_eff;
                                float gsoil_eff_out, rsoil_eff_out;
                                int status;

                                depac_wrapper(
                                        compnam,
                                        day_of_year,
                                        lat,
                                        temperature,
                                        ustar[ij],
                                        glrad,
                                        sinphi,
                                        rh,
                                        lai[ij],
                                        sai,
                                        0,  // nwet = 0 for dry soil
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
                                        catm,
                                        &rc_eff,
                                        &gsoil_eff_out,
                                        &rsoil_eff_out,
                                        pressure
                                            );

                                        if (status == STATUS_OK) {
                                            vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rsoil_eff_out);
                                        }
                            }
                    }
                    else if (lu_type == "wet") {
                        // Wet surfaces handling (both vegetation and soil)
                        for (int j=jstart; j<jend; ++j)
                            for (int i=istart; i<iend; ++i) {
                                const int ij = i + j*jj;

                                if (fraction[ij] < (TF)1e-12)
                                    continue;

                                char compnam[4] = "NH3";
                                float rc_tot, ccomp_tot, rc_eff;
                                float gsoil_eff_out, rsoil_eff_out;
                                int status;

                                if (c_veg[ij] > 0) {
                                    // NEW: Added same LAI-based determination for wet vegetation
                                    int local_lu;
                                    TF local_sai;
                                    if (lai[ij] <= 3.5) {
                                        local_lu = 1;  // grass
                                        local_sai = lai[ij];  // For grass, SAI = LAI
                                    } else {
                                        local_lu = 4;  // coniferous forest
                                        local_sai = lai[ij] + 1.0;  // For forest, add stem area
                                    }

                                    // Wet vegetation case
                                    const TF rb = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[0];

                                    depac_wrapper(
                                            compnam,
                                            day_of_year,
                                            lat,
                                            temperature,
                                            ustar[ij],
                                            glrad,
                                            sinphi,
                                            rh,
                                            lai[ij],
                                            //sai,
                                            local_sai,        // CHANGED: Use calculated SAI
                                            1,  // nwet = 1 for wet conditions
                                                //lu,
                                            local_lu,         // CHANGED: Use LAI-determined land use type
                                            iratns,
                                            &rc_tot,
                                            &ccomp_tot,
                                            hlaw,
                                            react,
                                            &status,
                                            c_ave_prev_nh3,
                                            ra[ij],
                                            rb,
                                            catm,
                                            &rc_eff,
                                            &gsoil_eff_out,
                                            &rsoil_eff_out,
                                            pressure
                                                );

                                            if (status == STATUS_OK) {
                                                vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rc_eff);
                                            }
                                }
                                else {
                                    // Wet soil case
                                    const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[0];

                                    depac_wrapper(
                                            compnam,
                                            day_of_year,
                                            lat,
                                            temperature,
                                            ustar[ij],
                                            glrad,
                                            sinphi,
                                            rh,
                                            lai[ij],
                                            sai,
                                            1,  // nwet = 1 for wet conditions
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
                                            catm,
                                            &rc_eff,
                                            &gsoil_eff_out,
                                            &rsoil_eff_out,
                                            pressure
                                                );

                                            if (status == STATUS_OK) {
                                                vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rsoil_eff_out);
                                            }
                                }
                            }
                    }
                }
}

template<typename TF>
Deposition<TF>::Deposition(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_deposition = inputin.get_item<bool>("deposition", "swdeposition", "", false);

    // Added: Initialize DEPAC parameters for NH3 deposition

    // Radiation parameters
    glrad = inputin.get_item<TF>("deposition", "glrad", "", (TF)600.0);                // Global radiation (W/m2)
    sinphi = inputin.get_item<TF>("deposition", "sinphi", "", (TF)0.75);               // Sine of solar elevation:  Solar elevation angle at noon ≈ 48.5° ==>  sinphi = sin(48.5°) ≈ 0.75

    // Meteorological parameters
    temperature = inputin.get_item<TF>("deposition", "temperature", "", (TF)293.15);  // Air temperature (K)
    rh = inputin.get_item<TF>("deposition", "rh", "", (TF)50.0);                      // Relative humidity (%)

    // Surface parameters
    sai = inputin.get_item<TF>("deposition", "sai", "", (TF)6.0);                     // Stem area index (m2/m2)
    lat = inputin.get_item<TF>("deposition", "lat", "", (TF)51.0);                    // Latitude (degrees)

    // Time and surface condition parameters
    day_of_year = inputin.get_item<int>("deposition", "day_of_year", "", 217);        // Day of year: 18 August
    nwet = inputin.get_item<int>("deposition", "nwet", "", 0);                        // Surface wetness indicator
    lu = inputin.get_item<int>("deposition", "lu", "", 4);                            // Land use type

    // NH3-specific parameters
    iratns = inputin.get_item<int>("deposition", "iratns", "", 2);                    // NH3 compensation point option
    hlaw = inputin.get_item<TF>("deposition", "hlaw", "", (TF)6.1e4);                 //rmes = 1/(henry/3000.+100.*react)  ! Wesely '89, eq. 6
    react = inputin.get_item<TF>("deposition", "react", "", (TF)0.0);                 // Reactivity factor
    c_ave_prev_nh3 = inputin.get_item<TF>("deposition", "c_ave_prev_nh3", "", (TF)5.0); // Previous NH3 concentration (μg/m3)
    catm = inputin.get_item<TF>("deposition", "catm", "", (TF)0.76);                  // Atmospheric NH3 concentration (μg/m3) (0.76 μg/m3 ~ 1 ppb)
    pressure = inputin.get_item<TF>("deposition", "pressure", "", (TF)1.013e5);  // Default sea level pressure

}


template <typename TF>
Deposition<TF>::~Deposition()
{
}


template <typename TF>
void Deposition<TF>::init(Input& inputin)
{
    // Always read the default deposition velocities. They are needed by 
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

    auto& gd = grid.get_grid_data();

    // Create surface tiles for deposition:
    for (auto& name : deposition_tile_names)
        deposition_tiles.emplace(name, Deposition_tile<TF>{});

    for (auto& tile : deposition_tiles)
    {
        // tile.second.vdo3.resize(gd.ijcells);
        // tile.second.vdno.resize(gd.ijcells);
        // tile.second.vdno2.resize(gd.ijcells);
        // tile.second.vdhno3.resize(gd.ijcells);
        // tile.second.vdh2o2.resize(gd.ijcells);
        // tile.second.vdrooh.resize(gd.ijcells);
        // tile.second.vdhcho.resize(gd.ijcells);
        tile.second.vdnh3.resize(gd.ijcells);  // Added NH3
    }

    deposition_tiles.at("veg" ).long_name = "vegetation";
    deposition_tiles.at("soil").long_name = "bare soil";
    deposition_tiles.at("wet" ).long_name = "wet skin";
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

    // Modified arrays for NH3 only
    rmes     = {(TF)0.0};  // NH3
    rsoil    = {(TF)400.0};  // NH3
    rcut     = {(TF)1e5};  // NH3
    rws      = {(TF)2000.0};  // NH3
    rwat     = {(TF)2000.0};  // NH3
    diff     = {(TF)0.13};  // NH3
    diff_scl = {(TF)1.6};  // NH3
    henry    = {(TF)2e4};  // NH3
    f0       = {(TF)1.0};  // NH3

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
            "vdnh3_soil", "vdnh3_wet", "vdnh3_veg"};  // Modified for NH3 only
        cross_list = cross.get_enabled_variables(allowed_crossvars);
    }
}

template <typename TF>
void Deposition<TF>::update_time_dependent(
        Timeloop<TF>& timeloop,
        Boundary<TF>& boundary,
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

    // get information from lsm:
    auto& tiles = boundary.get_tiles();
    auto& lai = boundary.get_lai();
    auto& water_mask = boundary.get_water_mask();
    auto& c_veg = boundary.get_c_veg();

    // calculate deposition per tile:
    for (auto& tile : tiles)
    {
        calc_deposition_per_tile(
                tile.first,
                // deposition_tiles.at(tile.first).vdo3.data(),
                // deposition_tiles.at(tile.first).vdno.data(),
                // deposition_tiles.at(tile.first).vdno2.data(),
                // deposition_tiles.at(tile.first).vdhno3.data(),
                // deposition_tiles.at(tile.first).vdh2o2.data(),
                // deposition_tiles.at(tile.first).vdrooh.data(),
                // deposition_tiles.at(tile.first).vdhcho.data(),
                deposition_tiles.at(tile.first).vdnh3.data(),  // Added NH3
                lai.data(),
                c_veg.data(),
                tile.second.rs.data(),
                tiles.at("veg").rs.data(),
                tile.second.ra.data(),
                tile.second.ustar.data(),
                tile.second.fraction.data(),
                //rmes.data(), rsoil.data(), rcut.data(),
                //rws.data(), rwat.data(),
                diff_scl.data(),   
                // Added: DEPAC parameters
                glrad,          // Global radiation
                sinphi,         // Sine of solar elevation
                temperature,    // Air temperature
                rh,            // Relative humidity
                sai,           // Stem area index
                lat,           // Latitude
                day_of_year,   // Day of year
                nwet,          // Surface wetness
                lu,            // Land use type
                iratns,        // NH3 compensation point option
                hlaw,          // Henry's law constant
                react,         // Reactivity factor
                c_ave_prev_nh3, // Previous NH3 concentration
                catm,          // Atmospheric NH3 concentration
                pressure,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
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

