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
            float p,  // Added pressure parameter
            float *gw_out,           // Add these parameters
            float *gstom_out,
            float *cw_out,
            float *cstom_out,
            float *csoil_out,
            bool use_input_ccomp  // Added flag to use input ccomp_tot

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
                const bool sw_override_ccomp,        // NEW parameter
                const TF ccomp_override_value,       // NEW parameter
                Deposition_tile_map<TF>& deposition_tiles, // NEW parameter
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
                    const TF c_ug = TF(1.0e9) * rho[kstart] * xmnh3 * xmair_i;   // mol/mol to ug/m3


                    // Define component name outside the loops (doesn't change)
                    char compnam[4] = "NH3";

                    if (lu_type == "veg") {
                        // Vegetation tile handling
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
                                // NEW: Automatic determination of land use type and SAI based on LAI
                                // This allows DEPAC to use different parameters for grass vs forest
                                int local_lu;
                                TF local_sai;

                                if (lai[ij] <= 3.5) {
                                    local_lu = 1;  // grass
                                    local_sai = lai[ij];  // For grass, SAI = LAI
                                } else {
                                    local_lu = 5;  // deciduous forest
                                    local_sai = lai[ij] + 1.0;  // For forest, add stem area
                                }

                                //// Calculate SAI based on the land use type
                                //if (local_lu == 4 || local_lu == 5 || local_lu == 17 || local_lu == 18) {  // Forest types
                                //    local_sai = lai[ij] + 1.0;
                                //} else if (local_lu == 3) {  // Permanent crops
                                //    local_sai = lai[ij] + 0.5;
                                //} else {  // Default case (includes grass)
                                //    local_sai = lai[ij];
                                //}

                                // Keep IFS Ra and use vegetation Rb scaling
                                const TF rb = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[0];

                                //const TF nh3_ugm3 = nh3_concentration[ijk] * xmnh3 / 22.414 * 1.0e9; //mol/mol to ug/m3 conversion(STP)
                                const TF nh3_ugm3 = nh3_concentration[ijk] * c_ug; // mol/mol to ug/m3 conversion


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

                                // Call DEPAC wrapper for dry vegetation
                                //char compnam[4] = "NH3";
                                float rc_tot, ccomp_tot=0.0, rc_eff;
                                float gsoil_eff_out, rsoil_eff_out;
                                float gw_out, gstom_out;            // Added: conductance variables
                                float cw_out, cstom_out, csoil_out; // Added: compensation point variables
                                int status;

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
                                        temperature -273.15,
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
                                        c_ave_prev_nh3 * c_ug,
                                        ra[ij],
                                        rb,
                                        nh3_ugm3,
                                        //catm,
                                        &rc_eff,
                                        &gsoil_eff_out,
                                        &rsoil_eff_out,
                                        pressure,
                                        &gw_out,            // Added output variable
                                        &gstom_out,         // Added output variable
                                        &cw_out,            // Added output variable
                                        &cstom_out,         // Added output variable
                                        &csoil_out,
                                        use_input_ccomp  // Pass the flag
                                            );

                                        // Calculate deposition velocity using resistance analogy
                                        if (status == STATUS_OK) {

                                            // Store ccomp_tot value
                                            deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot;

                                            // Store resistances (inverting conductances)
                                            deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? 
                                                (TF)1.0 / gw_out : (TF)9999.0;

                                            deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? 
                                                (TF)1.0 / gstom_out : (TF)9999.0;

                                            deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? 
                                                (TF)1.0 / gsoil_eff_out : (TF)9999.0;

                                            // Store compensation points directly
                                            deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out;
                                            deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out;
                                            deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out;
                                            deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot;
                                            deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff;


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

                                //const TF nh3_ugm3 = nh3_concentration[ijk] * xmnh3 / 22.414 * 1.0e9; // mol/mol to ug/m3 conversion(STP)
                                const TF nh3_ugm3 = nh3_concentration[ijk] * c_ug; // mol/mol to ug/m3 conversion

                                // Call DEPAC wrapper for dry soil
                                //char compnam[4] = "NH3";
                                float rc_tot, ccomp_tot=0.0, rc_eff;
                                float gsoil_eff_out, rsoil_eff_out;
                                float gw_out, gstom_out;            // Added: conductance variables
                                float cw_out, cstom_out, csoil_out; // Added: compensation point variables
                                int status;

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
                                        temperature -273.15,
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
                                        c_ave_prev_nh3 * c_ug,
                                        ra[ij],
                                        rb,
                                        nh3_ugm3,
                                        //catm,
                                        &rc_eff,
                                        &gsoil_eff_out,
                                        &rsoil_eff_out,
                                        pressure,
                                        &gw_out,            // Added output variable
                                        &gstom_out,         // Added output variable
                                        &cw_out,            // Added output variable
                                        &cstom_out,         // Added output variable
                                        &csoil_out,
                                        use_input_ccomp  // Pass the flag
                                            );

                                if (status == STATUS_OK) {

                                    // Store ccomp_tot value
                                    deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot;

                                    // Store resistances (inverting conductances)
                                    deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? 
                                        (TF)1.0 / gw_out : (TF)9999.0;

                                    deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? 
                                        (TF)1.0 / gstom_out : (TF)9999.0;

                                    deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? 
                                        (TF)1.0 / gsoil_eff_out : (TF)9999.0;

                                    // Store compensation points directly
                                    deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out;
                                    deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out;
                                    deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out;
                                    deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot;
                                    deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff;


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

                                const TF nh3_ugm3 = nh3_concentration[ijk] * c_ug; // mol/mol to ug/m3 conversion

                                // std::cout << "VEG tile: i=" << i << ", j=" << j << ", ijk=" << ijk << std::endl;
                                // std::cout << "  NH3 conc = " << nh3_concentration[ijk]
                                //     << ", glrad = " << glrad
                                //     << ", rh = " << rh
                                //     << ", sinphi = " << sinphi << std::endl;

                                if (fraction[ij] < (TF)1e-12)
                                    continue;

                                //char compnam[4] = "NH3";
                                float rc_tot, ccomp_tot=0.0, rc_eff;
                                float gsoil_eff_out, rsoil_eff_out;
                                float gw_out, gstom_out;            // Added: conductance variables
                                float cw_out, cstom_out, csoil_out; // Added: compensation point variables
                                int status;

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

                                    //const TF nh3_ugm3 = nh3_concentration[ijk] * xmnh3 / 22.414 * 1.0e9; // mol/mol to ug/m3 conversion(STP)
                                    //const TF nh3_ugm3 = nh3_concentration[ijk] * c_ug; // mol/mol to ug/m3 conversion


                                    depac_wrapper(
                                            compnam,
                                            day_of_year,
                                            lat,
                                            temperature -273.15,
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
                                            c_ave_prev_nh3 * c_ug,
                                            ra[ij],
                                            rb,
                                            nh3_ugm3,
                                            //catm,
                                            &rc_eff,
                                            &gsoil_eff_out,
                                            &rsoil_eff_out,
                                            pressure,
                                            &gw_out,            // Added output variable
                                            &gstom_out,         // Added output variable
                                            &cw_out,            // Added output variable
                                            &cstom_out,         // Added output variable
                                            &csoil_out,
                                            use_input_ccomp  // Pass the flag
                                                );

                                            if (status == STATUS_OK) {

                                                // Store ccomp_tot value
                                                deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot;

                                                // Store resistances (inverting conductances)
                                                deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? 
                                                    (TF)1.0 / gw_out : (TF)9999.0;

                                                deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? 
                                                    (TF)1.0 / gstom_out : (TF)9999.0;

                                                deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? 
                                                    (TF)1.0 / gsoil_eff_out : (TF)9999.0;

                                                // Store compensation points directly
                                                deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out;
                                                deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out;
                                                deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out;
                                                deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot;
                                                deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff;


                                                vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rc_eff);
                                            }
                                }
                                else {
                                    // Wet soil case
                                    const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[0];

                                    //const TF nh3_ugm3 = nh3_concentration[ijk] * xmnh3 / 22.414 * 1.0e9; //mol/mol to ug/m3 conversion(STP)
                                    //const TF nh3_ugm3 = nh3_concentration[ijk] * c_ug; // mol/mol to ug/m3 conversion

                                    depac_wrapper(
                                            compnam,
                                            day_of_year,
                                            lat,
                                            temperature -273.15,
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
                                            c_ave_prev_nh3 * c_ug,
                                            ra[ij],
                                            rb,
                                            nh3_ugm3,
                                            //catm,
                                            &rc_eff,
                                            &gsoil_eff_out,
                                            &rsoil_eff_out,
                                            pressure,
                                            &gw_out,            // Added output variable
                                            &gstom_out,         // Added output variable
                                            &cw_out,            // Added output variable
                                            &cstom_out,         // Added output variable
                                            &csoil_out,
                                            use_input_ccomp  // Pass the flag
                                                );

                                    if (status == STATUS_OK) {
                                        // Store ccomp_tot value
                                        deposition_tiles.at(lu_type).ccomp_tot.data()[ij] = ccomp_tot;

                                        // Store resistances (inverting conductances)
                                        deposition_tiles.at(lu_type).cw.data()[ij] = (gw_out > 0.0) ? 
                                            (TF)1.0 / gw_out : (TF)9999.0;

                                        deposition_tiles.at(lu_type).cstom.data()[ij] = (gstom_out > 0.0) ? 
                                            (TF)1.0 / gstom_out : (TF)9999.0;

                                        deposition_tiles.at(lu_type).csoil_eff.data()[ij] = (gsoil_eff_out > 0.0) ? 
                                            (TF)1.0 / gsoil_eff_out : (TF)9999.0;

                                        // Store compensation points directly
                                        deposition_tiles.at(lu_type).cw_out.data()[ij] = cw_out;
                                        deposition_tiles.at(lu_type).cstom_out.data()[ij] = cstom_out;
                                        deposition_tiles.at(lu_type).csoil_out.data()[ij] = csoil_out;
                                        deposition_tiles.at(lu_type).rc_tot.data()[ij] = rc_tot;
                                        deposition_tiles.at(lu_type).rc_eff.data()[ij] = rc_eff;


                                        vdnh3[ij] = (TF)1.0 / (ra[ij] + rb + rsoil_eff_out);
                                    }
                                }
                            }
                    }
                }
}


template<typename TF>
Deposition<TF>::Deposition(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, 
        Radiation<TF>& radiationin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), radiation(radiationin)
{
    sw_deposition = inputin.get_item<bool>("deposition", "swdeposition", "", false);

    // Added: Initialize DEPAC parameters for NH3 deposition

    // Get start hour from input

    // Radiation parameters
    //glrad = inputin.get_item<TF>("deposition", "glrad", "", (TF)400.0);                // Global radiation (W/m2)
    //sinphi = inputin.get_item<TF>("deposition", "sinphi", "", (TF)0.75);               // Sine of solar elevation:  Solar elevation angle at noon ≈ 48.5° ==>  sinphi = sin(48.5°) ≈ 0.75

    // Meteorological parameters
    //temperature = inputin.get_item<TF>("deposition", "temperature", "", (TF)293.15);  // Air temperature (K)
    //rh = inputin.get_item<TF>("deposition", "rh", "", (TF)50.0);                      // Relative humidity (%)

    // Surface parameters
    //sai = inputin.get_item<TF>("deposition", "sai", "", (TF)6.0);                     // Stem area index (m2/m2)
    //lat = inputin.get_item<TF>("deposition", "lat", "", (TF)52.3);                    // Latitude (degrees)

    // Time and surface condition parameters
    //day_of_year = inputin.get_item<int>("deposition", "day_of_year", "", 264);        // Day of year: 20 September
    nwet = inputin.get_item<int>("deposition", "nwet", "", 0);                        // Surface wetness indicator
                                                                                      //lu = inputin.get_item<int>("deposition", "lu", "", 5);                            // Land use type

                                                                                      // NH3-specific parameters
    iratns = inputin.get_item<int>("deposition", "iratns", "", 2);                    // NH3 compensation point option
    hlaw = inputin.get_item<TF>("deposition", "hlaw", "", (TF)6.1e4);                 //rmes = 1/(henry/3000.+100.*react)  ! Wesely '89, eq. 6
    react = inputin.get_item<TF>("deposition", "react", "", (TF)0.0);                 // Reactivity factor
    c_ave_prev_nh3 = inputin.get_item<TF>("deposition", "c_ave_prev_nh3", "", (TF)6.735e-9); // Previous NH3 concentration (mol/mol, then it converts to ug/m3)
                                                                                             //catm = inputin.get_item<TF>("deposition", "catm", "", (TF)0.76);                  // Atmospheric NH3 concentration (μg/m3) (0.76 μg/m3 ~ 1 ppb)
    pressure = inputin.get_item<TF>("deposition", "pressure", "", (TF)1.013e5);  // Default sea level pressure

    sw_override_ccomp = inputin.get_item<bool>("deposition", "sw_override_ccomp", "", false);
    ccomp_override_value = inputin.get_item<TF>("deposition", "ccomp_override_value", "", TF(0.0));

    //// Debug print
    //if (sw_override_ccomp) {
    //    master.print_message("DEPAC: Compensation point override ENABLED. Using value: %f\n", ccomp_override_value);
    //} else {
    //    master.print_message("DEPAC: Compensation point override DISABLED. Using calculated values.\n");
    //}

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
        tile.second.vdnh3.resize(gd.ijcells);  // allocate memory for arrays
        tile.second.ra.resize(gd.ijcells);
        tile.second.obuk.resize(gd.ijcells);
        tile.second.ustar.resize(gd.ijcells);
        tile.second.ccomp_tot.resize(gd.ijcells);
        tile.second.cw.resize(gd.ijcells);
        tile.second.cstom.resize(gd.ijcells);
        tile.second.csoil_eff.resize(gd.ijcells);

        // Resize new compensation point arrays
        tile.second.cw_out.resize(gd.ijcells);
        tile.second.cstom_out.resize(gd.ijcells);
        tile.second.csoil_out.resize(gd.ijcells);
        tile.second.rc_tot.resize(gd.ijcells);
        tile.second.rc_eff.resize(gd.ijcells);
    }
    // Initialize grid-mean arrays
    ra_mean.resize(gd.ijcells);
    ccomp_mean.resize(gd.ijcells);
    cw_mean.resize(gd.ijcells);
    cstom_mean.resize(gd.ijcells);
    csoil_eff_mean.resize(gd.ijcells);
    cw_out_mean.resize(gd.ijcells);
    cstom_out_mean.resize(gd.ijcells);
    csoil_out_mean.resize(gd.ijcells);
    rc_tot_mean.resize(gd.ijcells);
    rc_eff_mean.resize(gd.ijcells);

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
            "vdnh3_soil", "vdnh3_wet", "vdnh3_veg",
            "ra_soil", "ra_wet", "ra_veg",
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
            "rc_eff", "rc_eff_veg", "rc_eff_soil", "rc_eff_wet"
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

    // Get current model time
    const TF model_time = timeloop.get_time();

    // Calculate actual time of day 
    const TF actual_time = t0 + model_time;

    const std::vector<TF>& rho = thermo.get_basestate_vector("rho");

    // Get day of year from Timeloop
    day_of_year = int(timeloop.calc_day_of_year());

    // Get latitude from Grid
    lat = gd.lat;

    const int year = timeloop.get_year();
    const TF seconds_after_midnight = TF(timeloop.calc_hour_of_day() * 3600);
    TF azimuth;

    // Calculate sinphi using the radiation function
    std::tie(sinphi, azimuth) = Radiation_rrtmgp_functions::calc_cos_zenith_angle(
            lat, gd.lon, day_of_year, seconds_after_midnight, year);

    //master.print_message("DEBUG: Time step %f, Hour of day %f, sinphi (from radiation) = %f\n", 
    //        timeloop.get_time(), 
    //        timeloop.calc_hour_of_day(),
    //        sinphi);

    //master.print_message("DEBUG: About to access radiation, hour = %f\n", timeloop.calc_hour_of_day());


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


    auto tmp2 = fields.get_tmp();
    if (tmp2 && tmp2->fld.data()) {
        thermo.get_thermo_field(*tmp2, "T", true, false);
        // Get temperature at proper surface level
        temperature = tmp2->fld.data()[gd.kstart*gd.ijcells];
        fields.release_tmp(tmp2);
    }

    auto tmp1 = fields.get_tmp();
    if (tmp1 && tmp1->fld.data()) {
        thermo.get_thermo_field(*tmp1, "rh", true, false);
        // Get RH at proper surface level
        rh = tmp1->fld.data()[gd.kstart*gd.ijcells] * 100.0;
        fields.release_tmp(tmp1);
    }

    // get information from lsm:
    auto& tiles = boundary.get_tiles();
    auto& lai = boundary.get_lai();
    auto& water_mask = boundary.get_water_mask();
    auto& c_veg = boundary.get_c_veg();
    //auto& sw_flux_dn = radiation.get_c_veg(); //how to call a variable from a remote scheme (maarten)

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
        calc_deposition_per_tile(
                master,
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
                fields.sp.at("nh3")->fld.data(),  // Pass NH3 concentration directly from Fields
                                                  //rmes.data(), rsoil.data(), rcut.data(),
                                                  //rws.data(), rwat.data(),
                diff_scl.data(),   
                rho.data(),
                // Added: DEPAC parameters
                glrad,          // Now using calculated time-dependent radiation
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
                sw_override_ccomp,              // NEW argument
                ccomp_override_value,           // NEW argument
                deposition_tiles,               // NEW argument
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

    get_tiled_mean(ra_mean.data(), "ra", (TF)1.0, tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
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
        else if (name == "ccomp_tot")
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
    else if (name == "ra") {
        fld_veg  = deposition_tiles.at("veg").ra.data();
        fld_soil = deposition_tiles.at("soil").ra.data();
        fld_wet  = deposition_tiles.at("wet").ra.data();
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

