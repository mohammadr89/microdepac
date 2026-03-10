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
 * === INPUTS (*.ini file) ===
 * [chemistry]
 * swchemistry  = boolean : Enable/disable chemistry module (default: false)
 * lifetime     = float   : Tracer decay timescale [s] (default: 1e30)
<<<<<<< Updated upstream
 * z_target     = float   : Target height for concentration extrapolation [m] (required)
 * rsl_ratio    = float   : Roughness sublayer ratio for MO validity (default: 20.0)
=======
 * c_extrap_diff : Concentration difference between z_target and first grid level [ppb]
 * rsl_ratio    = float   : Roughness sublayer ratio (default: 20.0, only used if sw_const_ref_height=false)
 * sw_const_ref_height = boolean : Use constant reference height (default: true)
 * z_fixed            = float   : Fixed reference height [m] (default: 20, used if sw_const_ref_height=true)
 * sw_adapt_ref_height = boolean : Use adaptive extrapolation method (default: true)
 * sw_use_lowest_levels = boolean : Use kstart/kstart+1 for cstar (true) or bracketing levels (false) (default: true)
>>>>>>> Stashed changes
 * 
 * === INPUTS (NetCDF file: timedep_chem group) ===
 * Photolysis rates: jo31d, jh2o2, jno2, jno3, jn2o5, jch2or, jch2om, jch3o2h
 * Emissions: emi_isop, emi_no
 * All with time dimension: time_chem*
 * 
 * === OUTPUTS (Statistics) ===
 * vdnh3              : NH3 deposition velocity [m s-1]
 * flux_inst          : Instantaneous NH3 flux [kg m-2 s-1]
 * total_flux_mol_ha  : NH3 cumulative flux (total) [mol ha-1]
 * cstar1             : Concentration scaling (with stability) [mol mol-1]
 * cstar2             : Concentration scaling (neutral) [mol mol-1]
 * c_diff             : Stability effect on extrapolated concentration [ppb]
 * c_target           : NH3 at target height (optimal) [mol mol-1]
<<<<<<< Updated upstream
 * c_diff_flux        : Concentration difference flux [kg m-2 s-1]
=======
>>>>>>> Stashed changes
 * chem_budget        : Chemistry budget per layer [molecules cm-3 s-1]
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
#include "chemistry.h"
#include "constants.h"
#include "timeloop.h"
#include "deposition.h"
#include "boundary.h"
#include "cross.h"
#include <vector> 
#include "monin_obukhov.h"


/**
 * Calculate general factor for any two heights
 */
template<typename TF>
inline TF calc_factor(const TF z1, const TF z2, const TF L)
{
    // Calculate stability corrections
    TF psi1 = TF(0), psi2 = TF(0);
    
    if (std::abs(L) > TF(1e-15))
    {
        namespace most = Monin_obukhov;
        
        // Calculate z/L ratios
        const TF z1_over_L_raw = z1 / L;
        const TF z2_over_L_raw = z2 / L;
        
        // Apply MicroHH's capping: z/L between -10,000 and +10
        const TF z1_over_L = std::min(std::max(z1_over_L_raw, Constants::zL_min<TF>), 
                                       Constants::zL_max<TF>);
        const TF z2_over_L = std::min(std::max(z2_over_L_raw, Constants::zL_min<TF>), 
                                       Constants::zL_max<TF>);
        
        psi1 = (z1_over_L <= TF(0)) ? 
            most::psih_unstable(z1_over_L) : most::psih_stable(z1_over_L);
        psi2 = (z2_over_L <= TF(0)) ? 
            most::psih_unstable(z2_over_L) : most::psih_stable(z2_over_L);
    }
    
    return (std::log(z2/z1) - psi2 + psi1) / Constants::kappa<TF>;
}

namespace
{
    double CFACTOR;                          /* Conversion factor for concentration units */
    std::pair<std::string, int> check_for_unique_time_dim(const std::map<std::string, int>& dims)
    {
        // Check for the existence of a unique time dimension.
        bool only_one_time_dim = false;
        std::string time_dim;
        int time_dim_length = 0;

        for (auto i : dims)
        {
            if (i.first.substr(0, 9) == "time_chem")
            {
                if (only_one_time_dim)
                    throw std::runtime_error("More than one time dimensions in input");
                else
                {
                    only_one_time_dim = true;
                    time_dim = i.first;
                    time_dim_length = i.second;
                }
            }
        }

        return std::make_pair(time_dim, time_dim_length);
    }

    // Updated PSS function with concentration scaling calculations
    /**
     * Find closest grid point below target height
     */
    template<typename TF>
    inline int find_ref_below(const TF* const z, const TF z_target, 
                             const int kstart, const int kend)
    {
        int k_ref = kstart;
        for (int k = kstart; k < kend - 1; ++k)
        {
            if (z[k] <= z_target) k_ref = k;
            else break;
        }
        return k_ref;
    }
    
    template<typename TF>
    void pss(
        TF* restrict tnh3,
        const TF* const restrict nh3,
        const TF* const restrict jval,
        const TF* const restrict emval,
        const TF* const restrict vdnh3,
        const TF* const restrict tprof,
        const TF* const restrict qprof,
        const TF* const restrict dzi,
        const TF* const restrict rhoref,
        const TF* const restrict z,       // Add grid height levels
        const TF* const restrict obuk,    // Add Obukhov length (from surface)
        const TF* const restrict z0m,     // Add surface roughness length
        TF* restrict rfa,
        TF* restrict flux_nh3, 
        TF* restrict flux_inst, 
        TF* restrict total_flux_nh3,
        TF* restrict cstar1,  
        TF* restrict cstar2,  
<<<<<<< Updated upstream
        TF* restrict c_target, 
=======
        TF* restrict c_target,
        TF* restrict c_extrap_diff,
        TF* restrict c_diff,     
>>>>>>> Stashed changes
        TF& trfa,
        const TF dt,
        const TF sdt,
        const TF lifetime,
        // const TF z_target,
<<<<<<< Updated upstream
        const TF rsl_ratio,
=======
        const bool sw_const_ref_height,
        const TF z_fixed,              
        const bool sw_adapt_ref_height,
        const bool sw_use_lowest_levels,
>>>>>>> Stashed changes
        const int istart, const int iend,
        const int jstart, const int jend,
        const int kstart, const int kend,
        const int jstride, const int kstride,
        const TF dx,
        const TF dy)
    {
        const TF xmh2o = 18.015265;
        const TF xmnh3 = 17.031;
        const TF xmh2o_i = TF(1) / xmh2o;
        const TF xmair = 28.9647;       // Molar mass of dry air  [kg kmol-1]
        const TF xmair_i = TF(1) / xmair;
        const TF Na = 6.02214086e23; // Avogadros number [molecules mol-1]
    
        // Update the time integration of the reaction fluxes with the full timestep on first RK3 step
        if (abs(sdt/dt - 1./3.) < 1e-5) trfa += dt;
    
        // Import Monin_obukhov namespace for stability functions
        namespace most = Monin_obukhov;
    
        for (int k=kstart; k<kend; ++k)
        {
            const TF C_M = TF(1e-3) * rhoref[k] * Na * xmair_i;   // molecules/cm3 for chmistry!
    
            // From ppb (units mixing ratio) to molecules/cm3 --> changed: now mol/mol unit for transported tracers:
            const TF CFACTOR = C_M;
            const TF sdt_cfac_i = TF(1) / (sdt * CFACTOR);
            const TF lti = TF(1)/lifetime;  // 1/s
            TF decay;
            for (int j=jstart; j<jend; ++j)
            #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    const int ij = i + j*jstride;
    
                    // kg/kg --> molH2O/molAir --*C_M--> molecules/cm3 limit to 1 molecule/cm3 to avoid error usr_HO2_HO2
                    // const TF C_H2O = std::max(qt[ijk] * xmair * C_M * xmh2o_i, TF(1));
                    // const TF TEMP = temp[ijk];
    
                    if (k==kstart)
                    {
                        // Add new concentration scaling calculations
                        // Get concentrations at two vertical levels (kstart and kstart+1)
                        const int ijk1 = i + j*jstride + kstart*kstride;
                        const int ijk2 = i + j*jstride + (kstart+1)*kstride;
                        
                        const TF c_1 = nh3[ijk1];
                        const TF c_2 = nh3[ijk2];
                        
                        // Heights at the two levels
                        const TF z_1 = z[kstart];
                        const TF z_2 = z[kstart+1];
                        
                        // Obukhov length from surface
                        const TF L = obuk[ij];
                        
                        const TF gradient_factor = calc_factor(z_1, z_2, L);       // For calculating c* from gradient
                        const TF neutral_factor = calc_factor(z_1, z_2, TF(1e30)); // For neutral atmosphere (no stability)
                        
                        // Calculate cstar2 using simple logarithmic formula (no stability correction)
                        if (std::abs(neutral_factor) > TF(1e-15))
                        {
                            cstar2[ij] = +1.0 * (c_2 - c_1) / neutral_factor;
                        }
                        else
                        {
                            cstar2[ij] = TF(0.0);  // Avoid division by zero
                        }
                        
<<<<<<< Updated upstream
                        // Calculate cstar1 with stability correction
                        if (std::abs(gradient_factor) > TF(1e-15))
                        {
                            cstar1[ij] = +1.0 * (c_2 - c_1) / gradient_factor;
                        }
                        else
                        {
                            cstar1[ij] = TF(0.0);  // Add this else clause for safety
=======
                        // ========================================
                        // STEP 2: Check conditions for SIMPLIFIED method
                        // ========================================
                        // if (z_1 >= z_target || !sw_adapt_ref_height)
                        // {
                            // SIMPLIFIED METHOD: Use first grid level directly
                            // Condition 1: First grid is already above target (regardless of sw_adapt_ref_height)
                            // OR
                            // Condition 2: Adaptive method is disabled

                        bool use_simplified;
                        if (!sw_adapt_ref_height)
                            use_simplified = true;
                        else if (sw_const_ref_height)
                            use_simplified = false;  // always extrapolate to z_fixed for fair comparison
                        else
                            use_simplified = (z_1 >= z_target);
                        
                        if (use_simplified)
                        {
                            c_target[ij]      = c_1;
                            cstar1[ij]        = TF(0);
                            cstar2[ij]        = TF(0);
                            c_diff[ij]        = TF(0);
                            c_extrap_diff[ij] = TF(0);

                            // Calculate instantaneous flux using first grid level concentration (SIMPLIFIED)
                            flux_inst[ij] = (-1.0) * vdnh3[ij] * nh3[ijk] * rhoref[k] * xmair_i * xmnh3; // [kg(NH3) m-2 s-1]
                            
                            decay = vdnh3[ij]*dzi[k] + lti;   // 1/s (SIMPLIFIED - no concentration scaling)
                        }
                        else
                        {
                            // ========================================
                            // STEP 3: ADAPTIVE METHOD with extrapolation
                            // Only reaches here if: z_1 < z_target AND sw_adapt_ref_height == true
                            // ========================================
                            
                            // Get second level concentration
                            const int ijk2 = i + j*jstride + (kstart+1)*kstride;
                            const TF c_2 = nh3[ijk2];
                            
                            // Heights at the two levels
                            const TF z_2 = z[kstart+1];
                            
                            // Obukhov length from surface
                            const TF L = obuk[ij];
                            
                            // ========================================
                            // Determine which levels to use for cstar
                            // ========================================
                            int k_below, k_above;
                            if (sw_use_lowest_levels)
                            {
                                // Always use the two lowest levels
                                k_below = kstart;
                                k_above = kstart + 1;
                            }
                            else
                            {
                                // Use levels bracketing z_target
                                k_below = kstart;
                                k_above = kstart + 1;
                                for (int k_grid = kstart; k_grid < kend - 1; ++k_grid)
                                {
                                    if (z[k_grid] < z_target)
                                    {
                                        k_below = k_grid;
                                        k_above = k_grid + 1;
                                    }
                                    else break;
                                }
                            }
                            
                            
                            const int ijk_below = i + j*jstride + k_below*kstride;
                            const int ijk_above = i + j*jstride + k_above*kstride;
                            const TF c_below = nh3[ijk_below];
                            const TF c_above = nh3[ijk_above];

                            // ========================================
                            // Check if a grid level sits exactly at z_target
                            // ========================================
                            if (std::abs(z[k_below] - z_target) < TF(0.5) || std::abs(z[k_above] - z_target) < TF(0.5))
                            {
                                int k_exact = (std::abs(z[k_above] - z_target) < std::abs(z[k_below] - z_target))
                                              ? k_above : k_below;
                                const int ijk_exact = i + j*jstride + k_exact*kstride;
                                c_target[ij]      = nh3[ijk_exact];
                                cstar1[ij]        = TF(0);
                                cstar2[ij]        = TF(0);
                                c_diff[ij]        = TF(0);
                                c_extrap_diff[ij] = (c_target[ij] - c_1) * TF(1e9); 
                                flux_inst[ij] = (-1.0) * vdnh3[ij] * c_target[ij] * rhoref[k] * xmair_i * xmnh3;
                                if (std::abs(nh3[ijk]) > TF(1e-15))
                                    decay = (vdnh3[ij] * dzi[k] * c_target[ij] / nh3[ijk]) + lti;
                                else
                                    decay = lti;
                            }
                            else
                            {
                                // Calculate cstar1 (with stability) and cstar2 (neutral)
                                const TF gradient_factor = calc_factor(z[k_below], z[k_above], L);
                                const TF neutral_factor  = calc_factor(z[k_below], z[k_above], TF(1e30));
                                
                                cstar1[ij] = (std::abs(gradient_factor) > TF(1e-15)) ?
                                    (c_above - c_below) / gradient_factor : TF(0);
                                
                                cstar2[ij] = (std::abs(neutral_factor) > TF(1e-15)) ?
                                    (c_above - c_below) / neutral_factor  : TF(0);
                                
                                // Extrapolate c_target to z_target using cstar1
                                const TF scaling_factor = calc_factor(z[k_below], z_target, L);
                                c_target[ij] = c_below + cstar1[ij] * scaling_factor;

                                c_extrap_diff[ij] = (c_target[ij] - c_1) * TF(1e9);

                                
                                // c_diff: effect of stability correction [ppb]
                                const TF scaling_factor_neutral = calc_factor(z[k_below], z_target, TF(1e30));
                                c_diff[ij] = (cstar1[ij] - cstar2[ij]) * scaling_factor_neutral * TF(1e9);
                                
                                // Calculate flux using c_target
                                flux_inst[ij] = (-1.0) * vdnh3[ij] * c_target[ij] * rhoref[k] * xmair_i * xmnh3;
                                
                                // Calculate decay with concentration scaling
                                if (std::abs(nh3[ijk]) > TF(1e-15))
                                    decay = (vdnh3[ij] * dzi[k] * c_target[ij] / nh3[ijk]) + lti;
                                else
                                    decay = lti;
                            }
>>>>>>> Stashed changes
                        }

                        const TF z_target = rsl_ratio * z0m[ij];

                        const TF scaling_factor = calc_factor(z_1, z_target, L);
                        c_target[ij] = c_1 + cstar1[ij] * scaling_factor;
                        
                        // Check if first grid level is high enough above roughness for MO theory
                        TF concentration_for_flux;
                        if (z_1 < rsl_ratio * z0m[ij])
                        {
                            // First grid level is too close to surface (inside roughness sublayer)
                            // MO theory not valid at z_1, so extrapolate upward to z_target
                            concentration_for_flux = c_target[ij];
                            
                            // // Debug print for first grid point only
                            // if (i == istart && j == jstart)
                            // {
                            //     std::printf("Time step: z_1=%.3f < 20*z0m=%.3f -> Using c_target=%.6e\n", z_1, rsl_ratio*z0m[ij], c_target[ij]);
                            // }
                        }
                        else
                        {
                            // First grid level is high enough above surface (outside roughness sublayer)
                            // MO theory already valid at z_1, use concentration directly
                            c_target[ij] = c_1;  // Update c_target for consistency
                            concentration_for_flux = c_1;
                            
                            // // Debug print for first grid point only
                            // if (i == istart && j == jstart)
                            // {
                            //     std::printf("Time step: z_1=%.3f >= 20*z0m=%.3f -> Using c_1=%.6e\n", z_1, rsl_ratio*z0m[ij], c_1);
                            // }
                        }

                        // Calculate and accumulate flux for this RK3 step
                        // Note: flux is accumulated (+=) and scaled by sdt
                        
                        flux_inst[ij] = (-1.0) * vdnh3[ij] * concentration_for_flux * rhoref[k] * xmair_i * xmnh3; // [kg(NH3) m-2 s-1]
    
                        // accumulate over sub-timestep for statistics (gets reset periodically)
                        TF flux = flux_inst[ij] * sdt; // [kg m⁻² s⁻¹] × [s] = [kg m⁻²] Scale by timestep for accumulation     
                        flux_nh3[ij] += flux;        // For period statistics - Accumulate [kg m⁻²]
                        
                        // accumulate total flux (never gets reset)
                        total_flux_nh3[ij] += flux;  // [kg m⁻²]
                        
                        // decay = vdnh3[ij]*dzi[k] + lti;   // 1/s
                        if (std::abs(nh3[ijk]) > TF(1e-15)) //to prevent division by zero
                        {
                            decay = (vdnh3[ij] * dzi[k] * concentration_for_flux / nh3[ijk]) + lti;   // 1/s
                        }
                        else
                        {
                        decay = lti; // 1/s
                        }
                    }
                    else
                    { 
                        decay = lti; // 1/s
                    }

                    // update tendencies:
                    tnh3[ijk] -= decay*nh3[ijk];
    
                    // Get statistics for reaction fluxes:
                    //if (abs(sdt/dt - 1./3.) < 1e-5)
                    //{
                    //    for (int l=0; l<NREACT; ++l)
                    //        rfa[(k-kstart)*NREACT+l] +=  RF[l]*dt;    // take the first evaluation in the RK3 steps, but with full time step.
                    //}
    
                    //  Reculculate tendency and add to the tendency of the transported tracers:
                } // i
        } // k
    }
}

template<typename TF>
// Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Radiation<TF>& radiationin, Input& inputin):
//     master(masterin), grid(gridin), fields(fieldsin), radiation(radiationin), field3d_operators(master, grid, fields) // added Radiation<TF>& radiationin parameter and member variable
Chemistry<TF>::Chemistry(
    Master& masterin, 
    Grid<TF>& gridin, 
    Fields<TF>& fieldsin,      // ← Parameter comes in
    Radiation<TF>& radiationin, 
    Input& inputin)
    : master(masterin)          // ← Step 1: master member is initialized
    , grid(gridin)              // ← Step 2: grid member is initialized  
    , fields(fieldsin)          // ← Step 3: fields member is initialized FROM fieldsin
    // fields (is the member variable) - fieldsin (is the parameter)
    , radiation(radiationin)    // ← Step 4: radiation member is initialized
    , field3d_operators(master, grid, fields)  // ← Step 5: Uses members (already set!) // added Radiation<TF>& radiationin parameter and member variable
{
    // Rest of the constructor remains the same
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();
    sw_chemistry = inputin.get_item<bool>("chemistry", "swchemistry", "", false);
    //lifetime     = inputin.get_item<TF>("chemistry", "lifetime", "", (TF)72000);  // seconds (20 hour default)
    lifetime     = inputin.get_item<TF>("chemistry", "lifetime", "", (TF)1e30);  // seconds
    // z_target = inputin.get_item<TF>("chemistry", "z_target", "");

    // Roughness sublayer ratio (Basu & Lacser 2017, DOI: 10.1007/s10546-016-0225-y)
    // Recommended: 50 (z1 > 50*z0), EPA: 20, Wind engineering: 30
    rsl_ratio = inputin.get_item<TF>("chemistry", "rsl_ratio", "", (TF)20.0);

<<<<<<< Updated upstream
    master.print_message("Lifetime of the tracer:  = %13.5e s \n", lifetime);
    // master.print_message("Target height z_target = %13.5e m \n", z_target);
    master.print_message("Roughness sublayer ratio (rsl_ratio) = %13.5e \n", rsl_ratio);
=======
    sw_const_ref_height = inputin.get_item<bool>("chemistry", "sw_const_ref_height", "", true);
    z_fixed = inputin.get_item<TF>("chemistry", "z_fixed", "", (TF)20);
    sw_adapt_ref_height = inputin.get_item<bool>("chemistry", "sw_adapt_ref_height", "", true);
    sw_use_lowest_levels = inputin.get_item<bool>("chemistry", "sw_use_lowest_levels", "", true);
    master.print_message("Use lowest levels for cstar = %s \n", sw_use_lowest_levels ? "true" : "false");

    master.print_message("Lifetime of the tracer:  = %13.5e s \n", lifetime);
    // master.print_message("Target height z_target = %13.5e m \n", z_target);
    master.print_message("Roughness sublayer ratio (rsl_ratio) = %13.5e \n", rsl_ratio);
    master.print_message("Use constant reference height = %s \n", sw_const_ref_height ? "true" : "false");
    if (sw_const_ref_height)
        master.print_message("  Fixed reference height (z_fixed) = %13.5e m \n", z_fixed);
    master.print_message("Use adaptive reference height method = %s \n", sw_adapt_ref_height ? "true" : "false");
>>>>>>> Stashed changes

    if (!sw_chemistry)
        return;
    //deposition = std::make_shared<Deposition <TF>>(masterin, gridin, fieldsin, radiationin, inputin);

    deposition = std::make_shared<Deposition<TF>>(masterin, gridin, fieldsin, radiationin, *this, inputin);
    // *this passes the current Chemistry object as a reference to the Deposition constructor.
    // With *this: Deposition can access Chemistry's methods and data
    // Can call chemistry.get_vd("nh3")
    // Can read chemistry.lifetime
    // Can interact bidirectionally with Chemistry
    // Without *this: Deposition is independent
    // Cannot directly access Chemistry object
    // One-way relationship: Chemistry knows about Deposition, but not vice versa

    // deposition = std::make_shared<Deposition<TF>>(
    //     master,      // ← Member variable (already initialized in Step 1)
    //     grid,        // ← Member variable (already initialized in Step 2)
    //     fields,      // ← Member variable (already initialized in Step 3)
    //     radiation,   // ← Member variable (already initialized in Step 4)
    //     *this, 
    //     inputin);    // ← Still the parameter (not stored as member)
}

template <typename TF>
Chemistry<TF>::~Chemistry()
{
}

template<typename TF>
void Chemistry<TF>::exec_stats(const int iteration, const double time, Stats<TF>& stats)
{
    if (!sw_chemistry or stats.get_switch())
        return;

    const TF no_offset = 0.;
    const TF no_threshold = 0.;
    auto& gd = grid.get_grid_data();

    const TF NREACT = TF(1);
    if (iteration != 0)   // this does not make sense for first step = t=0.
    {
        // add deposition velocities to statistics:
        stats.calc_stats_2d("vdnh3"   , vdnh3,   no_offset);

        // Unit conversion constants
        const TF xmnh3 = 17.031;                        // Molar mass NH3 [g mol⁻¹]
        const TF xmnh3_i = TF(1.0) / xmnh3;               // [mol g⁻¹]
        const TF m2_to_ha = TF(1.0e4);                  // [m² ha⁻¹] 
        const TF s_to_year = TF(365.25 * 24 * 3600);    // [s yr⁻¹]
        // Combined conversion factor: [kg m⁻²] → [mol ha⁻¹ yr⁻¹]
        const TF conversion_factor = xmnh3_i * m2_to_ha * s_to_year * 1.0e3;

        // Convert flux_nh3 on-the-fly for statistics (keep original in kg m⁻²)
        std::vector<TF> flux_nh3_mol_ha_yr(gd.ijcells);
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ij = i + j*gd.jstride;
                // Convert periodic flux: [kg m⁻²] over period → [mol ha⁻¹ yr⁻¹]
                flux_nh3_mol_ha_yr[ij] = (flux_nh3[ij] / trfa) * conversion_factor;
            }

        // for (int j=gd.jstart; j<gd.jend; ++j)
        //     for (int i=gd.istart; i<gd.iend; ++i)
        //     {
        //         const int ij = i + j*gd.jstride;
        //         flux_nh3[ij] /= trfa;
        //     } 

        stats.calc_stats_2d("flux_nh3", flux_nh3_mol_ha_yr, no_offset); //added for nh3_flux
        stats.calc_stats_2d("flux_inst", flux_inst, no_offset); // added for instantaneous deposition flux of NH3

        // calculate total flux statistics (cumulative) in [mol ha⁻¹]
        std::vector<TF> total_flux_mol_ha(gd.ijcells);
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ij = i + j*gd.jstride;
                total_flux_mol_ha[ij] = total_flux_nh3[ij] * xmnh3_i * m2_to_ha * 1.0e3;  // [kg m⁻²] → [mol ha⁻¹]
            }
        
        stats.calc_stats_2d("c_extrap_diff", c_extrap_diff, no_offset);
        stats.calc_stats_2d("total_flux_mol_ha", total_flux_mol_ha, no_offset);  // Total [mol ha⁻¹]

        // Calculate statistics for new variables
        stats.calc_stats_2d("cstar1", cstar1, no_offset);
        stats.calc_stats_2d("cstar2", cstar2, no_offset);
        stats.calc_stats_2d("c_target", c_target, no_offset);
        stats.calc_stats_2d("c_diff_flux", c_diff_flux, no_offset);


        // Reset the periodic flux after saving to stats
        // Reset ONLY the periodic flux (NOT the total)
        trfa = 0;
        std::fill(flux_nh3.begin(), flux_nh3.end(), TF(0));
        // NOTE: total_flux_nh3 is NOT reset - it keeps accumulating

        // sum of all PEs:
        // printf("trfa: %13.4e iteration: %i time: %13.4e \n", trfa,iteration,time);
        master.sum(rfa.data(),NREACT*gd.ktot);
        // for (int l=0; l<NREACT*gd.ktot; ++l)
        //     rfa[l] /= (trfa*gd.itot*gd.jtot);    // mean over the horizontal plane in molecules/(cm3 * s)

        // Put the data into the NetCDF file.
        const std::vector<int> time_index{statistics_counter};

        // Write the time and iteration number.
        m.time_var->insert(time     , time_index);
        m.iter_var->insert(iteration, time_index);

        const std::vector<int> time_rfaz_index = {statistics_counter, 0};

        m.profs.at("chem_budget").data = rfa;

        const int ksize = NREACT*gd.ktot;
        std::vector<int> time_rfaz_size  = {1, ksize};
        std::vector<TF> prof_nogc(
                m.profs.at("chem_budget").data.begin() ,
                m.profs.at("chem_budget").data.begin() + ksize);

        m.profs.at("chem_budget").ncvar.insert(prof_nogc, time_rfaz_index, time_rfaz_size);

        // Synchronize the NetCDF file.
        m.data_file->sync();
        // Increment the statistics index.
        ++statistics_counter;

    }

    // (re-)intialize statistics
    // for (int l=0; l<NREACT*gd.ktot; ++l)
    //     rfa[l] = 0.0;
    // trfa = (TF) 0.0;
}


template <typename TF>
void Chemistry<TF>::init(Input& inputin)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    statistics_counter = 0;

    // Initialize existing arrays
    vdnh3.resize(gd.ijcells);
    flux_nh3.resize(gd.ijcells);
    std::fill(flux_nh3.begin(), flux_nh3.end(), TF(0));
    flux_inst.resize(gd.ijcells);
    std::fill(flux_inst.begin(), flux_inst.end(), TF(0));
    total_flux_nh3.resize(gd.ijcells);
    std::fill(total_flux_nh3.begin(), total_flux_nh3.end(), TF(0));

    // Initialize new arrays for concentration scaling with zeros
    cstar1.resize(gd.ijcells);
    std::fill(cstar1.begin(), cstar1.end(), TF(0));  // Initialize with zeros

    cstar2.resize(gd.ijcells);
    std::fill(cstar2.begin(), cstar2.end(), TF(0));  // Initialize with zeros
    
    // Only one concentration array needed now (optimal method)
    c_target.resize(gd.ijcells);
    std::fill(c_target.begin(), c_target.end(), TF(0));

<<<<<<< Updated upstream
    c_diff_flux.resize(gd.ijcells);
    std::fill(c_diff_flux.begin(), c_diff_flux.end(), TF(0));
=======
    c_extrap_diff.resize(gd.ijcells);
    std::fill(c_extrap_diff.begin(), c_extrap_diff.end(), TF(0));

    c_diff.resize(gd.ijcells);
    std::fill(c_diff.begin(), c_diff.end(), TF(0));
>>>>>>> Stashed changes

    // Initialize deposition routine
    deposition->init(inputin);

    // Fill deposition with standard values
    std::fill(vdnh3.begin(), vdnh3.end(), deposition->get_vd("nh3"));

    // master.print_message("Deposition arrays initialized, e.g. with vdnh3 = %13.5e m/s \n", deposition-> get_vd("nh3"));
}

// Fixed exec method to properly access boundary information


template <typename TF>
void Chemistry<TF>::create(
        const Timeloop<TF>& timeloop, std::string sim_name, Netcdf_handle& input_nc,
        Stats<TF>& stats, Cross<TF>& cross)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();
    int iotime = timeloop.get_iotime();

    Netcdf_group& group_nc = input_nc.get_group("timedep_chem");
    int time_dim_length;
    std::string time_dim;

    for (std::string varname : jname)    // check dimensions:
    {
        std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
        std::pair<std::string, int> unique_time = check_for_unique_time_dim(dims);
        time_dim = unique_time.first;
        time_dim_length = unique_time.second;
        time.resize(time_dim_length);
    }

    for (std::string varname : ename)    // check dimension also of emissions
    {
        std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
        std::pair<std::string, int> unique_time = check_for_unique_time_dim(dims);
        time_dim = unique_time.first;
        time_dim_length = unique_time.second;
        time.resize(time_dim_length);
    }

    jo31d.resize(time_dim_length);
    jh2o2.resize(time_dim_length);
    jno2.resize(time_dim_length);
    jno3.resize(time_dim_length);
    jn2o5.resize(time_dim_length);
    jch2or.resize(time_dim_length);
    jch2om.resize(time_dim_length);
    jch3o2h.resize(time_dim_length);
    emi_isop.resize(time_dim_length);
    emi_no.resize(time_dim_length);

    group_nc.get_variable(time, time_dim, {0}, {time_dim_length});
    group_nc.get_variable(jo31d, jname[0],  {0}, {time_dim_length});
    group_nc.get_variable(jh2o2, jname[1],  {0}, {time_dim_length});
    group_nc.get_variable(jno2, jname[2],  {0}, {time_dim_length});
    group_nc.get_variable(jno3, jname[3],  {0}, {time_dim_length});
    group_nc.get_variable(jn2o5, jname[4],  {0}, {time_dim_length});
    group_nc.get_variable(jch2or, jname[5],  {0}, {time_dim_length});
    group_nc.get_variable(jch2om, jname[6],  {0}, {time_dim_length});
    group_nc.get_variable(jch3o2h, jname[7],  {0}, {time_dim_length});
    group_nc.get_variable(emi_isop, ename[0],  {0}, {time_dim_length});
    group_nc.get_variable(emi_no,   ename[1],  {0}, {time_dim_length});

    // Store output of averaging.
    const TF NREACT = TF(1);
    //rfa.resize(NREACT*gd.ktot);
    // for (int l=0;l<NREACT*gd.ktot;++l)
    //     rfa[l] = 0.0;
    // trfa = (TF)0.0;
    qprof.resize(gd.kcells);
    tprof.resize(gd.kcells);

    if (stats.get_switch())
    {
        // Stats:
        const std::string group_name = "default";
        const std::vector<std::string> stat_op_def = {"mean", "2", "3", "4", "w", "grad", "diff", "flux", "path"};
        const std::vector<std::string> stat_op_w = {"mean", "2", "3", "4"};
        const std::vector<std::string> stat_op_p = {"mean", "2", "w", "grad"};

        std::stringstream filename;
        filename << sim_name << "." << "chemistry" << "." << std::setfill('0') << std::setw(7) << iotime << ".nc";

        // Create new NetCDF file in Mask<TF> m
        m.data_file = std::make_unique<Netcdf_file>(master, filename.str(), Netcdf_mode::Create);

        // Create dimensions.
        m.data_file->add_dimension("z", gd.kmax);
        m.data_file->add_dimension("zh", gd.kmax+1);
        m.data_file->add_dimension("rfaz", NREACT*gd.ktot);
        m.data_file->add_dimension("ijcells",gd.ijcells);
        m.data_file->add_dimension("time");

        // Create variables belonging to dimensions.
        Netcdf_handle& iter_handle =
            m.data_file->group_exists("default") ? m.data_file->get_group("default") : m.data_file->add_group("default");

        m.iter_var = std::make_unique<Netcdf_variable<int>>(iter_handle.add_variable<int>("iter", {"time"}));
        m.iter_var->add_attribute("units", "-");
        m.iter_var->add_attribute("long_name", "Iteration number");

        m.time_var = std::make_unique<Netcdf_variable<TF>>(m.data_file->template add_variable<TF>("time", {"time"}));
        if (timeloop.has_utc_time())
            m.time_var->add_attribute("units", "seconds since " + timeloop.get_datetime_utc_start_string());
        else
            m.time_var->add_attribute("units", "seconds since start");
        m.time_var->add_attribute("long_name", "Time");

        Netcdf_variable<TF> z_var = m.data_file->template add_variable<TF>("z", {"z"});
        z_var.add_attribute("units", "m");
        z_var.add_attribute("long_name", "Full level height");

        Netcdf_variable<TF> zh_var = m.data_file->template add_variable<TF>("zh", {"zh"});
        zh_var.add_attribute("units", "m");
        zh_var.add_attribute("long_name", "Half level height");

        std::string name = "chem_budget";
        std::string longname = "chemistry budget per layer";
        std::string unit = "molecules cm-3 s-1";
        Netcdf_variable<TF> rfaz_var = m.data_file->template add_variable<TF>("rfaz", {"rfaz"});
        rfaz_var.add_attribute("units", unit);
        rfaz_var.add_attribute("long_name", longname);

        // add a profile of reaction rates x z
        Level_type level =  Level_type::Full;

        Netcdf_handle& handle =
            m.data_file->group_exists("default") ? m.data_file->get_group("default") : m.data_file->add_group("default");
        Prof_var<TF> tmp{handle.add_variable<TF>(name, {"time", "rfaz"}), std::vector<TF>(gd.ktot*NREACT), level};
        m.profs.emplace(
                std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(std::move(tmp)));

        m.profs.at(name).ncvar.add_attribute("units", unit);
        m.profs.at(name).ncvar.add_attribute("long_name", longname);

        // Save the grid variables.
        std::vector<TF> z_nogc (gd.z. begin() + gd.kstart, gd.z. begin() + gd.kend  );
        std::vector<TF> zh_nogc(gd.zh.begin() + gd.kstart, gd.zh.begin() + gd.kend+1);
        z_var .insert( z_nogc, {0});
        zh_var.insert(zh_nogc, {0});

        // Synchronize the NetCDF file.
        m.data_file->sync();

        m.nmask. resize(gd.kcells);
        m.nmaskh.resize(gd.kcells);

        // add the deposition-velocity timeseries in deposition group statistics
        const std::string group_named = "deposition";

        // used in chemistry:
        stats.add_time_series("vdnh3", "NH3 deposition velocity", "m s-1", group_named);
        //stats.add_time_series("flux_nh3", "NH3 surface flux", "mol(NH3) m-2 s-1", group_named);
        stats.add_time_series("cstar1", "C*_1 concentration scaling parameter", "mol mol-1", group_named);
        stats.add_time_series("cstar2", "C*_2 concentration scaling parameter", "mol mol-1", group_named);
        stats.add_time_series("c_target", "NH3 concentration at target height (optimal)", "mol mol-1", group_named);
<<<<<<< Updated upstream
        stats.add_time_series("c_diff_flux", "Concentration difference flux (c_target - c_grid_closest) × 10^9 × rho × conversion", "kg m-2 s-1", group_named);
=======
        stats.add_time_series("c_extrap_diff", "Concentration difference between z_target and first grid level", "ppb", group_named);
        stats.add_time_series("c_diff", "Stability effect on extrapolated concentration (cstar1-cstar2)*factor", "ppb", group_named);
>>>>>>> Stashed changes
        stats.add_time_series("total_flux_mol_ha", "NH3 total cumulative flux", "mol ha-1", group_named);
    }

    // add cross-sections
    if (cross.get_switch())
    {
        //std::vector<std::string> allowed_crossvars = {"vdnh3"};
<<<<<<< Updated upstream
        std::vector<std::string> allowed_crossvars = {"vdnh3","flux_nh3","flux_inst","total_flux_mol_ha","cstar1","cstar2","c_grid_closest","c_target","c_diff_flux"};
=======
        std::vector<std::string> allowed_crossvars = {"vdnh3","flux_nh3","flux_inst","total_flux_mol_ha","cstar1","cstar2","c_target","c_extrap_diff","c_diff"};

>>>>>>> Stashed changes
        cross_list = cross.get_enabled_variables(allowed_crossvars);

        // `deposition->create()` only creates cross-sections.
        deposition->create(stats, cross);
    }
}

template<typename TF>
void Chemistry<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    const TF no_offset = TF(0);

    for (auto& name : cross_list)
    {
        if (name == "vdnh3")
            cross.cross_plane(vdnh3.data(), no_offset, name, iotime);
        else if (name == "flux_nh3") //added for nh3_flux
        {
            // Convert on-the-fly for cross-sections
            std::vector<TF> temp_flux_mol_ha_yr(gd.ijcells);
            const TF xmnh3_i = 1.0/17.031;
            const TF conversion_factor = xmnh3_i * 1.0e4 * (365.25 * 24 * 3600) * 1.0e3;
            for (int ij = 0; ij < gd.ijcells; ++ij)
                temp_flux_mol_ha_yr[ij] = (flux_nh3[ij] / trfa) * conversion_factor;
            cross.cross_plane(temp_flux_mol_ha_yr.data(), no_offset, name, iotime);
        }
        else if (name == "total_flux_mol_ha")
        {
            // Convert on-the-fly for cross-sections
            std::vector<TF> temp_mol_ha(gd.ijcells);
            const TF kg_to_mol_ha = (1.0/17.031) * 1.0e4 * 1.0e3;
            for (int ij = 0; ij < gd.ijcells; ++ij)
                temp_mol_ha[ij] = total_flux_nh3[ij] * kg_to_mol_ha;
            cross.cross_plane(temp_mol_ha.data(), no_offset, name, iotime);
        }
        else if (name == "flux_inst")  //added for instantaneous deposition flux of NH3
            cross.cross_plane(flux_inst.data(), no_offset, name, iotime);
        else if (name == "cstar1")
            cross.cross_plane(cstar1.data(), no_offset, name, iotime);
        else if (name == "cstar2")
            cross.cross_plane(cstar2.data(), no_offset, name, iotime);
        else if (name == "c_target")
            cross.cross_plane(c_target.data(), no_offset, name, iotime);
<<<<<<< Updated upstream
        else if (name == "c_diff_flux")
            cross.cross_plane(c_diff_flux.data(), no_offset, name, iotime);
=======
        else if (name == "c_extrap_diff")
            cross.cross_plane(c_extrap_diff.data(), no_offset, name, iotime);
        else if (name == "c_diff")
            cross.cross_plane(c_diff.data(), no_offset, name, iotime);
>>>>>>> Stashed changes
    }

    // see if to write per tile:
    deposition->exec_cross(cross, iotime);
}

template <typename TF>
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary, Thermo<TF>& thermo) // Added Thermo parameter
{
    if (!sw_chemistry)
        return;

    Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(time);

    jval[0] = ifac.fac0 * jo31d[ifac.index0] + ifac.fac1 * jo31d[ifac.index1];
    jval[1] = ifac.fac0 * jh2o2[ifac.index0] + ifac.fac1 * jh2o2[ifac.index1];
    jval[2] = ifac.fac0 * jno2[ifac.index0]  + ifac.fac1 * jno2[ifac.index1];
    jval[3] = ifac.fac0 * jno3[ifac.index0] + ifac.fac1 * jno3[ifac.index1];
    jval[4] = ifac.fac0 * jn2o5[ifac.index0] + ifac.fac1 * jn2o5[ifac.index1];
    jval[5] = ifac.fac0 * jch2or[ifac.index0] + ifac.fac1 * jch2or[ifac.index1];
    jval[6] = ifac.fac0 * jch2om[ifac.index0] + ifac.fac1 * jch2om[ifac.index1];
    jval[7] = ifac.fac0 * jch3o2h[ifac.index0] + ifac.fac1 * jch3o2h[ifac.index1];
    emval[0] = ifac.fac0 * emi_isop[ifac.index0] + ifac.fac1 * emi_isop[ifac.index1];
    emval[1] = ifac.fac0 * emi_no[ifac.index0] + ifac.fac1 * emi_no[ifac.index1];

    deposition->update_time_dependent(
            timeloop,
            boundary,
            thermo,
            //fields.sp.at("nh3")->fld.data(),  // Pass NH3 concentration
            vdnh3.data());
}

#ifndef USECUDA

template<typename TF>
void Chemistry<TF>::calc_c_target(Boundary<TF>& boundary)
{
    if (!sw_chemistry)
        return;
        
    auto& gd = grid.get_grid_data();
    // const TF target_height = z_target;
    
    // Constants for flux calculation (same as in pss function)
    const TF xmair = 28.9647;       // Molar mass of dry air  [kg kmol-1]
    const TF xmair_i = TF(1) / xmair;
    const TF xmnh3 = 17.031;
    const int k = gd.kstart;        // Surface level index

    // Get roughness length from boundary
    const std::vector<TF>& z0m = boundary.get_z0m();
    
    // Get Obukhov length
    const auto& tiles = boundary.get_tiles();
    std::vector<TF> obuk(gd.ijcells, TF(0));
    
    for (const auto& tile : tiles) {
        const auto& tile_data = tile.second;
        for (int ij = 0; ij < gd.ijcells; ++ij) {
            obuk[ij] += tile_data.fraction[ij] * tile_data.obuk[ij];
        }
    }
    
    // Calculate concentrations using FIXED method (existing cstar1)
    for (int j = gd.jstart; j < gd.jend; ++j) {
        for (int i = gd.istart; i < gd.iend; ++i) {
            const int ij = i + j * gd.jstride;
            // const int ijk1 = i + j * gd.jstride + gd.kstart * gd.ijcells;
            const int ijk1 = i + j * gd.jstride + gd.kstart * gd.kstride;

            // Calculate z_target for this grid point
            const TF z_target = rsl_ratio * z0m[ij];
            
            // Get concentration at first grid level
            const TF c_1 = fields.sp.at("nh3")->fld[ijk1];
            
            // Get heights
            const TF z_1 = gd.z[gd.kstart];
            
            // Use pre-calculated cstar1 from pss function
            const TF cstar_fixed = cstar1[ij];
            
            // Calculate scaling factor from z_1 to target height
            const TF scaling_factor = calc_factor(z_1, z_target, obuk[ij]);
            
            // Calculate c_grid_closest: Find closest grid point to target height
            int k_closest = gd.kstart;
            TF min_distance = std::abs(gd.z[gd.kstart] - z_target);
            
            for (int k_grid = gd.kstart; k_grid < gd.kend; ++k_grid) {
                TF distance = std::abs(gd.z[k_grid] - z_target);
                if (distance < min_distance) {
                    min_distance = distance;
                    k_closest = k_grid;
                }
            }
            
            // Get actual simulated concentration at closest grid point
            //const int ijk_closest = i + j * gd.jstride + k_closest * gd.ijcells;
            const int ijk_closest = i + j * gd.jstride + k_closest * gd.kstride;

            c_grid_closest[ij] = fields.sp.at("nh3")->fld[ijk_closest];

            const TF concentration_diff = c_target[ij] - c_grid_closest[ij];
            c_diff_flux[ij] = concentration_diff * TF(1e9) * fields.rhoref[k] * xmair_i * xmnh3;
        }
    }
}

// what the following function (exec) does:
//     Gets atmospheric conditions (temperature, humidity, wind)
//     Calculates chemical reactions (like NH₃ decay)
//     Computes surface fluxes (how much NH₃ deposits to ground)
//     Updates concentrations for the next timestep
//     Calls our new calc_c_target() to get target height concentrations
template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo, Boundary<TF>& boundary, double sdt, double dt) // Added Boundary parameter to access Obukhov length and z0m
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    auto tmp = fields.get_tmp();
    thermo.get_thermo_field(*tmp, "T", true, false);

    // Calculate the mean temperature and water vapor mixing ratio
    field3d_operators.calc_mean_profile(tprof.data(), tmp->fld.data());
    qprof = fields.sp.at("qt")->fld_mean;
    
    // Get roughness length and Obukhov length from boundary
    const std::vector<TF>& z0m = boundary.get_z0m();
    const auto& tiles = boundary.get_tiles();
    std::vector<TF> obuk(gd.ijcells, TF(0));
    
    for (const auto& tile : tiles) {
        const auto& tile_data = tile.second;
        for (int ij = 0; ij < gd.ijcells; ++ij) {
            obuk[ij] += tile_data.fraction[ij] * tile_data.obuk[ij];
        }
    }

    pss<TF>(
            fields.st.at("nh3")->fld.data(),
            fields.sp.at("nh3")->fld.data(),
            jval, emval,
            vdnh3.data(),
            tprof.data(),
            qprof.data(),
            gd.dzi.data(),
            fields.rhoref.data(),
            gd.z.data(),                   
            obuk.data(),                   // CHANGE TO obuk.data()
            z0m.data(),                    
            rfa.data(),
            flux_nh3.data(),
            flux_inst.data(),
            total_flux_nh3.data(), 
            cstar1.data(),               
            cstar2.data(),
            c_target.data(),
<<<<<<< Updated upstream
=======
            c_extrap_diff.data(),
            c_diff.data(),
>>>>>>> Stashed changes
            trfa,
            dt, sdt, lifetime,
            // z_target,
            rsl_ratio,
<<<<<<< Updated upstream
=======
            sw_const_ref_height,
            z_fixed,
            sw_adapt_ref_height,
            sw_use_lowest_levels,
>>>>>>> Stashed changes
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells,
            gd.dx, gd.dy);

    calc_c_target(boundary);


    fields.release_tmp(tmp);
}

#endif

template class Chemistry<double>;
//:template class Chemistry<float>;



