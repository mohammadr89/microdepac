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
 * c_extrap_diff : Concentration difference between z_target and first grid level [ppb]
 * rsl_ratio    = float   : Roughness sublayer ratio (default: 20.0, only used if sw_const_ref_height=false)
 * sw_const_ref_height = boolean : Use constant reference height (default: true)
 * z_fixed            = float   : Fixed reference height [m] (default: 20, used if sw_const_ref_height=true)
 * sw_adapt_ref_height = boolean : Use adaptive extrapolation method (default: true)
 * sw_use_lowest_levels = boolean : Use kstart/kstart+1 for cstar (true) or bracketing levels (false) (default: true)
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
 * chem_budget        : Chemistry budget per layer [molecules cm-3 s-1]
 * T_target           : Temperature at target height [K]
 * rho_target         : Air density at target height [kg m-3]
 */

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

template<typename TF>
inline TF calc_factor(const TF z1, const TF z2, const TF L)
{
    TF psi1 = TF(0), psi2 = TF(0);
    
    if (std::abs(L) > TF(1e-15))
    {
        namespace most = Monin_obukhov;
        const TF z1_over_L_raw = z1 / L;
        const TF z2_over_L_raw = z2 / L;
        const TF z1_over_L = std::min(std::max(z1_over_L_raw, Constants::zL_min<TF>), Constants::zL_max<TF>);
        const TF z2_over_L = std::min(std::max(z2_over_L_raw, Constants::zL_min<TF>), Constants::zL_max<TF>);
        psi1 = (z1_over_L <= TF(0)) ? most::psih_unstable(z1_over_L) : most::psih_stable(z1_over_L);
        psi2 = (z2_over_L <= TF(0)) ? most::psih_unstable(z2_over_L) : most::psih_stable(z2_over_L);
    }
    return (std::log(z2/z1) - psi2 + psi1) / Constants::kappa<TF>;
}

namespace
{
    double CFACTOR;
    std::pair<std::string, int> check_for_unique_time_dim(const std::map<std::string, int>& dims)
    {
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

    template<typename TF>
    inline int find_ref_below(const TF* const z, const TF z_target, const int kstart, const int kend)
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
    void pss
        (
        TF* restrict tnh3,
        const TF* const restrict nh3,
        const TF* const restrict jval,
        const TF* const restrict emval,
        const TF* const restrict vdnh3,
        const TF* const restrict tprof,
        const TF* const restrict qprof,
        const TF* const restrict dzi,
        const TF* const restrict rhoref,
        const TF* const restrict z,
        const TF* const restrict obuk,
        const TF* const restrict z0m,
        TF* restrict rfa,
        TF* restrict flux_nh3, 
        TF* restrict flux_inst, 
        TF* restrict total_flux_nh3,
        TF* restrict cstar1,  
        TF* restrict cstar2,  
        TF* restrict c_target,
        TF* restrict c_extrap_diff,
        TF* restrict c_diff,
        TF* restrict T_target,
        TF* restrict rho_target,
        const TF* const restrict T_fld,
        TF& trfa,
        const TF dt,
        const TF sdt,
        const TF lifetime,
        const TF rsl_ratio,
        const bool sw_const_ref_height,
        const TF z_fixed,              
        const bool sw_adapt_ref_height,
        const bool sw_use_lowest_levels,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int kstart, const int kend,
        const int jstride, const int kstride,
        const TF dx,
        const TF dy
        )

    {
        const TF xmh2o = 18.015265;
        const TF xmnh3 = 17.031;
        const TF xmh2o_i = TF(1) / xmh2o;
        const TF xmair = 28.9647;
        const TF xmair_i = TF(1) / xmair;
        const TF Na = 6.02214086e23;
    
        // Update the time integration of the reaction fluxes with the full timestep on first RK3 step
        if (abs(sdt/dt - 1./3.) < 1e-5) trfa += dt;
    
        namespace most = Monin_obukhov;

        for (int k=kstart; k<kend; ++k)
        {
            const TF C_M = TF(1e-3) * rhoref[k] * Na * xmair_i;
            const TF CFACTOR = C_M;
            const TF sdt_cfac_i = TF(1) / (sdt * CFACTOR);
            const TF lti = TF(1)/lifetime;
            TF decay;
            for (int j=jstart; j<jend; ++j)
            #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    const int ij = i + j*jstride;
    
                    if (k==kstart)
                    {
                        const int ijk1 = i + j*jstride + kstart*kstride;
                        const TF c_1 = nh3[ijk1];
                        const TF z_1 = z[kstart];
                        
                        // ========================================
                        // STEP 1: Determine reference height
                        // ========================================
                        TF z_target;
                        if (sw_const_ref_height)
                        {
                            z_target = z_fixed;
                        }
                        else
                        {
                            z_target = rsl_ratio * z0m[ij];
                        }
                        
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
                            T_target[ij]      = T_fld[i + j*jstride + kstart*kstride];
                            rho_target[ij]    = rhoref[kstart];

                            // Calculate instantaneous flux using first grid level concentration (SIMPLIFIED)
                            flux_inst[ij] = (-1.0) * vdnh3[ij] * nh3[ijk] * rhoref[k] * xmair_i * xmnh3; // [kg(NH3) m-2 s-1]
                            decay = vdnh3[ij]*dzi[k] + lti;
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
                                T_target[ij]      = T_fld[i + j*jstride + k_exact*kstride];
                                rho_target[ij]    = rhoref[k_exact];
                                const TF rhoref_exact = rhoref[k_exact];
                                flux_inst[ij] = (-1.0) * vdnh3[ij] * c_target[ij] * rhoref_exact * xmair_i * xmnh3;
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
                                const TF scaling_factor_neutral = calc_factor(z[k_below], z_target, TF(1e30));
                                c_diff[ij] = (cstar1[ij] - cstar2[ij]) * scaling_factor_neutral * TF(1e9);
                                
                                // Calculate flux using c_target
                                const TF rhoref_target = rhoref[k_below]
                                    + (z_target - z[k_below]) / (z[k_above] - z[k_below])
                                    * (rhoref[k_above] - rhoref[k_below]);
                                rho_target[ij] = rhoref_target;
                                T_target[ij] = T_fld[i + j*jstride + k_below*kstride]
                                    + (z_target - z[k_below]) / (z[k_above] - z[k_below])
                                    * (T_fld[i + j*jstride + k_above*kstride] - T_fld[i + j*jstride + k_below*kstride]);
                                flux_inst[ij] = (-1.0) * vdnh3[ij] * c_target[ij] * rhoref_target * xmair_i * xmnh3;
                                
                                // Calculate decay with concentration scaling
                                if (std::abs(nh3[ijk]) > TF(1e-15))
                                    decay = (vdnh3[ij] * dzi[k] * c_target[ij] / nh3[ijk]) + lti;
                                else
                                    decay = lti;
                            }
                        }
                        // Accumulate fluxes 
                        TF flux = flux_inst[ij] * sdt;
                        flux_nh3[ij] += flux;
                        total_flux_nh3[ij] += flux;
                    }
                    else
                    { 
                        decay = lti;
                    }
                    // Update tendencies
                    tnh3[ijk] -= decay*nh3[ijk];
                }
        }
    }
}

template<typename TF>
Chemistry<TF>::Chemistry
    (Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Radiation<TF>& radiationin, Input& inputin) : 
    master(masterin), grid(gridin), fields(fieldsin), radiation(radiationin), field3d_operators(master, grid, fields)
{
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();
    sw_chemistry = inputin.get_item<bool>("chemistry", "swchemistry", "", false);
    lifetime     = inputin.get_item<TF>("chemistry", "lifetime", "", (TF)1e30);

    // Roughness sublayer ratio (Basu & Lacser 2017, DOI: 10.1007/s10546-016-0225-y)
    rsl_ratio = inputin.get_item<TF>("chemistry", "rsl_ratio", "", (TF)20.0);
    sw_const_ref_height = inputin.get_item<bool>("chemistry", "sw_const_ref_height", "", true);
    z_fixed = inputin.get_item<TF>("chemistry", "z_fixed", "", (TF)20);
    sw_adapt_ref_height = inputin.get_item<bool>("chemistry", "sw_adapt_ref_height", "", true);
    sw_use_lowest_levels = inputin.get_item<bool>("chemistry", "sw_use_lowest_levels", "", true);

    if (!sw_chemistry)
        return;

    deposition = std::make_shared<Deposition<TF>>(masterin, gridin, fieldsin, radiationin, *this, inputin);
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
    if (iteration != 0)
    {
        stats.calc_stats_2d("vdnh3"   , vdnh3,   no_offset);
        const TF xmnh3 = 17.031;
        const TF xmnh3_i = TF(1.0) / xmnh3;
        const TF m2_to_ha = TF(1.0e4);
        const TF s_to_year = TF(365.25 * 24 * 3600);
        // Combined conversion factor: [kg m⁻²] → [mol ha⁻¹ yr⁻¹]
        const TF conversion_factor = xmnh3_i * m2_to_ha * s_to_year * 1.0e3;

        std::vector<TF> flux_nh3_mol_ha_yr(gd.ijcells);
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ij = i + j*gd.jstride;
                flux_nh3_mol_ha_yr[ij] = (flux_nh3[ij] / trfa) * conversion_factor;
            }

        stats.calc_stats_2d("flux_nh3", flux_nh3_mol_ha_yr, no_offset);
        stats.calc_stats_2d("flux_inst", flux_inst, no_offset);

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
        stats.calc_stats_2d("cstar1", cstar1, no_offset);
        stats.calc_stats_2d("cstar2", cstar2, no_offset);
        stats.calc_stats_2d("c_target", c_target, no_offset);
        stats.calc_stats_2d("c_diff", c_diff, no_offset);
        stats.calc_stats_2d("T_target", T_target, no_offset);
        stats.calc_stats_2d("rho_target", rho_target, no_offset);

        // Reset the periodic flux after saving to stats
        trfa = 0;
        std::fill(flux_nh3.begin(), flux_nh3.end(), TF(0));
        const std::vector<int> time_index{statistics_counter};
        m.time_var->insert(time     , time_index);
        m.iter_var->insert(iteration, time_index);
        const std::vector<int> time_rfaz_index = {statistics_counter, 0};
        m.profs.at("chem_budget").data = rfa;
        const int ksize = NREACT*gd.ktot;
        std::vector<int> time_rfaz_size  = {1, ksize};
        std::vector<TF> prof_nogc (m.profs.at("chem_budget").data.begin() ,m.profs.at("chem_budget").data.begin() + ksize);
        m.profs.at("chem_budget").ncvar.insert(prof_nogc, time_rfaz_index, time_rfaz_size);
        m.data_file->sync();
        ++statistics_counter;
    }
}

template <typename TF>
void Chemistry<TF>::init(Input& inputin)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();
    statistics_counter = 0;
    vdnh3.resize(gd.ijcells);
    flux_nh3.resize(gd.ijcells);
    std::fill(flux_nh3.begin(), flux_nh3.end(), TF(0));
    flux_inst.resize(gd.ijcells);
    std::fill(flux_inst.begin(), flux_inst.end(), TF(0));
    total_flux_nh3.resize(gd.ijcells);
    std::fill(total_flux_nh3.begin(), total_flux_nh3.end(), TF(0));
    cstar1.resize(gd.ijcells);
    std::fill(cstar1.begin(), cstar1.end(), TF(0));
    cstar2.resize(gd.ijcells);
    std::fill(cstar2.begin(), cstar2.end(), TF(0));
    c_target.resize(gd.ijcells);
    std::fill(c_target.begin(), c_target.end(), TF(0));
    c_extrap_diff.resize(gd.ijcells);
    std::fill(c_extrap_diff.begin(), c_extrap_diff.end(), TF(0));
    c_diff.resize(gd.ijcells);
    std::fill(c_diff.begin(), c_diff.end(), TF(0));
    T_target.resize(gd.ijcells);
    std::fill(T_target.begin(), T_target.end(), TF(0));
    rho_target.resize(gd.ijcells);
    std::fill(rho_target.begin(), rho_target.end(), TF(0));
    deposition->init(inputin);
    std::fill(vdnh3.begin(), vdnh3.end(), deposition->get_vd("nh3"));
}

template <typename TF>
void Chemistry<TF>::create
    (
    const Timeloop<TF>& timeloop, std::string sim_name, Netcdf_handle& input_nc,
    Stats<TF>& stats, Cross<TF>& cross
    )
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();
    int iotime = timeloop.get_iotime();

    Netcdf_group& group_nc = input_nc.get_group("timedep_chem");
    int time_dim_length;
    std::string time_dim;

    for (std::string varname : jname)
    {
        std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
        std::pair<std::string, int> unique_time = check_for_unique_time_dim(dims);
        time_dim = unique_time.first;
        time_dim_length = unique_time.second;
        time.resize(time_dim_length);
    }

    for (std::string varname : ename)
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

    const TF NREACT = TF(1);
    qprof.resize(gd.kcells);
    tprof.resize(gd.kcells);

    if (stats.get_switch())
    {
        const std::string group_name = "default";
        const std::vector<std::string> stat_op_def = {"mean", "2", "3", "4", "w", "grad", "diff", "flux", "path"};
        const std::vector<std::string> stat_op_w = {"mean", "2", "3", "4"};
        const std::vector<std::string> stat_op_p = {"mean", "2", "w", "grad"};
        std::stringstream filename;
        filename << sim_name << "." << "chemistry" << "." << std::setfill('0') << std::setw(7) << iotime << ".nc";
        m.data_file = std::make_unique<Netcdf_file>(master, filename.str(), Netcdf_mode::Create);
        m.data_file->add_dimension("z", gd.kmax);
        m.data_file->add_dimension("zh", gd.kmax+1);
        m.data_file->add_dimension("rfaz", NREACT*gd.ktot);
        m.data_file->add_dimension("ijcells",gd.ijcells);
        m.data_file->add_dimension("time");
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
        Level_type level =  Level_type::Full;
        Netcdf_handle& handle =
            m.data_file->group_exists("default") ? m.data_file->get_group("default") : m.data_file->add_group("default");
        Prof_var<TF> tmp{handle.add_variable<TF>(name, {"time", "rfaz"}), std::vector<TF>(gd.ktot*NREACT), level};
        m.profs.emplace(
                std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(std::move(tmp)));
        m.profs.at(name).ncvar.add_attribute("units", unit);
        m.profs.at(name).ncvar.add_attribute("long_name", longname);
        std::vector<TF> z_nogc (gd.z. begin() + gd.kstart, gd.z. begin() + gd.kend  );
        std::vector<TF> zh_nogc(gd.zh.begin() + gd.kstart, gd.zh.begin() + gd.kend+1);
        z_var .insert( z_nogc, {0});
        zh_var.insert(zh_nogc, {0});
        m.data_file->sync();
        m.nmask. resize(gd.kcells);
        m.nmaskh.resize(gd.kcells);
        const std::string group_named = "deposition";
        stats.add_time_series("vdnh3", "NH3 deposition velocity", "m s-1", group_named);
        stats.add_time_series("cstar1", "C*_1 concentration scaling parameter", "mol mol-1", group_named);
        stats.add_time_series("cstar2", "C*_2 concentration scaling parameter", "mol mol-1", group_named);
        stats.add_time_series("c_target", "NH3 concentration at target height (optimal)", "mol mol-1", group_named);
        stats.add_time_series("c_extrap_diff", "Concentration difference between z_target and first grid level", "ppb", group_named);
        stats.add_time_series("c_diff", "Stability effect on extrapolated concentration (cstar1-cstar2)*factor", "ppb", group_named);
        stats.add_time_series("total_flux_mol_ha", "NH3 total cumulative flux", "mol ha-1", group_named);
        stats.add_time_series("T_target", "Temperature at target height", "K", group_named);
        stats.add_time_series("rho_target", "Air density at target height", "kg m-3", group_named);
    }

    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = {"vdnh3","flux_nh3","flux_inst","total_flux_mol_ha","cstar1","cstar2","c_target","c_extrap_diff","c_diff","T_target","rho_target"};
        cross_list = cross.get_enabled_variables(allowed_crossvars);
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
        else if (name == "flux_nh3")
        {
            std::vector<TF> temp_flux_mol_ha_yr(gd.ijcells);
            const TF xmnh3_i = 1.0/17.031;
            const TF conversion_factor = xmnh3_i * 1.0e4 * (365.25 * 24 * 3600) * 1.0e3;
            for (int ij = 0; ij < gd.ijcells; ++ij)
                temp_flux_mol_ha_yr[ij] = (flux_nh3[ij] / trfa) * conversion_factor;
            cross.cross_plane(temp_flux_mol_ha_yr.data(), no_offset, name, iotime);
        }
        else if (name == "total_flux_mol_ha")
        {
            std::vector<TF> temp_mol_ha(gd.ijcells);
            const TF kg_to_mol_ha = (1.0/17.031) * 1.0e4 * 1.0e3;
            for (int ij = 0; ij < gd.ijcells; ++ij)
                temp_mol_ha[ij] = total_flux_nh3[ij] * kg_to_mol_ha;
            cross.cross_plane(temp_mol_ha.data(), no_offset, name, iotime);
        }
        else if (name == "flux_inst")
            cross.cross_plane(flux_inst.data(), no_offset, name, iotime);
        else if (name == "cstar1")
            cross.cross_plane(cstar1.data(), no_offset, name, iotime);
        else if (name == "cstar2")
            cross.cross_plane(cstar2.data(), no_offset, name, iotime);
        else if (name == "c_target")
            cross.cross_plane(c_target.data(), no_offset, name, iotime);
        else if (name == "c_extrap_diff")
            cross.cross_plane(c_extrap_diff.data(), no_offset, name, iotime);
        else if (name == "c_diff")
            cross.cross_plane(c_diff.data(), no_offset, name, iotime);
        else if (name == "T_target")
            cross.cross_plane(T_target.data(), no_offset, name, iotime);
        else if (name == "rho_target")
            cross.cross_plane(rho_target.data(), no_offset, name, iotime);
    }

    deposition->exec_cross(cross, iotime);
}

template <typename TF>
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary, Thermo<TF>& thermo)
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
            vdnh3.data());
}

#ifndef USECUDA

template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo, Boundary<TF>& boundary, double sdt, double dt)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();
    auto tmp = fields.get_tmp();
    thermo.get_thermo_field(*tmp, "T", true, false);
    field3d_operators.calc_mean_profile(tprof.data(), tmp->fld.data());
    qprof = fields.sp.at("qt")->fld_mean;
    const std::vector<TF>& z0m = boundary.get_z0m();
    const auto& tiles = boundary.get_tiles();
    std::vector<TF> obuk(gd.ijcells, TF(0));
    
    for (const auto& tile : tiles) {
        const auto& tile_data = tile.second;
        for (int ij = 0; ij < gd.ijcells; ++ij) {
            obuk[ij] += tile_data.fraction[ij] * tile_data.obuk[ij];
        }
    }

    pss<TF>
        (
        fields.st.at("nh3")->fld.data(),
        fields.sp.at("nh3")->fld.data(),
        jval, emval,
        vdnh3.data(),
        tprof.data(),
        qprof.data(),
        gd.dzi.data(),
        fields.rhoref.data(),
        gd.z.data(),                   
        obuk.data(),
        z0m.data(),                    
        rfa.data(),
        flux_nh3.data(),
        flux_inst.data(),
        total_flux_nh3.data(), 
        cstar1.data(),               
        cstar2.data(),
        c_target.data(),
        c_extrap_diff.data(),
        c_diff.data(),
        T_target.data(),
        rho_target.data(),
        tmp->fld.data(),
        trfa,
        dt, sdt, lifetime,
        rsl_ratio,
        sw_const_ref_height,
        z_fixed,
        sw_adapt_ref_height,
        sw_use_lowest_levels,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells,
        gd.dx, gd.dy
        );

    fields.release_tmp(tmp);
}

#endif

template class Chemistry<double>;
