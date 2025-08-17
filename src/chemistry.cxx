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

#include "boundary.h"
#include "boundary_surface.h"
#include "boundary_surface_bulk.h"
#include "boundary_surface_lsm.h"
#include "chemistry.h"
#include "constants.h"
#include "cross.h"
#include "deposition.h"
#include "fields.h"
#include "grid.h"
#include "master.h"
#include "monin_obukhov.h"
#include "netcdf_interface.h"
#include "stats.h"
#include "thermo.h"
#include "timeloop.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <utility>
#include <vector> 

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
        
        const TF z1_over_L = z1 / L;
        const TF z2_over_L = z2 / L;
        
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
    
    /**
     * APPROACH B: Use existing cstar1 with optimal reference
     */
    template<typename TF>
    TF calc_approach_B(
        const TF cstar1_val,
        const TF* const nh3,
        const TF* const z,
        const TF z_target,
        const TF L,
        const int i, const int j,
        const int kstart, const int kend,
        const int jstride, const int kstride)
    {
        // Find optimal reference height
        const int k1 = find_ref_below(z, z_target, kstart, kend);
        
        // Get reference concentration and height
        const int ijk1 = i + j * jstride + k1 * kstride;
        const TF c1 = nh3[ijk1];
        const TF z1 = z[k1];
        
        // Calculate factor using general function
        const TF scaling_factor = calc_factor(z1, z_target, L);      // For scaling to target height
        
        return c1 + cstar1_val * scaling_factor;
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
        TF* restrict cstar,
        TF* restrict c_target, 
        TF& trfa,
        const TF dt,
        const TF sdt,
        const TF lifetime,
        const TF* const restrict z_ref_field,  // ADD: adaptive reference heights
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
        //if (abs(sdt/dt - 1./3.) < 1e-5) trfa += dt;
        trfa += sdt;
    
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
                        
                        // Calculate cstar with stability correction (remove cstar2)
                        if (std::abs(gradient_factor) > TF(1e-15))
                        {
                            cstar[ij] = +1.0 * (c_2 - c_1) / gradient_factor;
                        }
                        else
                        {
                            cstar[ij] = TF(0.0);
                        }
                        
                        // Use adaptive reference height for each grid point
                        const TF local_z_ref = z_ref_field[ij];
                        const TF scaling_factor = calc_factor(z_1, local_z_ref, L);
                        c_target[ij] = c_1 + cstar[ij] * scaling_factor;
    
                        // Calculate and accumulate flux for this RK3 step
                        // Note: flux is accumulated (+=) and scaled by sdt
    
                        // // Calculate instantaneous flux first (original method)
                        // flux_inst[ij] = (-1.0) * vdnh3[ij] * nh3[ijk] * rhoref[k] * xmair_i * xmnh3; // [kg(NH3) m-2 s-1]
                        flux_inst[ij] = (-1.0) * vdnh3[ij] * c_target[ij] * rhoref[k] * xmair_i * xmnh3; // [kg(NH3) m-2 s-1]
    
                        // accumulate over sub-timestep for statistics (gets reset periodically)
                        TF flux = flux_inst[ij] * sdt; // [kg m⁻² s⁻¹] × [s] = [kg m⁻²] Scale by timestep for accumulation     
                        flux_nh3[ij] += flux;        // For period statistics - Accumulate [kg m⁻²]
                        
                        // accumulate total flux (never gets reset)
                        total_flux_nh3[ij] += flux;  // [kg m⁻²]
    
                        // decay = vdnh3[ij]*dzi[k] + lti;   // 1/s
                        if (std::abs(nh3[ijk]) > TF(1e-15)) //to prevent division by zero
                        {
                            decay = (vdnh3[ij] * dzi[k] * c_target[ij] / nh3[ijk]) + lti;   // 1/s
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
Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Radiation<TF>& radiationin, Input& inputin):
    master(masterin), grid(gridin), fields(fieldsin), radiation(radiationin), field3d_operators(master, grid, fields)
{
    // Rest of the constructor remains the same
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();
    sw_chemistry = inputin.get_item<bool>("chemistry", "swchemistry", "", false);
    //lifetime     = inputin.get_item<TF>("chemistry", "lifetime", "", (TF)72000);  // seconds (20 hour default)
    lifetime     = inputin.get_item<TF>("chemistry", "lifetime", "", (TF)1e30);  // seconds

    master.print_message("Lifetime of the tracer:  = %13.5e s \n", lifetime);
    if (!sw_chemistry)
        return;
    //deposition = std::make_shared<Deposition <TF>>(masterin, gridin, fieldsin, radiationin, inputin);

    deposition = std::make_shared<Deposition<TF>>(master, grid, fields, radiation, *this, inputin);
    //The *this passes the current Chemistry object as a reference to the Deposition constructor.

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
        
        stats.calc_stats_2d("total_flux_mol_ha", total_flux_mol_ha, no_offset);  // Total [mol ha⁻¹]

        // Calculate statistics for new variables
        stats.calc_stats_2d("cstar", cstar, no_offset);
        stats.calc_stats_2d("c_ref_grid", c_ref_grid, no_offset);
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
    cstar.resize(gd.ijcells);
    std::fill(cstar.begin(), cstar.end(), TF(0));
    
    c_ref_grid.resize(gd.ijcells);
    std::fill(c_ref_grid.begin(), c_ref_grid.end(), TF(0));

    // Only one concentration array needed now (optimal method)
    c_target.resize(gd.ijcells);
    std::fill(c_target.begin(), c_target.end(), TF(0));

    c_diff_flux.resize(gd.ijcells);
    std::fill(c_diff_flux.begin(), c_diff_flux.end(), TF(0));

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
        stats.add_time_series("cstar", "C* concentration scaling parameter with stability correction", "mol mol-1", group_named);
        stats.add_time_series("c_ref_grid", "NH3 concentration at closest grid point to adaptive reference height", "mol mol-1", group_named);
        stats.add_time_series("c_target", "NH3 concentration at adaptive reference height (optimal)", "mol mol-1", group_named);
        stats.add_time_series("c_diff_flux", "Concentration difference flux (c_target - c_ref_grid) × 10^9 × rho × conversion", "kg m-2 s-1", group_named);
        stats.add_time_series("total_flux_mol_ha", "NH3 total cumulative flux", "mol ha-1", group_named);
    }

    // add cross-sections
    if (cross.get_switch())
    {
        //std::vector<std::string> allowed_crossvars = {"vdnh3"};
        std::vector<std::string> allowed_crossvars = {"vdnh3","flux_nh3","flux_inst","total_flux_mol_ha","cstar","c_ref_grid","c_target","c_diff_flux"};
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
        else if (name == "cstar")
            cross.cross_plane(cstar.data(), no_offset, name, iotime);
        else if (name == "c_ref_grid")
            cross.cross_plane(c_ref_grid.data(), no_offset, name, iotime);
        else if (name == "c_target")
            cross.cross_plane(c_target.data(), no_offset, name, iotime);
        else if (name == "c_diff_flux")
            cross.cross_plane(c_diff_flux.data(), no_offset, name, iotime);
    }

    // see if to write per tile:
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
            //fields.sp.at("nh3")->fld.data(),  // Pass NH3 concentration
            vdnh3.data());
}

#ifndef USECUDA

template<typename TF>
void Chemistry<TF>::calc_c_ref_grid(Boundary<TF>& boundary, const TF* const z_ref_field, const TF* const z_ref_level_field)
{
    if (!sw_chemistry)
        return;
        
    auto& gd = grid.get_grid_data();
    
    // Constants for flux calculation
    const TF xmair = 28.9647;
    const TF xmair_i = TF(1) / xmair;
    const TF xmnh3 = 17.031;
    const int k = gd.kstart;
    
    for (int j = gd.jstart; j < gd.jend; ++j) {
        for (int i = gd.istart; i < gd.iend; ++i) {
            const int ij = i + j * gd.jstride;
            
            // Use pre-calculated grid level index (no search needed!)
            const int k_ref = static_cast<int>(z_ref_level_field[ij]);
            
            // Get concentration at the exact reference level
            const int ijk_ref = i + j * gd.jstride + k_ref * gd.kstride;
            c_ref_grid[ij] = fields.sp.at("nh3")->fld[ijk_ref];

            // Calculate difference between scaled and grid concentrations
            const TF concentration_diff = c_target[ij] - c_ref_grid[ij];
            c_diff_flux[ij] = concentration_diff * TF(1e9) * fields.rhoref[k] * xmair_i * xmnh3;
        }
    }
}

// what the following function (exec) does:
//     Gets atmospheric conditions (temperature, humidity, wind)
//     Calculates chemical reactions (like NH₃ decay)
//     Computes surface fluxes (how much NH₃ deposits to ground)
//     Updates concentrations for the next timestep
//     Calls our new calc_c20m() to get 20m concentrations
template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo, Boundary<TF>& boundary, double sdt, double dt)
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


// Get adaptive reference heights from boundary layer
// Cast to specific boundary type to access getter methods
const std::vector<TF>* z_ref_field_ptr = nullptr;
const std::vector<TF>* z_ref_level_field_ptr = nullptr;

if (auto* lsm_boundary = dynamic_cast<Boundary_surface_lsm<TF>*>(&boundary)) {
    z_ref_field_ptr = &(lsm_boundary->get_z_ref_field());
    z_ref_level_field_ptr = &(lsm_boundary->get_z_ref_level_field());
}
else if (auto* surface_boundary = dynamic_cast<Boundary_surface<TF>*>(&boundary)) {
    z_ref_field_ptr = &(surface_boundary->get_z_ref_field());
    z_ref_level_field_ptr = &(surface_boundary->get_z_ref_level_field());
}
else if (auto* bulk_boundary = dynamic_cast<Boundary_surface_bulk<TF>*>(&boundary)) {
    z_ref_field_ptr = &(bulk_boundary->get_z_ref_field());
    z_ref_level_field_ptr = &(bulk_boundary->get_z_ref_level_field());
}
else {
    throw std::runtime_error("Unknown boundary type in chemistry exec()");
}

const std::vector<TF>& z_ref_field = *z_ref_field_ptr;
const std::vector<TF>& z_ref_level_field = *z_ref_level_field_ptr;

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
            obuk.data(),
            z0m.data(),                    
            rfa.data(),
            flux_nh3.data(),
            flux_inst.data(),
            total_flux_nh3.data(), 
            cstar.data(),
            c_target.data(),
            trfa,
            dt, sdt, lifetime,
            z_ref_field.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells,
            gd.dx, gd.dy);
    
    calc_c_ref_grid(boundary, z_ref_field.data(), z_ref_level_field.data());

    fields.release_tmp(tmp);
}

#endif

template class Chemistry<double>;
//:template class Chemistry<float>;



