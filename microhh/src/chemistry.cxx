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
                TF* restrict rfa,
                TF* restrict flux_nh3,  // New parameter for flux
                TF* restrict flux_inst, // Add this line for instantaneous flux
                TF& trfa,
                const TF dt,
                const TF sdt,
                const TF lifetime,
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
                                    // Calculate and accumulate flux for this RK3 step
                                    // Note: flux is accumulated (+=) and scaled by sdt

                                    //TF flux = (-1.0) * vdnh3[ij] * nh3[ijk] * rhoref[k] * xmair_i * xmnh3 * sdt; // [kg(NH3) m-2 s-1]
                                    //TF flux[ij] = (-1.0) * 1.0e3 * vdnh3[ij] * nh3[ijk] * rhoref[k] * xmair_i * sdt; // [mol(NH3) m-2 s-1 * sdt!!!]
                                    //flux_nh3[ij] = (-1.0) * 1.0e3 * vdnh3[ij] * nh3[ijk] * rhoref[k] * xmair_i * sdt; // [mol(NH3) m-2 s-1 * sdt!!!]


                                    // Calculate instantaneous flux first
                                    flux_inst[ij] = (-1.0) * vdnh3[ij] * nh3[ijk] * rhoref[k] * xmair_i * xmnh3; // [kg(NH3) m-2 s-1]

                                    // Then calculate accumulated flux using the instantaneous value
                                    TF flux = flux_inst[ij] * sdt; // Scale by timestep for accumulation

                                    //TF flux = (-1.0) * vdnh3[ij] * nh3[ijk] * rhoref[k] * xmair_i * xmnh3 * sdt; // [kg(NH3) m-2 s-1] 
                                    flux_nh3[ij] += flux;        // For period statistics
                                    decay = vdnh3[ij]*dzi[k] + lti;   // 1/s
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
Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(master, grid, fields)
{
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();

    sw_chemistry = inputin.get_item<bool>("chemistry", "swchemistry", "", false);
    //lifetime     = inputin.get_item<TF>("chemistry", "lifetime", "", (TF)72000);  // seconds (20 hour default)
    lifetime     = inputin.get_item<TF>("chemistry", "lifetime", "", (TF)1e30);  // seconds 
    master.print_message("Lifetime of the tracer:  = %13.5e s \n", lifetime);

    if (!sw_chemistry)
        return;

    deposition = std::make_shared<Deposition <TF>>(masterin, gridin, fieldsin, inputin);
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


        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ij = i + j*gd.jstride;
                flux_nh3[ij] /= trfa;
            } 

        stats.calc_stats_2d("flux_nh3", flux_nh3, no_offset); //added for nh3_flux
        stats.calc_stats_2d("flux_inst", flux_inst, no_offset); // added for instantaneous deposition flux of NH3


        // Reset the periodic flux after saving to stats
        trfa = 0;
        std::fill(flux_nh3.begin(), flux_nh3.end(), TF(0));

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

    // initialize 2D deposition arrays:
    vdnh3.resize(gd.ijcells);

    // added for nh3_flux (initialize nh3_flux arrays)
    flux_nh3.resize(gd.ijcells);
    std::fill(flux_nh3.begin(), flux_nh3.end(), TF(0));


    // added for instantaneous deposition flux of NH3 
    flux_inst.resize(gd.ijcells);
    std::fill(flux_inst.begin(), flux_inst.end(), TF(0));

    // initialize deposition routine:
    deposition-> init(inputin);

    // fill deposition with standard values:
    std::fill(vdnh3.begin(), vdnh3.end(), deposition-> get_vd("nh3"));

    master.print_message("Deposition arrays initialized, e.g. with vdnh3 = %13.5e m/s \n", deposition-> get_vd("nh3"));
}

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
    }

    // add cross-sections
    if (cross.get_switch())
    {
        //std::vector<std::string> allowed_crossvars = {"vdnh3"};
        std::vector<std::string> allowed_crossvars = {"vdnh3","flux_nh3","flux_inst"};
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
            cross.cross_plane(flux_nh3.data(), no_offset, name, iotime);
        else if (name == "flux_inst")  //added for instantaneous deposition flux of NH3
            cross.cross_plane(flux_inst.data(), no_offset, name, iotime);
    }

    // see if to write per tile:
    deposition->exec_cross(cross, iotime);
}

template <typename TF>
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary)
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
            vdnh3.data());
}


#ifndef USECUDA
template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo,double sdt,double dt)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    auto tmp = fields.get_tmp();
    thermo.get_thermo_field(*tmp, "T", true, false);

    // Calculate the mean temperature and water vapor mixing ratio.
    field3d_operators.calc_mean_profile(tprof.data(), tmp->fld.data());
    qprof = fields.sp.at("qt")->fld_mean;
    //field3d_operators.calc_mean_profile(qprof.data(), fields.sp.at("qt")->fld.data());


    pss<TF>(
            fields.st.at("nh3")->fld.data(), fields.sp.at("nh3")->fld.data(),
            jval, emval,
            vdnh3.data(),
            tprof.data(),
            qprof.data(),
            gd.dzi.data(),
            fields.rhoref.data(),
            rfa.data(),
            flux_nh3.data(),  //added for nh3_flux
            flux_inst.data(), //added for instantaneous deposition flux of NH3
            trfa, //accumulated for the time
            dt, sdt, lifetime,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells,
            gd.dx, gd.dy);    // Add these parameters

    fields.release_tmp(tmp);

    //isop_stat<TF>(
    //      fields.st.at("isop")->fld.data(), fields.sp.at("isop")->fld.data(),
    //      gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
    //      gd.icells, gd.ijcells);
}
#endif

template class Chemistry<double>;
//:template class Chemistry<float>;

