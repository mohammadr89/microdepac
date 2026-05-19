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

#ifndef CHEMISTRY_H
#define CHEMISTRY_H
#include <vector>
#include <string>
#include <map>
#include "timeloop.h"
#include "boundary.h"
#include "stats.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Soil_grid;
template<typename> class Fields;
template<typename> class Stats;
template<typename> class Deposition;

enum class Chemistry_type {disabled, enabled, simple};

template<typename TF>
class Chemistry
{
    public:

    Chemistry(Master&, Grid<TF>&, Fields<TF>&, Radiation<TF>&, Input&);
    ~Chemistry();

    void init(Input&);
    void create(const Timeloop<TF>&, std::string, Netcdf_handle&, Stats<TF>&, Cross<TF>&);
    void update_time_dependent(Timeloop<TF>&, Boundary<TF>&, Thermo<TF>&);
    void exec(Thermo<TF>&, Boundary<TF>&, double, double);
    void exec_stats(const int, const double, Stats<TF>&);
    void exec_cross(Cross<TF>&, unsigned long);
    void calc_c_target(Boundary<TF>&);
    const std::vector<TF>& get_c_target() const { return c_target; }
    const std::vector<TF>& get_rho_target() const { return rho_target; }

    protected:

    std::vector<std::string> cross_list;

    private:

    Master& master;
    Grid<TF>& grid;
    Fields<TF>& fields;
    Radiation<TF>& radiation;
    bool sw_chemistry;
    TF lifetime;
    std::vector<TF> c_target;
    TF rsl_ratio; // Roughness sublayer ratio (alpha = z*/h_c) - Basu & Lacser (2017); reference height z_ref >= rsl_ratio * z0m
    bool sw_const_ref_height;
    TF z_fixed;
    bool sw_adapt_ref_height;
    bool sw_use_lowest_levels;
    bool sw_force_ref_height;
    Field3d_operators<TF> field3d_operators;
    std::shared_ptr<Deposition<TF>> deposition;
    Mask<TF> m;     // borrow from Stats to gather chemistry statistics
    int statistics_counter;
    std::vector<std::string> jname={"jo31d","jh2o2","jno2","jno3","jn2o5","jch2or","jch2om","jch3o2h"};
    std::vector<std::string> ename={"emi_isop","emi_no"};
    TF jval[8];   // time-interpolated value to pass to the chemistry routine
    TF emval[2];
    std::vector<TF> time;
    std::vector<TF> jo31d;
    std::vector<TF> jh2o2;
    std::vector<TF> jno2;
    std::vector<TF> jno3;
    std::vector<TF> jn2o5;
    std::vector<TF> jch2or;
    std::vector<TF> jch2om;
    std::vector<TF> jch3o2h;
    std::vector<TF> emi_isop;
    std::vector<TF> emi_no;
    std::vector<TF> rfa;
    std::vector<TF> rka;
    std::vector<TF> qprof;
    std::vector<TF> tprof;
    std::vector<TF> flux_nh3;
    std::vector<TF> flux_inst;
    std::vector<TF> total_flux_nh3;
    TF trfa;
    std::vector<TF> vdnh3;
    const std::string tend_name = "chemistry";
    const std::string tend_longname = "Chemistry";
    std::vector<TF> cstar1;
    std::vector<TF> cstar2;
    std::vector<TF> c_diff;
    std::vector<TF> c_extrap_diff;
    std::vector<TF> T_target;
    std::vector<TF> rho_target;
};
#endif
