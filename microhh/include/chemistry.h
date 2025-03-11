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

/**
 * Class that creates a chemistry term for scalars.
 */

enum class Chemistry_type {disabled, enabled, simple};

template<typename TF>
class Chemistry
{
    public:
        Chemistry(Master&, Grid<TF>&, Fields<TF>&, Radiation<TF>&, Input&);

        ~Chemistry();                                       ///< Destructor  of the chemistry class.

        void init(Input&);                 ///< Initialize the arrays that contain the profiles.
        void create(const Timeloop<TF>&, std::string, Netcdf_handle&, Stats<TF>&, Cross<TF>&);   ///< Read the profiles of the forces from the input.
        void update_time_dependent(Timeloop<TF>&, Boundary<TF>&, Thermo<TF>&); ///< Update the time dependent parameters.
                                                                               //void update_time_dependent(Timeloop<TF>&,Boundary<TF>&); ///< Update the time dependent parameters.
        void exec(Thermo<TF>&,double,double);     ///< Add the tendencies belonging to the chemistry processes.
        void exec_stats(const int, const double, Stats<TF>&);   /// calculate statistics
        void exec_cross(Cross<TF>&, unsigned long);

    protected:
        // Cross sections
        std::vector<std::string> cross_list;         // List of active cross variables

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Radiation<TF>& radiation; // Add this line


        bool sw_chemistry;
        TF lifetime; // lifetime of species in seconds

        Field3d_operators<TF> field3d_operators;

        std::shared_ptr<Deposition<TF>> deposition;

        Mask<TF> m;     // borrow from Stats to gather statistics chemistry
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
        std::vector<TF> flux_nh3;  //added for nh3_flux 
        std::vector<TF> flux_inst; //added for instantaneous flux
        TF trfa;

        // vectors to contain calculated deposition velocities (m/s)
        std::vector<TF> vdnh3;

        const std::string tend_name = "chemistry";
        const std::string tend_longname = "Chemistry";
};
#endif


