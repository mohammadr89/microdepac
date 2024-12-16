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

#ifndef DEPOSITION_H
#define DEPOSITION_H

#include <vector>
#include <string>
#include <map>
#include "timeloop.h"
#include "boundary_surface_lsm.h"
#include "boundary.h"


class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;
template<typename> class Boundary_surface_lsm;

/**
 * Class that creates a deposition liniked to the chemistry
 */

enum class Deposition_type {disabled, enabled, simple, average};

template<typename TF>
struct Deposition_tile
{
    std::string long_name;    // Descriptive name of tile
    // Land surface
    std::vector<TF> vdnh3;     // deposition velocity of ozone (m s-1)
};

template<typename TF>
using Deposition_tile_map = std::map<std::string, Deposition_tile<TF>>;

template<typename TF>
class Deposition
{
    public:
        Deposition(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the chemistry class.
        ~Deposition();                                       ///< Destructor  of the chemistry class.

        void init(Input&);                 ///< Initialize the arrays that contain the profiles.
        void create(Stats<TF>&, Cross<TF>&);
        void update_time_dependent(Timeloop<TF>&, Boundary<TF>&,
             TF*); ///< Update the time dependent deposition parameters.

        const TF get_vd(const std::string&) const;                  ///< get the standard vd value (nh3, no, no2, ..)
        void get_tiled_mean(TF*, std::string, TF, const TF*, const TF*, const TF*);
        void update_vd_water(TF*, std::string, const TF*, const TF*, const int*, const TF*, const TF*);
        void exec_cross(Cross<TF>&, unsigned long);
        void spatial_avg_vd(TF*);                           //MAQ_AV_21042022+ added spatial_avg_vd

    protected:
        std::vector<std::string> cross_list;         // List of active cross variables

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_deposition;

        std::shared_ptr<Boundary_surface_lsm<TF>> boundary_surface_lsm;

        TF deposition_var;   // here put the local vars
        TF henry_so2;
        TF rsoil_so2;
        TF rwat_so2;
        TF rws_so2;
        //TF lai;
        std::vector<TF> rmes;
        std::vector<TF> rsoil;
        std::vector<TF> rcut;
        std::vector<TF> rws;
        std::vector<TF> rwat;
        std::vector<TF> diff;
        std::vector<TF> diff_scl;
        std::vector<TF> henry;
        std::vector<TF> f0;

        TF vd_nh3;

        std::vector<std::string> deposition_tile_names {"veg", "soil" ,"wet"};
        Deposition_tile_map<TF> deposition_tiles;

        const std::string tend_name = "deposition";
        const std::string tend_longname = "Deposition";
};
#endif
