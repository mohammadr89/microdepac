/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#include <iostream>
#include <algorithm>

#include "aerosol.h"
#include "timeloop.h"
#include "input.h"
#include "grid.h"
#include "netcdf_interface.h"
#include "stats.h"
#include "thermo.h"
#include "fields.h"
#include "timedep.h"
#include "Gas_concs.h"
#include "Array.h"

using Aerosol_concs = Gas_concs;

template<typename TF>
Aerosol<TF>::Aerosol(
    Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    // Read `.ini` settings.
    sw_aerosol = inputin.get_item<bool>("aerosol", "swaerosol", "", 0);
    sw_timedep = inputin.get_item<bool>("aerosol", "swtimedep", "", 0);
}

template <typename TF>
Aerosol<TF>::~Aerosol()
{
}

template <typename TF>
void Aerosol<TF>::init()
{
    // Allocate (`.resize`) arrays.
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

    aermr01.resize(gd.kcells);
    aermr02.resize(gd.kcells);
    aermr03.resize(gd.kcells);
    aermr04.resize(gd.kcells);
    aermr05.resize(gd.kcells);
    aermr06.resize(gd.kcells);
    aermr07.resize(gd.kcells);
    aermr08.resize(gd.kcells);
    aermr09.resize(gd.kcells);
    aermr10.resize(gd.kcells);
    aermr11.resize(gd.kcells);

}

template <typename TF>
void Aerosol<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    // Read input from NetCDF and prepare statistics output.
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

    if (sw_timedep)
    {
        // create time dependent profiles
        const TF offset = 0;
        std::string timedep_dim = "time_rad";

        tdep_aermr01 = std::make_unique<Timedep<TF>>(master, grid, "aermr01", sw_timedep);
        tdep_aermr01->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr02 = std::make_unique<Timedep<TF>>(master, grid, "aermr02", sw_timedep);
        tdep_aermr02->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr03 = std::make_unique<Timedep<TF>>(master, grid, "aermr03", sw_timedep);
        tdep_aermr03->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr04 = std::make_unique<Timedep<TF>>(master, grid, "aermr04", sw_timedep);
        tdep_aermr04->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr05 = std::make_unique<Timedep<TF>>(master, grid, "aermr05", sw_timedep);
        tdep_aermr05->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr06 = std::make_unique<Timedep<TF>>(master, grid, "aermr06", sw_timedep);
        tdep_aermr06->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr07 = std::make_unique<Timedep<TF>>(master, grid, "aermr07", sw_timedep);
        tdep_aermr07->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr08 = std::make_unique<Timedep<TF>>(master, grid, "aermr08", sw_timedep);
        tdep_aermr08->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr09 = std::make_unique<Timedep<TF>>(master, grid, "aermr09", sw_timedep);
        tdep_aermr09->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr10 = std::make_unique<Timedep<TF>>(master, grid, "aermr10", sw_timedep);
        tdep_aermr10->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr11 = std::make_unique<Timedep<TF>>(master, grid, "aermr11", sw_timedep);
        tdep_aermr11->create_timedep_prof(input_nc, offset, timedep_dim);
    }
    else
    {
        // Read NetCDF input.
        Netcdf_group& group_nc = input_nc.get_group("init");
        group_nc.get_variable(aermr01, "aermr01", {0}, {gd.ktot});
        std::rotate(aermr01.rbegin(), aermr01.rbegin()+gd.kstart, aermr01.rend());
        group_nc.get_variable(aermr02, "aermr02", {0}, {gd.ktot});
        std::rotate(aermr02.rbegin(), aermr02.rbegin()+gd.kstart, aermr02.rend());
        group_nc.get_variable(aermr03, "aermr03", {0}, {gd.ktot});
        std::rotate(aermr03.rbegin(), aermr03.rbegin()+gd.kstart, aermr03.rend());
        group_nc.get_variable(aermr04, "aermr04", {0}, {gd.ktot});
        std::rotate(aermr04.rbegin(), aermr04.rbegin()+gd.kstart, aermr04.rend());
        group_nc.get_variable(aermr05, "aermr05", {0}, {gd.ktot});
        std::rotate(aermr05.rbegin(), aermr05.rbegin()+gd.kstart, aermr05.rend());
        group_nc.get_variable(aermr06, "aermr06", {0}, {gd.ktot});
        std::rotate(aermr06.rbegin(), aermr06.rbegin()+gd.kstart, aermr06.rend());
        group_nc.get_variable(aermr07, "aermr07", {0}, {gd.ktot});
        std::rotate(aermr07.rbegin(), aermr07.rbegin()+gd.kstart, aermr07.rend());
        group_nc.get_variable(aermr08, "aermr08", {0}, {gd.ktot});
        std::rotate(aermr08.rbegin(), aermr08.rbegin()+gd.kstart, aermr08.rend());
        group_nc.get_variable(aermr09, "aermr09", {0}, {gd.ktot});
        std::rotate(aermr09.rbegin(), aermr09.rbegin()+gd.kstart, aermr09.rend());
        group_nc.get_variable(aermr10, "aermr10", {0}, {gd.ktot});
        std::rotate(aermr10.rbegin(), aermr10.rbegin()+gd.kstart, aermr10.rend());
        group_nc.get_variable(aermr11, "aermr11", {0}, {gd.ktot});
        std::rotate(aermr11.rbegin(), aermr11.rbegin()+gd.kstart, aermr11.rend());
    }
}

#ifndef USECUDA
template <typename TF>
void Aerosol<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_aerosol)
        return;

    if (sw_timedep)
    {
        tdep_aermr01 ->update_time_dependent_prof(aermr01, timeloop);
        tdep_aermr02 ->update_time_dependent_prof(aermr02, timeloop);
        tdep_aermr03 ->update_time_dependent_prof(aermr03, timeloop);
        tdep_aermr04 ->update_time_dependent_prof(aermr04, timeloop);
        tdep_aermr05 ->update_time_dependent_prof(aermr05, timeloop);
        tdep_aermr06 ->update_time_dependent_prof(aermr06, timeloop);
        tdep_aermr07 ->update_time_dependent_prof(aermr07, timeloop);
        tdep_aermr08 ->update_time_dependent_prof(aermr08, timeloop);
        tdep_aermr09 ->update_time_dependent_prof(aermr09, timeloop);
        tdep_aermr10 ->update_time_dependent_prof(aermr10, timeloop);
        tdep_aermr11 ->update_time_dependent_prof(aermr11, timeloop);
    }
}
#endif

#ifndef USECUDA
template<typename TF>
void Aerosol<TF>::get_radiation_fields(Aerosol_concs& aerosol_concs)
{
    auto& gd = grid.get_grid_data();

    Array<Float,2> aermr01_a({1, gd.ktot});
    Array<Float,2> aermr02_a({1, gd.ktot});
    Array<Float,2> aermr03_a({1, gd.ktot});
    Array<Float,2> aermr04_a({1, gd.ktot});
    Array<Float,2> aermr05_a({1, gd.ktot});
    Array<Float,2> aermr06_a({1, gd.ktot});
    Array<Float,2> aermr07_a({1, gd.ktot});
    Array<Float,2> aermr08_a({1, gd.ktot});
    Array<Float,2> aermr09_a({1, gd.ktot});
    Array<Float,2> aermr10_a({1, gd.ktot});
    Array<Float,2> aermr11_a({1, gd.ktot});

    for (int k=gd.kstart; k<gd.kend; ++k)
    {
        const int k_nogc = k-gd.kgc+1;

        aermr01_a({1, k_nogc}) = aermr01[k];
        aermr02_a({1, k_nogc}) = aermr02[k];
        aermr03_a({1, k_nogc}) = aermr03[k];
        aermr04_a({1, k_nogc}) = aermr04[k];
        aermr05_a({1, k_nogc}) = aermr05[k];
        aermr06_a({1, k_nogc}) = aermr06[k];
        aermr07_a({1, k_nogc}) = aermr07[k];
        aermr08_a({1, k_nogc}) = aermr08[k];
        aermr09_a({1, k_nogc}) = aermr09[k];
        aermr10_a({1, k_nogc}) = aermr10[k];
        aermr11_a({1, k_nogc}) = aermr11[k];
    }

    aerosol_concs.set_vmr("aermr01", aermr01_a);
    aerosol_concs.set_vmr("aermr02", aermr02_a);
    aerosol_concs.set_vmr("aermr03", aermr03_a);
    aerosol_concs.set_vmr("aermr04", aermr04_a);
    aerosol_concs.set_vmr("aermr05", aermr05_a);
    aerosol_concs.set_vmr("aermr06", aermr06_a);
    aerosol_concs.set_vmr("aermr07", aermr07_a);
    aerosol_concs.set_vmr("aermr08", aermr08_a);
    aerosol_concs.set_vmr("aermr09", aermr09_a);
    aerosol_concs.set_vmr("aermr10", aermr10_a);
    aerosol_concs.set_vmr("aermr11", aermr11_a);
}
#endif


#ifdef FLOAT_SINGLE
template class Aerosol<float>;
#else
template class Aerosol<double>;
#endif
