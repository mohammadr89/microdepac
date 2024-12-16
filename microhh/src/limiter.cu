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

#include <algorithm>
#include <iostream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "limiter.h"
#include "tools.h"
#include "constants.h"

namespace
{
    template<typename TF>__global__
    void tendency_limiter(
            TF* const __restrict__ at,
            const TF* const __restrict__ a,
            const TF min_value,
            const TF dt, const TF dti,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            const TF a_new = a[ijk] + dt*at[ijk];
            at[ijk] += (a_new < min_value) ? (-a_new + min_value) * dti : TF(0.);
        }
    }
}

#ifdef USECUDA
template <typename TF>
void Limiter<TF>::exec(double dt, Stats<TF>& stats)
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dti = 1./dt;

    // Add epsilon, to make sure the final result ends just above zero.
    // NOTE: don't use `eps<TF>` here; `eps<float>` is too large
    //       as a lower limit for e.g. hydrometeors or chemical species.
    constexpr TF min_value = std::numeric_limits<double>::epsilon();

    for (auto& name : limit_list)
    {
        tendency_limiter<TF><<<gridGPU, blockGPU>>>(
            fields.at.at(name)->fld_g,
            fields.ap.at(name)->fld_g,
            min_value,
            dt, dti,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();

        stats.calc_tend(*fields.at.at(name), tend_name);
    }

    if (limit_sgstke)
    {
        tendency_limiter<TF><<<gridGPU, blockGPU>>>(
            fields.at.at("sgstke")->fld_g,
            fields.ap.at("sgstke")->fld_g,
            Constants::sgstke_min<TF>,
            dt, dti,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();

        stats.calc_tend(*fields.at.at("sgstke"), tend_name);
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Limiter<float>;
#else
template class Limiter<double>;
#endif
