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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "field3d_operators.h"
#include "timeloop.h"
#include "timedep.h"
#include "stats.h"
#include "finite_difference.h"
#include "constants.h"
#include "tools.h"
#include "boundary.h"
#include "thermo.h"
#include "force.h"

// Kernel/CUDA launcher:
#include "force_kernels.cuh"
#include "cuda_launcher.h"
#include "cuda_tiling.h"

using namespace Finite_difference::O2;

namespace
{
    template<typename TF> __global__
    void add_pressure_force_g(TF* const __restrict__ ut,
                       const TF fbody,
                       const int jj, const int kk,
                       const int istart, const int jstart, const int kstart,
                       const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] += fbody;
        }
    }

    template<typename TF> __global__
    void coriolis_4th_g(TF* const __restrict__ ut, TF* const __restrict__ vt,
                        TF* const __restrict__ u,  TF* const __restrict__ v,
                        TF* const __restrict__ ug, TF* const __restrict__ vg,
                        const TF fc, const TF ugrid, const TF vgrid,
                        const int jj, const int kk,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend)
    {
        using namespace Finite_difference::O4;

        const int i   = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j   = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k   = blockIdx.z + kstart;
        const int ii  = 1;
        const int ii2 = 2;
        const int jj2 = 2*jj;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] += fc * ( ( ci0<TF>*(ci0<TF>*v[ijk-ii2-jj ] + ci1<TF>*v[ijk-ii-jj ] + ci2<TF>*v[ijk-jj    ] + ci3<TF>*v[ijk+ii-jj  ])
                              + ci1<TF>*(ci0<TF>*v[ijk-ii2    ] + ci1<TF>*v[ijk-ii    ] + ci2<TF>*v[ijk       ] + ci3<TF>*v[ijk+ii     ])
                              + ci2<TF>*(ci0<TF>*v[ijk-ii2+jj ] + ci1<TF>*v[ijk-ii+jj ] + ci2<TF>*v[ijk+jj    ] + ci3<TF>*v[ijk+ii+jj  ])
                              + ci3<TF>*(ci0<TF>*v[ijk-ii2+jj2] + ci1<TF>*v[ijk-ii+jj2] + ci2<TF>*v[ijk+jj2   ] + ci3<TF>*v[ijk+ii+jj2 ]) )
                       + vgrid - vg[k] );

            vt[ijk] -= fc * ( ( ci0<TF>*(ci0<TF>*u[ijk-ii-jj2 ] + ci1<TF>*u[ijk-jj2   ] + ci2<TF>*u[ijk+ii-jj2] + ci3<TF>*u[ijk+ii2-jj2])
                              + ci1<TF>*(ci0<TF>*u[ijk-ii-jj  ] + ci1<TF>*u[ijk-jj    ] + ci2<TF>*u[ijk+ii-jj ] + ci3<TF>*u[ijk+ii2-jj ])
                              + ci2<TF>*(ci0<TF>*u[ijk-ii     ] + ci1<TF>*u[ijk       ] + ci2<TF>*u[ijk+ii    ] + ci3<TF>*u[ijk+ii2    ])
                              + ci3<TF>*(ci0<TF>*u[ijk-ii+jj  ] + ci1<TF>*u[ijk+jj    ] + ci2<TF>*u[ijk+ii+jj ] + ci3<TF>*u[ijk+ii2+jj ]) )
                       + ugrid - ug[k]);
        }
    }


    template<typename TF> __global__
    void nudging_tendency_g(
        TF* const __restrict__ st,
        const TF* const __restrict__ smn,
        const TF* const __restrict__ snudge,
        const TF* const __restrict__ nudge_fac,
        const int kstart, const int kend)
    {
        const int k = blockIdx.x*blockDim.x + threadIdx.x + kstart;

        if (k < kend)
            st[k] = -nudge_fac[k] * (smn[k]-snudge[k]);
    }


    template<typename TF>
    int calc_zi(
            const TF* const restrict fldmean,
            const int kstart, const int kend,
            const int plusminus)
    {
        TF maxgrad = 0.;
        TF grad = 0.;
        int kinv = kstart;
        for (int k=kstart+1; k<kend; ++k)
        {
            grad = plusminus * (fldmean[k] - fldmean[k-1]);
            if (grad > maxgrad)
            {
                maxgrad = grad;
                kinv = k;
            }
        }
        return kinv;
    }

    template<typename TF>
    void rescale_nudgeprof(TF* const restrict fldmean, const int kinv, const int kstart, const int kend)
    {
        for (int k=kstart+1; k<kinv; ++k)
            fldmean[k] = fldmean[kstart];

        for (int k=kinv+1; k<kend-2; ++k)
            fldmean[k] = fldmean[kend-1];
    }

    template<typename TF> __global__
    void calc_time_dependent_prof_g(
            TF* const __restrict__ prof, const TF* const __restrict__ data,
            const double fac0, const double fac1,
            const int index0,  const int index1,
            const int kmax,    const int kgc)
    {
        const int k = blockIdx.x*blockDim.x + threadIdx.x;
        const int kk = kmax;

        if (k < kmax)
            prof[k+kgc] = fac0*data[index0*kk+k] + fac1*data[index1*kk+k];
    }
} // end namespace

template<typename TF>
void Force<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const int nmemsize  = gd.kcells*sizeof(TF);

    if (swlspres == Large_scale_pressure_type::Geo_wind)
    {
        ug_g.allocate(gd.kcells);
        vg_g.allocate(gd.kcells);

        cuda_safe_call(cudaMemcpy(ug_g, ug.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(vg_g, vg.data(), nmemsize, cudaMemcpyHostToDevice));
    }

    if (swls == Large_scale_tendency_type::Enabled)
    {
        for (auto& it : lslist)
        {
            lsprofs_g.emplace(it, cuda_vector<TF>(gd.kcells));
            cuda_safe_call(cudaMemcpy(lsprofs_g.at(it), lsprofs.at(it).data(), nmemsize, cudaMemcpyHostToDevice));
        }
    }

    if (swnudge == Nudging_type::Enabled)
    {
        for (auto& it : nudgelist)
        {
            nudgeprofs_g.emplace(it, cuda_vector<TF>(gd.kcells));
            cuda_safe_call(cudaMemcpy(nudgeprofs_g.at(it), nudgeprofs.at(it).data(), nmemsize, cudaMemcpyHostToDevice));
        }
        nudge_factor_g.allocate(gd.kcells);
        nudge_tend_g.allocate(gd.kcells);
        cuda_safe_call(cudaMemcpy(nudge_factor_g, nudge_factor.data(), nmemsize, cudaMemcpyHostToDevice));
    }

    if (swwls == Large_scale_subsidence_type::Mean_field ||
        swwls == Large_scale_subsidence_type::Local_field)
    {
        wls_g.allocate(gd.kcells);
        cuda_safe_call(cudaMemcpy(wls_g, wls.data(), nmemsize, cudaMemcpyHostToDevice));
    }
}

template<typename TF>
void Force<TF>::clear_device()
{
    if (swlspres == Large_scale_pressure_type::Geo_wind)
    {
        for (auto& it : tdep_geo)
            it.second->clear_device();
    }

    if (swls == Large_scale_tendency_type::Enabled)
    {
        for (auto& it : tdep_ls)
            it.second->clear_device();
    }

    if (swnudge == Nudging_type::Enabled)
    {
        for (auto& it : tdep_nudge)
            it.second->clear_device();
    }

    if (swwls == Large_scale_subsidence_type::Mean_field ||
        swwls == Large_scale_subsidence_type::Local_field)
    {
        tdep_wls->clear_device();
    }
}

#ifdef USECUDA
template<typename TF>
void Force<TF>::exec(double dt, Thermo<TF>& thermo, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Grid layout for KL/CL launches over interior.
    Grid_layout grid_layout = {
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.istride,
            gd.jstride,
            gd.kstride};

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    if (swlspres == Large_scale_pressure_type::Fixed_flux)
    {
        auto tmp = fields.get_tmp_g();

        TF uavg  = field3d_operators.calc_mean_g(fields.mp.at("u")->fld_g);
        TF utavg = field3d_operators.calc_mean_g(fields.mt.at("u")->fld_g);

        fields.release_tmp_g(tmp);

        const TF fbody = (uflux - uavg - gd.utrans) / dt - utavg;

        add_pressure_force_g<TF><<<gridGPU, blockGPU>>>(
            fields.mt.at("u")->fld_g,
            fbody,
            gd.icells, gd.ijcells,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend);
        cuda_check_error();
        cudaDeviceSynchronize();

        stats.calc_tend(*fields.mt.at("u"), tend_name_pres);
    }
    else if (swlspres == Large_scale_pressure_type::Pressure_gradient)
    {
        const TF fbody = TF(-1.)*dpdx;
        add_pressure_force_g<TF><<<gridGPU, blockGPU>>>(
            fields.mt.at("u")->fld_g, fbody,
            gd.icells, gd.ijcells,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend);
        cuda_check_error();
        cudaDeviceSynchronize();

        stats.calc_tend(*fields.mt.at("u"), tend_name_pres);
    }
    else if (swlspres == Large_scale_pressure_type::Geo_wind)
    {
        TF fc_loc = fc;
        if (fc_loc < 0)
            fc_loc = 2. * Constants::e_rot<TF> * std::sin(gd.lat * TF(M_PI) / 180.);
        
        if (grid.get_spatial_order() == Grid_order::Second)
        {
            launch_grid_kernel<Force_kernels::coriolis_2nd_g<TF>>(
                    grid_layout,
                    fields.mt.at("u")->fld_g.view(),
                    fields.mt.at("v")->fld_g.view(),
                    fields.mp.at("u")->fld_g,
                    fields.mp.at("v")->fld_g,
                    ug_g, vg_g, fc,
                    gd.utrans, gd.vtrans);
        }
        else if (grid.get_spatial_order() == Grid_order::Fourth)
        {
            coriolis_4th_g<TF><<<gridGPU, blockGPU>>>(
                fields.mt.at("u")->fld_g, fields.mt.at("v")->fld_g,
                fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
                ug_g, vg_g, fc, gd.utrans, gd.vtrans,
                gd.icells, gd.ijcells,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend);
            cuda_check_error();
        }
        cudaDeviceSynchronize();

        stats.calc_tend(*fields.mt.at("u"), tend_name_cor);
        stats.calc_tend(*fields.mt.at("v"), tend_name_cor);
    }

    if (swls == Large_scale_tendency_type::Enabled)
    {
        for (auto& it : lslist)
        {
            launch_grid_kernel<Force_kernels::add_profile_g<TF>>(
                    grid_layout,
                    fields.at.at(it)->fld_g.view(),
                    lsprofs_g.at(it));

            cudaDeviceSynchronize();
            stats.calc_tend(*fields.at.at(it), tend_name_ls);
        }
    }

    if (swnudge == Nudging_type::Enabled)
    {
        for (auto& it : nudgelist)
        {
            auto it1 = std::find(scalednudgelist.begin(), scalednudgelist.end(), it);
            if (it1 != scalednudgelist.end())
            {
                cudaMemcpy(fields.ap.at(it)->fld_mean.data(), fields.ap.at(it)->fld_mean_g, gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost);
                const int kinv = thermo.get_bl_depth();
                rescale_nudgeprof(nudgeprofs.at(it).data(), kinv, gd.kstart, gd.kend);
                cudaMemcpy(nudgeprofs_g.at(it), nudgeprofs.at(it).data(), gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
            }

            // Calculate nudging tendency profile.
            const int blocki = 32;
            const int gridi  = gd.kmax/blocki + (gd.kmax%blocki > 0);

            dim3 gridGPU (gridi);
            dim3 blockGPU(blocki);

            // Calculate nudging tendency as single profile.
            nudging_tendency_g<TF><<<gridGPU, blockGPU>>>(
                    nudge_tend_g.view(),
                    fields.ap.at(it)->fld_mean_g,
                    nudgeprofs_g.at(it),
                    nudge_factor_g,
                    gd.kstart, gd.kend);

            // Add tendency profile to 3D tendency field.
            launch_grid_kernel<Force_kernels::add_profile_g<TF>>(
                    grid_layout,
                    fields.at.at(it)->fld_g.view(),
                    nudge_tend_g);

            cudaDeviceSynchronize();
            stats.calc_tend(*fields.at.at(it), tend_name_nudge);
        }
    }

    if (swwls == Large_scale_subsidence_type::Mean_field)
    {
        if (swwls_mom)
        {
            launch_grid_kernel<Force_kernels::advec_wls_2nd_mean_g<TF>>(
                    grid_layout,
                    fields.mt.at("u")->fld_g.view(),
                    fields.mp.at("u")->fld_mean_g,
                    wls_g, gd.dzhi_g);

            cudaDeviceSynchronize();
            stats.calc_tend(*fields.mt.at("u"), tend_name_subs);

            launch_grid_kernel<Force_kernels::advec_wls_2nd_mean_g<TF>>(
                    grid_layout,
                    fields.mt.at("v")->fld_g.view(),
                    fields.mp.at("v")->fld_mean_g,
                    wls_g, gd.dzhi_g);

            cudaDeviceSynchronize();
            stats.calc_tend(*fields.mt.at("v"), tend_name_subs);
        }

        for (auto& it : fields.st)
        {
            launch_grid_kernel<Force_kernels::advec_wls_2nd_mean_g<TF>>(
                    grid_layout,
                    fields.st.at(it.first)->fld_g.view(),
                    fields.sp.at(it.first)->fld_mean_g,
                    wls_g, gd.dzhi_g);

            cudaDeviceSynchronize();
            stats.calc_tend(*it.second, tend_name_subs);
        }
    }
    else if (swwls == Large_scale_subsidence_type::Local_field)
    {
        if (swwls_mom)
        {
            launch_grid_kernel<Force_kernels::advec_wls_2nd_local_g<TF>>(
                    grid_layout,
                    fields.mt.at("u")->fld_g.view(),
                    fields.mp.at("u")->fld_g,
                    wls_g, gd.dzhi_g);

            cudaDeviceSynchronize();
            stats.calc_tend(*fields.mt.at("u"), tend_name_subs);

            launch_grid_kernel<Force_kernels::advec_wls_2nd_local_g<TF>>(
                    grid_layout,
                    fields.mt.at("v")->fld_g.view(),
                    fields.mp.at("v")->fld_g,
                    wls_g, gd.dzhi_g);

            cudaDeviceSynchronize();
            stats.calc_tend(*fields.mt.at("v"), tend_name_subs);

            // Modified `grid_layout`, which starts at kstart+1 (for w subsidence).
            Grid_layout grid_layout_kp1 = {
                    gd.istart,   gd.iend,
                    gd.jstart,   gd.jend,
                    gd.kstart+1, gd.kend,
                    gd.istride,
                    gd.jstride,
                    gd.kstride};

            launch_grid_kernel<Force_kernels::advec_wls_2nd_local_w_g<TF>>(
                    grid_layout_kp1,
                    fields.mt.at("w")->fld_g.view(),
                    fields.mp.at("w")->fld_g,
                    wls_g, gd.dzi_g);

            cudaDeviceSynchronize();
            stats.calc_tend(*fields.mt.at("w"), tend_name_subs);
        }

        for (auto& it : fields.st)
        {
            launch_grid_kernel<Force_kernels::advec_wls_2nd_local_g<TF>>(
                    grid_layout,
                    fields.st.at(it.first)->fld_g.view(),
                    fields.sp.at(it.first)->fld_g,
                    wls_g, gd.dzhi_g);

            cudaDeviceSynchronize();
            stats.calc_tend(*it.second, tend_name_subs);
        }
    }
}
#endif

#ifdef USECUDA
template <typename TF>
void Force<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (swls == Large_scale_tendency_type::Enabled)
    {
        for (auto& it : tdep_ls)
            it.second->update_time_dependent_prof_g(lsprofs_g.at(it.first), timeloop);
    }

    if (swnudge == Nudging_type::Enabled)
    {
        for (auto& it : tdep_nudge)
            it.second->update_time_dependent_prof_g(nudgeprofs_g.at(it.first), timeloop);
    }

    if (swlspres == Large_scale_pressure_type::Geo_wind)
    {
        tdep_geo.at("u_geo")->update_time_dependent_prof_g(ug_g, timeloop);
        tdep_geo.at("v_geo")->update_time_dependent_prof_g(vg_g, timeloop);
    }

    if (swwls == Large_scale_subsidence_type::Mean_field ||
        swwls == Large_scale_subsidence_type::Local_field)
    {
        tdep_wls->update_time_dependent_prof_g(wls_g, timeloop);
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Force<float>;
#else
template class Force<double>;
#endif
