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
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "advec_2i62.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "advec_monotonic.h"

template<typename TF>
Advec_2i62<TF>::Advec_2i62(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Advec<TF>(masterin, gridin, fieldsin, inputin)
{
    fluxlimit_list = inputin.get_list<std::string>(
            "advec", "fluxlimit_list", "", std::vector<std::string>());

    const int igc = 3;
    const int jgc = 3;
    const int kgc = (fluxlimit_list.empty()) ? 1 : 2;
    grid.set_minimum_ghost_cells(igc, jgc, kgc);
}

template<typename TF>
Advec_2i62<TF>::~Advec_2i62() {}


namespace
{
    using namespace Finite_difference::O2;
    using namespace Finite_difference::O6;
    using namespace Advec_monotonic;

    template<typename TF>
    TF calc_cfl(
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzi,
            const TF dx, const TF dy,
            const TF dt, Master& master,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        const int jj1 = jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;

        const int kk1 = kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxi = 1./dx;
        const TF dyi = 1./dy;

        TF cfl = 0;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    cfl = std::max(cfl, std::abs(interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]))*dxi
                                      + std::abs(interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]))*dyi
                                      + std::abs(interp2(w[ijk], w[ijk+kk1]))*dzi[k]);
                }

        master.max(&cfl, 1);
        cfl = cfl*dt;
        return cfl;
    }

    template<typename TF>
    void advec_u(
            TF* const restrict ut,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzi,
            const TF dx, const TF dy,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        const int jj1 = jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;

        const int kk1 = kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        // Calculate horizontal terms
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    ut[ijk] +=
                            // u*du/dx
                            - ( interp2(u[ijk        ], u[ijk+ii1]) * interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3])
                              - interp2(u[ijk-ii1    ], u[ijk    ]) * interp6_ws(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) ) * dxi

                            // v*du/dy
                            - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp6_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3])
                              - interp2(v[ijk-ii1    ], v[ijk    ]) * interp6_ws(u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2]) ) * dyi

                            - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1])
                              - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k];
                }
    }

    template<typename TF>
    void advec_v(
            TF* const restrict vt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzi,
            const TF dx, const TF dy,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        const int jj1 = jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;

        const int kk1 = kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        // Calculate horizontal terms
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    vt[ijk] +=
                            // u*dv/dx
                            - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp6_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2], v[ijk+ii3])
                              - interp2(u[ijk    -jj1], u[ijk    ]) * interp6_ws(v[ijk-ii3], v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2]) ) * dxi

                            // v*dv/dy
                            - ( interp2(v[ijk        ], v[ijk+jj1]) * interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3])
                              - interp2(v[ijk-jj1    ], v[ijk    ]) * interp6_ws(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) ) * dyi

                            - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1])
                              - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k];
                }
    }

    template<typename TF>
    void advec_w(
            TF* const restrict wt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzhi,
            const TF dx, const TF dy,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        const int jj1 = jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;

        const int kk1 = kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        // Calculate horizontal terms
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    wt[ijk] +=
                            // u*dw/dx
                            - ( interp2(u[ijk+ii1-kk], u[ijk+ii1]) * interp6_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3])
                              - interp2(u[ijk    -kk], u[ijk    ]) * interp6_ws(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2]) ) * dxi

                            // v*dw/dy
                            - ( interp2(v[ijk+jj1-kk], v[ijk+jj1]) * interp6_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3])
                              - interp2(v[ijk    -kk], v[ijk    ]) * interp6_ws(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2]) ) * dyi

                            - ( rhoref[k  ] * interp2(w[ijk    ], w[ijk+kk1]) * interp2(w[ijk    ], w[ijk+kk1])
                              - rhoref[k-1] * interp2(w[ijk-kk1], w[ijk    ]) * interp2(w[ijk-kk1], w[ijk    ]) ) / rhorefh[k] * dzhi[k];
                }
    }

    template<typename TF>
    void advec_s(
            TF* const restrict st,
            const TF* const restrict s,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzi,
            const TF dx, const TF dy,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        const int jj1 = jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;

        const int kk1 = kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        // Calculate horizontal terms
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    st[ijk] +=
                            - ( u[ijk+ii1] * interp6_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2], s[ijk+ii3])
                              - u[ijk    ] * interp6_ws(s[ijk-ii3], s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2]) ) * dxi

                            - ( v[ijk+jj1] * interp6_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2], s[ijk+jj3])
                              - v[ijk    ] * interp6_ws(s[ijk-jj3], s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2]) ) * dyi

                            - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1])
                              - rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];
                }
    }

    template<typename TF>
    void advec_flux_u(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = interp2(w[ijk-ii], w[ijk]) * interp2(s[ijk-kk], s[ijk]);
                }
    }


    template<typename TF>
    void advec_flux_v(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = interp2(w[ijk-jj], w[ijk]) * interp2(s[ijk-kk], s[ijk]);
                }
    }

    template<typename TF>
    void advec_flux_w(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = w[ijk] * s[ijk];
                }
    }

    template<typename TF>
    void advec_flux_s(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = w[ijk] * interp2(s[ijk-kk], s[ijk]);
                }
    }
}

template<typename TF>
void Advec_2i62<TF>::create(Stats<TF>& stats)
{
    for (auto& s : fields.sp)
    {
        if (std::find(fluxlimit_list.begin(), fluxlimit_list.end(), s.first) != fluxlimit_list.end())
            sp_limit.push_back(s.first);
        else
            sp_no_limit.push_back(s.first);
    }

    stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);
    for (auto it : fields.st)
        stats.add_tendency(*it.second, "z", tend_name, tend_longname);
}

#ifndef USECUDA
template<typename TF>
double Advec_2i62<TF>::get_cfl(double dt)
{
    auto& gd = grid.get_grid_data();
    TF cfl = calc_cfl<TF>(
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            dt, master,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    // CFL is kept in double precision for time stepping accuracy
    return static_cast<double>(cfl);
}

template<typename TF>
unsigned long Advec_2i62<TF>::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    auto& gd = grid.get_grid_data();

    double cfl = calc_cfl<TF>(
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            dt, master,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    cfl = std::max(cflmin, cfl);
    return idt * cflmax / cfl;
}


template<typename TF>
void Advec_2i62<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    advec_u(fields.mt.at("u")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    advec_v(fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    advec_w(fields.mt.at("w")->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            gd.dzhi.data(), gd.dx, gd.dy,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    for (const std::string& s : sp_no_limit)
    {
        advec_s(fields.st.at(s)->fld.data(),
                fields.sp.at(s)->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(),
                gd.dzi.data(), gd.dx, gd.dy,
                fields.rhoref.data(), fields.rhorefh.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }

    for (const std::string& s : sp_limit)
    {
        advec_s_lim(
                fields.st.at(s)->fld.data(), fields.sp.at(s)->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzi.data(), gd.dx, gd.dy,
                fields.rhoref.data(), fields.rhorefh.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}
#endif

template<typename TF>
void Advec_2i62<TF>::get_advec_flux(
        Field3d<TF>& advec_flux, const Field3d<TF>& fld)
{
    auto& gd = grid.get_grid_data();

    if (fld.loc == gd.uloc)
    {
        advec_flux_u(
                advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else if (fld.loc == gd.vloc)
    {
        advec_flux_v(
                advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else if (fld.loc == gd.wloc)
    {
        advec_flux_w(
                advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else if (fld.loc == gd.sloc)
    {
        if (std::find(sp_limit.begin(), sp_limit.end(), fld.name) != sp_limit.end())
            advec_flux_s_lim(
                    advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        else
            advec_flux_s(
                    advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }
    else
        throw std::runtime_error("Advec_2i62 cannot deliver flux field at that location");
}

#ifdef FLOAT_SINGLE
template class Advec_2i62<float>;
#else
template class Advec_2i62<double>;
#endif
