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

#ifndef RADIATION_RRTMGP_H
#define RADIATION_RRTMGP_H

#include "radiation.h"
#include "field3d_operators.h"
#include "boundary_cyclic.h"

#include "Gas_concs.h"
#include "Gas_optics_rrtmgp.h"
#include "Source_functions.h"
#include "Cloud_optics.h"
#include "Aerosol_optics.h"
#include "Rte_lw.h"
#include "Rte_sw.h"
#include "types.h"


class Master;
class Input;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Thermo;
template<typename> class Timeloop;
template<typename TF> class Aerosol;
template<typename TF> class Background;

using Aerosol_concs = Gas_concs;
#ifdef USECUDA
using Aerosol_concs_gpu = Gas_concs_gpu;
#endif

template<typename TF>
class Radiation_rrtmgp : public Radiation<TF>
{
    public:
        Radiation_rrtmgp(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Radiation_rrtmgp() {}

        bool check_field_exists(const std::string& name)
        { throw std::runtime_error("\"check_field_exists()\" is not implemented in radiation_rrtmpg"); }

        void init(Timeloop<TF>&);
        void create(
                Input&, Netcdf_handle&, Thermo<TF>&,
                Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(Thermo<TF>&, double, Timeloop<TF>&, Stats<TF>&, Aerosol<TF>&, Background<TF>&, Microphys<TF>&);

        unsigned long get_time_limit(unsigned long);
        void update_time_dependent(Timeloop<TF>&);

        void get_radiation_field(Field3d<TF>&, const std::string&, Thermo<TF>&, Timeloop<TF>&)
        { throw std::runtime_error("\"get_radiation_field()\" is not implemented in radiation_rrtmpg"); }
        std::vector<TF>& get_surface_radiation(const std::string&);
        std::vector<TF>& get_surface_emissivity(const std::string&)
        { throw std::runtime_error("This radiation class cannot provide a surface emissivity field"); }
        std::vector<TF>& get_surface_albedo(const std::string&)
        { throw std::runtime_error("This radiation class cannot provide a surface albedo field"); }

        void exec_all_stats(
                Stats<TF>&, Cross<TF>&, Dump<TF>&, Column<TF>&,
                Thermo<TF>&, Timeloop<TF>&, const unsigned long, const int);
        void exec_individual_column_stats(
                Column<TF>&, Thermo<TF>&, Microphys<TF>&, Timeloop<TF>&, Stats<TF>&,
                Aerosol<TF>&, Background<TF>&);
        void exec_column(Column<TF>&, Thermo<TF>&, Timeloop<TF>&) {};

        #ifdef USECUDA
        TF* get_surface_radiation_g(const std::string&);
        void prepare_device();
        void clear_device();
        void forward_device() {};
        void backward_device() {};
        #endif

    private:
        using Radiation<TF>::swradiation;
        using Radiation<TF>::master;
        using Radiation<TF>::grid;
        using Radiation<TF>::fields;
        using Radiation<TF>::field3d_operators;

        Boundary_cyclic<TF> boundary_cyclic;

        void create_column(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&);
        void create_column_longwave(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&,
                const Gas_concs&);
        void create_column_shortwave(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&,
                const Gas_concs&);

        void read_background_profiles(
                Netcdf_handle&, const Gas_concs&);

        void create_solver(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&, Column<TF>&);
        void create_solver_longwave(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&, Column<TF>&,
                const Gas_concs&);
        void create_solver_shortwave(
                Input&, Netcdf_handle&, Thermo<TF>&, Stats<TF>&, Column<TF>&,
                const Gas_concs&);
        void create_diffuse_filter();


    void solve_shortwave_column(
            std::unique_ptr<Optical_props_arry>&,
            std::unique_ptr<Optical_props_2str>&,
            Array<Float,2>&, Array<Float,2>&,
            Array<Float,2>&, Array<Float,2>&,
            Array<Float,2>&, Array<Float,2>&, const Float,
            const Gas_concs&,
            const Gas_optics_rrtmgp&,
            const Array<Float,2>&,
            const Array<Float,2>&, const Array<Float,2>&,
            const Array<Float,2>&, const Array<Float,2>&,
            Aerosol_concs&,
            const Array<Float,1>&,
            const Array<Float,2>&, const Array<Float,2>&,
            const Float,
            const int
            );

        void solve_longwave_column(
                std::unique_ptr<Optical_props_arry>&,
                Array<Float,2>&, Array<Float,2>&, Array<Float,2>&,
                Array<Float,2>&, const Float,
                const Gas_concs&,
                const std::unique_ptr<Gas_optics_rrtmgp>&,
                const std::unique_ptr<Source_func_lw>&,
                const Array<Float,2>&,
                const Array<Float,2>&, const Array<Float,2>&,
                const Array<Float,2>&, const Array<Float,2>&,
                const Array<Float,1>&, const Array<Float,2>&,
                const int);

        void exec_longwave(
                Thermo<TF>&, Microphys<TF>&, Timeloop<TF>&, Stats<TF>&,
                Array<Float,2>&, Array<Float,2>&, Array<Float,2>&,
                const Array<Float,2>&, const Array<Float,2>&, const Array<Float,1>&,
                const Array<Float,2>&, const Array<Float,2>&, const Array<Float,2>&,
                const bool, const int);

        void exec_shortwave(
                Thermo<TF>&, Microphys<TF>&, Timeloop<TF>&, Stats<TF>&,
                Array<Float,2>&, Array<Float,2>&, Array<Float,2>&, Array<Float,2>&,
                Array<Float,1>&,
                const Array<Float,2>&, const Array<Float,2>&,
                const Array<Float,2>&, const Array<Float,2>&,
                const Array<Float,2>&, const Array<Float,2>&,
                const bool, const int);

        #ifdef USECUDA
        void exec_longwave(
                Thermo<TF>&, Microphys<TF>&, Timeloop<TF>&, Stats<TF>&,
                Array_gpu<Float,2>&, Array_gpu<Float,2>&, Array_gpu<Float,2>&,
                const Array_gpu<Float,2>&, const Array_gpu<Float,2>&, const Array_gpu<Float,1>&,
                const Array_gpu<Float,2>&, const Array_gpu<Float,2>&, const Array_gpu<Float,2>&,
                const bool, const int);

        void exec_shortwave(
                Thermo<TF>&, Microphys<TF>&, Timeloop<TF>&, Stats<TF>&,
                Array_gpu<Float,2>&, Array_gpu<Float,2>&, Array_gpu<Float,2>&, Array_gpu<Float,2>&,
                Array<Float, 1>&,
                const Array_gpu<Float,2>&, const Array_gpu<Float,2>&,
                const Array_gpu<Float,2>&, const Array_gpu<Float,2>&,
                const Array_gpu<Float,2>&, const Array_gpu<Float,2>&,
                const bool, const int);
        #endif

        bool is_day(const Float); // Switch between day/night, based on sza
        void set_sun_location(Timeloop<TF>&);
        void set_background_column_longwave(const TF);
        void set_background_column_shortwave(const TF);

        const std::string tend_name = "rad";
        const std::string tend_longname = "Radiation";

        bool sw_longwave;
        bool sw_shortwave;
        bool sw_clear_sky_stats;
        bool sw_fixed_sza;
        bool sw_aerosol;
        bool sw_delta_cloud;
        bool sw_delta_aer;

        bool swtimedep_background;
        bool swtimedep_aerosol;

        bool sw_homogenize_sfc_sw;
        bool sw_homogenize_sfc_lw;
        bool sw_homogenize_hr_sw;
        bool sw_homogenize_hr_lw;

        // Make sure that the sw radiation is tuned at the first `exec()`. This
        // ensures that sw is tuned for the full 3D field, and not for the column stats.
        bool sw_is_tuned = false;

        double dt_rad;
        unsigned long idt_rad;

        std::vector<std::string> crosslist;

        // RRTMGP related variables.
        Float tsi_scaling; // Total solar irradiance scaling factor.
        Float t_sfc;       // Surface absolute temperature in K.
        Float mu0;         // Cosine of solar zenith angle.
        Float Nc0;         // Total droplet number concentration.

        // The reference column for the full profile.
        Array<Float,2> lw_flux_dn_inc;
        Array<Float,2> sw_flux_dn_dir_inc;
        Array<Float,2> sw_flux_dn_dif_inc;

        int n_col;
        int n_lay_col;
        int n_lev_col;

        Array<Float,2> p_lay_col;
        Array<Float,2> t_lay_col;
        Array<Float,2> p_lev_col;
        Array<Float,2> t_lev_col;
        Array<Float,2> col_dry;

        // Fluxes of reference column
        Array<Float,2> lw_flux_up_col;
        Array<Float,2> lw_flux_dn_col;
        Array<Float,2> lw_flux_net_col;

        Array<Float,2> sw_flux_up_col;
        Array<Float,2> sw_flux_dn_col;
        Array<Float,2> sw_flux_dn_dir_col;
        Array<Float,2> sw_flux_net_col;
        Array<Float,1> aod550;
        int ibnd_550;

        Gas_concs gas_concs_col;

        Aerosol_concs aerosol_concs_col;
        Array<Float,2> rh_col;

        std::unique_ptr<Source_func_lw> sources_lw;
        std::unique_ptr<Optical_props_arry> optical_props_lw;
        std::unique_ptr<Optical_props_arry> optical_props_sw;
        std::unique_ptr<Optical_props_2str> aerosol_props_sw;

        // The full solver.
        Gas_concs gas_concs;
        Aerosol_concs aerosol_concs;

        std::unique_ptr<Gas_optics_rrtmgp> kdist_lw;
        std::unique_ptr<Gas_optics_rrtmgp> kdist_sw;

        std::unique_ptr<Cloud_optics> cloud_lw;
        std::unique_ptr<Cloud_optics> cloud_sw;

        std::unique_ptr<Aerosol_optics> aerosol_sw;

        // Surface fields that go into solver;
        TF emis_sfc_hom;
        TF sfc_alb_dir_hom;
        TF sfc_alb_dif_hom;

        Array<Float,2> emis_sfc;
        Array<Float,2> sfc_alb_dir;
        Array<Float,2> sfc_alb_dif;

        // Surface radiative fluxes CPU
        std::vector<Float> lw_flux_dn_sfc;
        std::vector<Float> lw_flux_up_sfc;

        std::vector<Float> sw_flux_dn_sfc;
        std::vector<Float> sw_flux_up_sfc;

        // Surface radiative fluxes GPU
        Float* lw_flux_dn_sfc_g;
        Float* lw_flux_up_sfc_g;

        Float* sw_flux_dn_sfc_g;
        Float* sw_flux_up_sfc_g;

        Float* sw_flux_dn_dir_inc_g;
        Float* sw_flux_dn_dif_inc_g;
        Float* lw_flux_dn_inc_g;

        // Surface diffuse radiation filtering
        bool sw_diffuse_filter;

        Float sigma_filter;
        Float sigma_filter_small;
        int n_filter_iterations;

        std::vector<Float> sw_flux_dn_dif_f;

        std::vector<Float> filter_kernel_x;
        std::vector<Float> filter_kernel_y;

        // timedependent gases
        std::map<std::string, Timedep<TF>*> tdep_gases;
        std::vector<std::string> gaslist;        ///< List of gases that have timedependent background profiles.
        std::map<std::string, std::vector<TF>> gasprofs; ///< Map of profiles with gases stored by its name.

        #ifdef USECUDA
        std::unique_ptr<Gas_concs_gpu> gas_concs_gpu;
        std::unique_ptr<Aerosol_concs_gpu> aerosol_concs_gpu;
        std::unique_ptr<Gas_optics_gpu> kdist_lw_gpu;
        std::unique_ptr<Cloud_optics_gpu> cloud_lw_gpu;
        std::unique_ptr<Gas_optics_gpu> kdist_sw_gpu;
        std::unique_ptr<Cloud_optics_gpu> cloud_sw_gpu;
        std::unique_ptr<Aerosol_optics_gpu> aerosol_sw_gpu;

        std::map<std::string, TF*> gasprofs_g;    ///< Map of profiles with gasses stored by its name.
        Float* aod550_g;

        Array_gpu<Float,2> emis_sfc_g;
        Array_gpu<Float,2> sfc_alb_dir_g;
        Array_gpu<Float,2> sfc_alb_dif_g;

        Rte_lw_gpu rte_lw_gpu;
        Rte_sw_gpu rte_sw_gpu;
        #endif
};
#endif
