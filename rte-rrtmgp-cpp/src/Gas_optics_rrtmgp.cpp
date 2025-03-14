/*
 * This file is part of a C++ interface to the Radiative Transfer for Energetics (RTE)
 * and Rapid Radiative Transfer Model for GCM applications Parallel (RRTMGP).
 *
 * The original code is found at https://github.com/earth-system-radiation/rte-rrtmgp.
 *
 * Contacts: Robert Pincus and Eli Mlawer
 * email: rrtmgp@aer.com
 *
 * Copyright 2015-2020,  Atmospheric and Environmental Research and
 * Regents of the University of Colorado.  All right reserved.
 *
 * This C++ interface can be downloaded from https://github.com/earth-system-radiation/rte-rrtmgp-cpp
 *
 * Contact: Chiel van Heerwaarden
 * email: chiel.vanheerwaarden@wur.nl
 *
 * Copyright 2020, Wageningen University & Research.
 *
 * Use and duplication is permitted under the terms of the
 * BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
 *
 */

#include <chrono>
#include <numeric>
#include <iomanip>
#include <cmath>
#include <boost/algorithm/string.hpp>

#include "Gas_concs.h"
#include "Gas_optics_rrtmgp.h"
#include "Array.h"
#include "Optical_props.h"
#include "Source_functions.h"

#include "rrtmgp_kernels.h"

#define restrict __restrict__

namespace
{
    int find_index(
            const Array<std::string,1>& data, const std::string& value)
    {
        auto it = std::find(data.v().begin(), data.v().end(), value);
        if (it == data.v().end())
            return -1;
        else
            return it - data.v().begin() + 1;
    }

    template<typename Float>
    void reduce_minor_arrays(
                const Gas_concs& available_gases,
                const Array<std::string,1>& gas_names,
                const Array<std::string,1>& gas_minor,
                const Array<std::string,1>& identifier_minor,
                const Array<Float,3>& kminor_atm,
                const Array<std::string,1>& minor_gases_atm,
                const Array<int,2>& minor_limits_gpt_atm,
                const Array<Bool,1>& minor_scales_with_density_atm,
                const Array<std::string,1>& scaling_gas_atm,
                const Array<Bool,1>& scale_by_complement_atm,
                const Array<int,1>& kminor_start_atm,

                Array<Float,3>& kminor_atm_red,
                Array<std::string,1>& minor_gases_atm_red,
                Array<int,2>& minor_limits_gpt_atm_red,
                Array<Bool,1>& minor_scales_with_density_atm_red,
                Array<std::string,1>& scaling_gas_atm_red,
                Array<Bool,1>& scale_by_complement_atm_red,
                Array<int,1>& kminor_start_atm_red)
    {
        int nm = minor_gases_atm.dim(1);
        int tot_g = 0;

        Array<Bool,1> gas_is_present({nm});

        for (int i=1; i<=nm; ++i)
        {
            const int idx_mnr = find_index(identifier_minor, minor_gases_atm({i}));

            // Search for
            std::string gas_minor_trimmed = gas_minor({idx_mnr});
            boost::trim(gas_minor_trimmed);

            gas_is_present({i}) = available_gases.exists(gas_minor_trimmed);
            if (gas_is_present({i}))
                tot_g += minor_limits_gpt_atm({2,i}) - minor_limits_gpt_atm({1,i}) + 1;
        }

        const int red_nm = std::accumulate(gas_is_present.v().begin(), gas_is_present.v().end(), 0);

        Array<Float,3> kminor_atm_red_t;

        if (red_nm == nm)
        {
            kminor_atm_red_t = kminor_atm;
            minor_gases_atm_red = minor_gases_atm;
            minor_limits_gpt_atm_red = minor_limits_gpt_atm;
            minor_scales_with_density_atm_red = minor_scales_with_density_atm;
            scaling_gas_atm_red = scaling_gas_atm;
            scale_by_complement_atm_red = scale_by_complement_atm;
            kminor_start_atm_red = kminor_start_atm;
        }
        else
        {
            // Use a lambda function as the operation has to be repeated many times.
            auto resize_and_set = [&](auto& a_red, const auto& a)
            {
                a_red.set_dims({red_nm});
                int counter = 1;
                for (int i=1; i<=gas_is_present.dim(1); ++i)
                {
                    if (gas_is_present({i}))
                    {
                       a_red({counter}) = a({i});
                       ++counter;
                    }
                }
            };

            resize_and_set(minor_gases_atm_red, minor_gases_atm);
            resize_and_set(minor_scales_with_density_atm_red, minor_scales_with_density_atm);
            resize_and_set(scaling_gas_atm_red, scaling_gas_atm);
            resize_and_set(scale_by_complement_atm_red, scale_by_complement_atm);
            resize_and_set(kminor_start_atm_red, kminor_start_atm);

            minor_limits_gpt_atm_red.set_dims({2, red_nm});
            kminor_atm_red_t.set_dims({tot_g, kminor_atm.dim(2), kminor_atm.dim(3)});

            int icnt = 0;
            int n_elim = 0;
            for (int i=1; i<=nm; ++i)
            {
                int ng = minor_limits_gpt_atm({2,i}) - minor_limits_gpt_atm({1,i}) + 1;
                if (gas_is_present({i}))
                {
                    ++icnt;
                    minor_limits_gpt_atm_red({1,icnt}) = minor_limits_gpt_atm({1,i});
                    minor_limits_gpt_atm_red({2,icnt}) = minor_limits_gpt_atm({2,i});
                    kminor_start_atm_red({icnt}) = kminor_start_atm({i}) - n_elim;

                    for (int j=1; j<=ng; ++j)
                        for (int i2=1; i2<=kminor_atm.dim(2); ++i2)
                            for (int i3=1; i3<=kminor_atm.dim(3); ++i3)
                                kminor_atm_red_t({kminor_start_atm_red({icnt})+j-1,i2,i3}) =
                                        kminor_atm({kminor_start_atm({i})+j-1,i2,i3});
                }
                else
                    n_elim += ng;
            }
        }

        // Reshape following the new ordering in v1.5.
        kminor_atm_red.set_dims({kminor_atm_red_t.dim(3), kminor_atm_red_t.dim(2), kminor_atm_red_t.dim(1)});
        for (int i3=1; i3<=kminor_atm_red.dim(3); ++i3)
            for (int i2=1; i2<=kminor_atm_red.dim(2); ++i2)
                for (int i1=1; i1<=kminor_atm_red.dim(1); ++i1)
                    kminor_atm_red({i1, i2, i3}) = kminor_atm_red_t({i3, i2, i1});
    }

    void create_idx_minor(
            const Array<std::string,1>& gas_names,
            const Array<std::string,1>& gas_minor,
            const Array<std::string,1>& identifier_minor,
            const Array<std::string,1>& minor_gases_atm,
            Array<int,1>& idx_minor_atm)
    {
        Array<int,1> idx_minor_atm_out({minor_gases_atm.dim(1)});

        for (int imnr=1; imnr<=minor_gases_atm.dim(1); ++imnr)
        {
            // Find identifying string for minor species in list of possible identifiers (e.g. h2o_slf)
            const int idx_mnr = find_index(identifier_minor, minor_gases_atm({imnr}));

            // Find name of gas associated with minor species identifier (e.g. h2o)
            idx_minor_atm_out({imnr}) = find_index(gas_names, gas_minor({idx_mnr}));
        }

        idx_minor_atm = idx_minor_atm_out;
    }

    void create_idx_minor_scaling(
            const Array<std::string,1>& gas_names,
            const Array<std::string,1>& scaling_gas_atm,
            Array<int,1>& idx_minor_scaling_atm)
    {
        Array<int,1> idx_minor_scaling_atm_out({scaling_gas_atm.dim(1)});

        for (int imnr=1; imnr<=scaling_gas_atm.dim(1); ++imnr)
            idx_minor_scaling_atm_out({imnr}) = find_index(gas_names, scaling_gas_atm({imnr}));

        idx_minor_scaling_atm = idx_minor_scaling_atm_out;
    }

    void create_key_species_reduce(
            const Array<std::string,1>& gas_names,
            const Array<std::string,1>& gas_names_red,
            const Array<int,3>& key_species,
            Array<int,3>& key_species_red,
            Array<Bool,1>& key_species_present_init)
    {
        const int np = key_species.dim(1);
        const int na = key_species.dim(2);
        const int nt = key_species.dim(3);

        key_species_red.set_dims({key_species.dim(1), key_species.dim(2), key_species.dim(3)});
        key_species_present_init.set_dims({gas_names.dim(1)});

        for (int i=1; i<=key_species_present_init.dim(1); ++i)
            key_species_present_init({i}) = 1;

        for (int ip=1; ip<=np; ++ip)
            for (int ia=1; ia<=na; ++ia)
                for (int it=1; it<=nt; ++it)
                {
                    const int ks = key_species({ip,ia,it});
                    if (ks != 0)
                    {
                        const int ksr = find_index(gas_names_red, gas_names({ks}));
                        key_species_red({ip,ia,it}) = ksr;
                        if (ksr == -1)
                            key_species_present_init({ks}) = 0;
                    }
                    else
                        key_species_red({ip,ia,it}) = ks;
                }
    }

    void check_key_species_present_init(
            const Array<std::string,1>& gas_names,
            const Array<Bool,1>& key_species_present_init
            )
    {
        for (int i=1; i<=key_species_present_init.dim(1); ++i)
        {
            if (key_species_present_init({i}) == 0)
            {
                std::string error_message = "Gas optics: required gas " + gas_names({i}) + " is missing";
                throw std::runtime_error(error_message);
            }
        }
    }

    void create_flavor(
            const Array<int,3>& key_species,
            Array<int,2>& flavor)
    {
        Array<int,2> key_species_list({2, key_species.dim(3)*2});

        // Prepare list of key species.
        int i = 1;
        for (int ibnd=1; ibnd<=key_species.dim(3); ++ibnd)
            for (int iatm=1; iatm<=key_species.dim(2); ++iatm)
            {
                key_species_list({1,i}) = key_species({1,iatm,ibnd});
                key_species_list({2,i}) = key_species({2,iatm,ibnd});
                ++i;
            }

        // Rewrite single key_species pairs.
        for (int i=1; i<=key_species_list.dim(2); ++i)
        {
            if ( key_species_list({1,i}) == 0 && key_species_list({2,i}) == 0 )
            {
                key_species_list({1,i}) = 2;
                key_species_list({2,i}) = 2;
            }
        }

        // Count unique key species pairs.
        int iflavor = 0;
        for (int i=1; i<=key_species_list.dim(2); ++i)
        {
            bool pair_exists = false;
            for (int ii=1; ii<=i-1; ++ii)
            {
                if ( (key_species_list({1,i}) == key_species_list({1,ii})) &&
                     (key_species_list({2,i}) == key_species_list({2,ii})) )
                {
                    pair_exists = true;
                    break;
                }
            }
            if (!pair_exists)
                ++iflavor;
        }

        // Fill flavors.
        flavor.set_dims({2,iflavor});
        iflavor = 0;
        for (int i=1; i<=key_species_list.dim(2); ++i)
        {
            bool pair_exists = false;
            for (int ii=1; ii<=i-1; ++ii)
            {
                if ( (key_species_list({1,i}) == key_species_list({1,ii})) &&
                     (key_species_list({2,i}) == key_species_list({2,ii})) )
                {
                    pair_exists = true;
                    break;
                }
            }
            if (!pair_exists)
            {
                ++iflavor;
                flavor({1,iflavor}) = key_species_list({1,i});
                flavor({2,iflavor}) = key_species_list({2,i});
            }
        }
    }

    int key_species_pair2flavor(
            const Array<int,2>& flavor,
            const Array<int,1>& key_species_pair)
    {
        // Search for match.
        for (int iflav=1; iflav<=flavor.dim(2); ++iflav)
        {
            if ( key_species_pair({1}) == flavor({1, iflav}) &&
                 key_species_pair({2}) == flavor({2, iflav}) )
                return iflav;
        }

        // No match found.
        return -1;
    }

    void create_gpoint_flavor(
            const Array<int,3>& key_species,
            const Array<int,1>& gpt2band,
            const Array<int,2>& flavor,
            Array<int,2>& gpoint_flavor)
    {
        const int ngpt = gpt2band.dim(1);
        gpoint_flavor.set_dims({2,ngpt});

        for (int igpt=1; igpt<=ngpt; ++igpt)
            for (int iatm=1; iatm<=2; ++iatm)
            {
                int pair_1 = key_species( {1, iatm, gpt2band({igpt})} );
                int pair_2 = key_species( {2, iatm, gpt2band({igpt})} );

                // Rewrite species pair.
                Array<int,1> rewritten_pair({2});
                if (pair_1 == 0 && pair_2 == 0)
                {
                    rewritten_pair({1}) = 2;
                    rewritten_pair({2}) = 2;
                }
                else
                {
                    rewritten_pair({1}) = pair_1;
                    rewritten_pair({2}) = pair_2;
                }

                // Write the output.
                gpoint_flavor({iatm,igpt}) = key_species_pair2flavor(
                        flavor, rewritten_pair);
            }
    }

    template<typename Float>
    void combine_abs_and_rayleigh_kernel(
            const int ncol, const int nlay, const int ngpt,
            const Float* __restrict__ tau, const Float* __restrict__ tau_rayleigh,
            Float* __restrict__ tau_out, Float* __restrict__ ssa_out)
    {
        for (int igpt=0; igpt<ngpt; ++igpt)
            for (int ilay=0; ilay<nlay; ++ilay)
                for (int icol=0; icol<ncol; ++icol)
                {
                    const int idx = icol + ilay*ncol + igpt*ncol*nlay;
                    const Float t = tau[idx] + tau_rayleigh[idx];
    
                    if (t > Float(2.) * std::numeric_limits<Float>::epsilon())
                        ssa_out[idx] = tau_rayleigh[idx] / t;
                    else
                        ssa_out[idx] = Float(0.);
    
                    tau_out[idx] = t;
                }
    }
}


// IMPLEMENTATION OF CLASS FUNCTIONS.
// Constructor of longwave variant.
Gas_optics_rrtmgp::Gas_optics_rrtmgp(
        const Gas_concs& available_gases,
        const Array<std::string,1>& gas_names,
        const Array<int,3>& key_species,
        const Array<int,2>& band2gpt,
        const Array<Float,2>& band_lims_wavenum,
        const Array<Float,1>& press_ref,
        const Float press_ref_trop,
        const Array<Float,1>& temp_ref,
        const Float temp_ref_p,
        const Float temp_ref_t,
        const Array<Float,3>& vmr_ref,
        const Array<Float,4>& kmajor,
        const Array<Float,3>& kminor_lower,
        const Array<Float,3>& kminor_upper,
        const Array<std::string,1>& gas_minor,
        const Array<std::string,1>& identifier_minor,
        const Array<std::string,1>& minor_gases_lower,
        const Array<std::string,1>& minor_gases_upper,
        const Array<int,2>& minor_limits_gpt_lower,
        const Array<int,2>& minor_limits_gpt_upper,
        const Array<Bool,1>& minor_scales_with_density_lower,
        const Array<Bool,1>& minor_scales_with_density_upper,
        const Array<std::string,1>& scaling_gas_lower,
        const Array<std::string,1>& scaling_gas_upper,
        const Array<Bool,1>& scale_by_complement_lower,
        const Array<Bool,1>& scale_by_complement_upper,
        const Array<int,1>& kminor_start_lower,
        const Array<int,1>& kminor_start_upper,
        const Array<Float,2>& totplnk,
        const Array<Float,4>& planck_frac,
        const Array<Float,3>& rayl_lower,
        const Array<Float,3>& rayl_upper) :
            Gas_optics(band_lims_wavenum, band2gpt),
            totplnk(totplnk)
{
    // Reshaping according to new dimension ordering since v1.5
    this->planck_frac.set_dims({planck_frac.dim(4), planck_frac.dim(2), planck_frac.dim(3), planck_frac.dim(1)});
    for (int i4=1; i4<=this->planck_frac.dim(4); ++i4)
        for (int i3=1; i3<=this->planck_frac.dim(3); ++i3)
            for (int i2=1; i2<=this->planck_frac.dim(2); ++i2)
                for (int i1=1; i1<=this->planck_frac.dim(1); ++i1)
                    this->planck_frac({i1, i2, i3, i4}) = planck_frac({i4, i2, i3, i1});

    // Initialize the absorption coefficient array, including Rayleigh scattering
    // tables if provided.
    init_abs_coeffs(
            available_gases,
            gas_names, key_species,
            band2gpt, band_lims_wavenum,
            press_ref, temp_ref,
            press_ref_trop, temp_ref_p, temp_ref_t,
            vmr_ref,
            kmajor, kminor_lower, kminor_upper,
            gas_minor,identifier_minor,
            minor_gases_lower, minor_gases_upper,
            minor_limits_gpt_lower,
            minor_limits_gpt_upper,
            minor_scales_with_density_lower,
            minor_scales_with_density_upper,
            scaling_gas_lower, scaling_gas_upper,
            scale_by_complement_lower,
            scale_by_complement_upper,
            kminor_start_lower,
            kminor_start_upper,
            rayl_lower, rayl_upper);

    // Temperature steps for Planck function interpolation.
    // Assumes that temperature minimum and max are the same for the absorption coefficient grid and the
    // Planck grid and the Planck grid is equally spaced.
    totplnk_delta = (temp_ref_max - temp_ref_min) / (totplnk.dim(1)-1);
}


// Constructor of the shortwave variant.
Gas_optics_rrtmgp::Gas_optics_rrtmgp(
        const Gas_concs& available_gases,
        const Array<std::string,1>& gas_names,
        const Array<int,3>& key_species,
        const Array<int,2>& band2gpt,
        const Array<Float,2>& band_lims_wavenum,
        const Array<Float,1>& press_ref,
        const Float press_ref_trop,
        const Array<Float,1>& temp_ref,
        const Float temp_ref_p,
        const Float temp_ref_t,
        const Array<Float,3>& vmr_ref,
        const Array<Float,4>& kmajor,
        const Array<Float,3>& kminor_lower,
        const Array<Float,3>& kminor_upper,
        const Array<std::string,1>& gas_minor,
        const Array<std::string,1>& identifier_minor,
        const Array<std::string,1>& minor_gases_lower,
        const Array<std::string,1>& minor_gases_upper,
        const Array<int,2>& minor_limits_gpt_lower,
        const Array<int,2>& minor_limits_gpt_upper,
        const Array<Bool,1>& minor_scales_with_density_lower,
        const Array<Bool,1>& minor_scales_with_density_upper,
        const Array<std::string,1>& scaling_gas_lower,
        const Array<std::string,1>& scaling_gas_upper,
        const Array<Bool,1>& scale_by_complement_lower,
        const Array<Bool,1>& scale_by_complement_upper,
        const Array<int,1>& kminor_start_lower,
        const Array<int,1>& kminor_start_upper,
        const Array<Float,1>& solar_source_quiet,
        const Array<Float,1>& solar_source_facular,
        const Array<Float,1>& solar_source_sunspot,
        const Float tsi_default,
        const Float mg_default,
        const Float sb_default,
        const Array<Float,3>& rayl_lower,
        const Array<Float,3>& rayl_upper) :
            Gas_optics(band_lims_wavenum, band2gpt)
{
    // Initialize the absorption coefficient array, including Rayleigh scattering
    // tables if provided.
    init_abs_coeffs(
            available_gases,
            gas_names, key_species,
            band2gpt, band_lims_wavenum,
            press_ref, temp_ref,
            press_ref_trop, temp_ref_p, temp_ref_t,
            vmr_ref,
            kmajor, kminor_lower, kminor_upper,
            gas_minor,identifier_minor,
            minor_gases_lower, minor_gases_upper,
            minor_limits_gpt_lower,
            minor_limits_gpt_upper,
            minor_scales_with_density_lower,
            minor_scales_with_density_upper,
            scaling_gas_lower, scaling_gas_upper,
            scale_by_complement_lower,
            scale_by_complement_upper,
            kminor_start_lower,
            kminor_start_upper,
            rayl_lower, rayl_upper);

    // Compute the solar source.
    this->solar_source_quiet = solar_source_quiet;
    this->solar_source_facular = solar_source_facular;
    this->solar_source_sunspot = solar_source_sunspot;

    this->solar_source.set_dims(solar_source_quiet.get_dims());

    set_solar_variability(mg_default, sb_default);
}


void Gas_optics_rrtmgp::init_abs_coeffs(
        const Gas_concs& available_gases,
        const Array<std::string,1>& gas_names,
        const Array<int,3>& key_species,
        const Array<int,2>& band2gpt,
        const Array<Float,2>& band_lims_wavenum,
        const Array<Float,1>& press_ref,
        const Array<Float,1>& temp_ref,
        const Float press_ref_trop,
        const Float temp_ref_p,
        const Float temp_ref_t,
        const Array<Float,3>& vmr_ref,
        const Array<Float,4>& kmajor,
        const Array<Float,3>& kminor_lower,
        const Array<Float,3>& kminor_upper,
        const Array<std::string,1>& gas_minor,
        const Array<std::string,1>& identifier_minor,
        const Array<std::string,1>& minor_gases_lower,
        const Array<std::string,1>& minor_gases_upper,
        const Array<int,2>& minor_limits_gpt_lower,
        const Array<int,2>& minor_limits_gpt_upper,
        const Array<Bool,1>& minor_scales_with_density_lower,
        const Array<Bool,1>& minor_scales_with_density_upper,
        const Array<std::string,1>& scaling_gas_lower,
        const Array<std::string,1>& scaling_gas_upper,
        const Array<Bool,1>& scale_by_complement_lower,
        const Array<Bool,1>& scale_by_complement_upper,
        const Array<int,1>& kminor_start_lower,
        const Array<int,1>& kminor_start_upper,
        const Array<Float,3>& rayl_lower,
        const Array<Float,3>& rayl_upper)
{
    // Which gases known to the gas optics are present in the host model (available_gases)?
    std::vector<std::string> gas_names_to_use;

    for (const std::string &s : gas_names.v())
    {
        if (available_gases.exists(s))
            gas_names_to_use.push_back(s);
    }

    // Now the number of gases is the union of those known to the k-distribution and provided
    // by the host model.
    const int n_gas = gas_names_to_use.size();
    Array<std::string,1> gas_names_this(std::move(gas_names_to_use), {n_gas});
    this->gas_names = gas_names_this;

    // Initialize the gas optics object, keeping only those gases known to the
    // gas optics and also present in the host model.
    // Add an offset to the indexing to interface the negative ranging of fortran.
    Array<Float,3> vmr_ref_red({vmr_ref.dim(1), n_gas + 1, vmr_ref.dim(3)});
    vmr_ref_red.set_offsets({0, -1, 0});

    // Gas 0 is used in single-key species method, set to 1.0 (col_dry)
    for (int i1=1; i1<=vmr_ref_red.dim(1); ++i1)
        for (int i3=1; i3<=vmr_ref_red.dim(3); ++i3)
            vmr_ref_red({i1, 0, i3}) = vmr_ref({i1, 1, i3});

    for (int i=1; i<=n_gas; ++i)
    {
        int idx = find_index(gas_names, this->gas_names({i}));
        for (int i1=1; i1<=vmr_ref_red.dim(1); ++i1)
            for (int i3=1; i3<=vmr_ref_red.dim(3); ++i3)
                vmr_ref_red({i1, i, i3}) = vmr_ref({i1, idx+1, i3}); // CvH: why +1?
    }

    this->vmr_ref = std::move(vmr_ref_red);

    // Reduce minor arrays so variables only contain minor gases that are available.
    // Reduce size of minor Arrays.
    Array<std::string, 1> minor_gases_lower_red;
    Array<std::string, 1> scaling_gas_lower_red;
    Array<std::string, 1> minor_gases_upper_red;
    Array<std::string, 1> scaling_gas_upper_red;

    reduce_minor_arrays(
            available_gases,
            gas_names,
            gas_minor, identifier_minor,
            kminor_lower,
            minor_gases_lower,
            minor_limits_gpt_lower,
            minor_scales_with_density_lower,
            scaling_gas_lower,
            scale_by_complement_lower,
            kminor_start_lower,
            this->kminor_lower,
            minor_gases_lower_red,
            this->minor_limits_gpt_lower,
            this->minor_scales_with_density_lower,
            scaling_gas_lower_red,
            this->scale_by_complement_lower,
            this->kminor_start_lower);

    reduce_minor_arrays(
            available_gases,
            gas_names,
            gas_minor,
            identifier_minor,
            kminor_upper,
            minor_gases_upper,
            minor_limits_gpt_upper,
            minor_scales_with_density_upper,
            scaling_gas_upper,
            scale_by_complement_upper,
            kminor_start_upper,
            this->kminor_upper,
            minor_gases_upper_red,
            this->minor_limits_gpt_upper,
            this->minor_scales_with_density_upper,
            scaling_gas_upper_red,
            this->scale_by_complement_upper,
            this->kminor_start_upper);

    // Arrays not reduced by the presence, or lack thereof, of a gas
    this->press_ref = press_ref;
    this->temp_ref = temp_ref;

    // Reshaping according to new dimension ordering since v1.5
    this->kmajor.set_dims({kmajor.dim(4), kmajor.dim(2), kmajor.dim(3), kmajor.dim(1)});
    for (int i4=1; i4<=this->kmajor.dim(4); ++i4)
        for (int i3=1; i3<=this->kmajor.dim(3); ++i3)
            for (int i2=1; i2<=this->kmajor.dim(2); ++i2)
                for (int i1=1; i1<=this->kmajor.dim(1); ++i1)
                    this->kmajor({i1, i2, i3, i4}) = kmajor({i4, i2, i3, i1});

    // Reshaping according to new 1.5 release.
    // Create a new vector that consists of rayl_lower and rayl_upper stored in one variable.
    if (rayl_lower.size() > 0)
    {
        this->krayl.set_dims({rayl_lower.dim(3), rayl_lower.dim(2), rayl_lower.dim(1), 2});
        for (int i3=1; i3<=this->krayl.dim(3); ++i3)
            for (int i2=1; i2<=this->krayl.dim(2); ++i2)
                for (int i1=1; i1<=this->krayl.dim(1); ++i1)
                {
                    this->krayl({i1, i2, i3, 1}) = rayl_lower({i3, i2, i1});
                    this->krayl({i1, i2, i3, 2}) = rayl_upper({i3, i2, i1});
                }
    }

    // ---- post processing ----
    //  creates log reference pressure
    this->press_ref_log = this->press_ref;
    for (int i1=1; i1<=this->press_ref_log.dim(1); ++i1)
        this->press_ref_log({i1}) = std::log(this->press_ref_log({i1}));

    // log scale of reference pressure
    this->press_ref_trop_log = std::log(press_ref_trop);

    // Get index of gas (if present) for determining col_gas
    create_idx_minor(
            this->gas_names, gas_minor, identifier_minor, minor_gases_lower_red, this->idx_minor_lower);
    create_idx_minor(
            this->gas_names, gas_minor, identifier_minor, minor_gases_upper_red, this->idx_minor_upper);

    // Get index of gas (if present) that has special treatment in density scaling
    create_idx_minor_scaling(
            this->gas_names, scaling_gas_lower_red, this->idx_minor_scaling_lower);
    create_idx_minor_scaling(
            this->gas_names, scaling_gas_upper_red, this->idx_minor_scaling_upper);

    // Create flavor list.
    // Reduce (remap) key_species list; checks that all key gases are present in incoming
    Array<int,3> key_species_red;
    Array<Bool,1> key_species_present_init;

    create_key_species_reduce(
            gas_names, this->gas_names, key_species, key_species_red, key_species_present_init);

    check_key_species_present_init(gas_names, key_species_present_init);

    // create flavor list
    create_flavor(key_species_red, this->flavor);

    // create gpoint flavor list
    create_gpoint_flavor(
            key_species_red, this->get_gpoint_bands(), this->flavor, this->gpoint_flavor);

    // minimum, maximum reference temperature, pressure -- assumes low-to-high ordering
    // for T, high-to-low ordering for p
    this->temp_ref_min = this->temp_ref({1});
    this->temp_ref_max = this->temp_ref({temp_ref.dim(1)});
    this->press_ref_min = this->press_ref({press_ref.dim(1)});
    this->press_ref_max = this->press_ref({1});

    // creates press_ref_log, temp_ref_delta
    this->press_ref_log_delta =
            (std::log(this->press_ref_min) - std::log(this->press_ref_max)) / (this->press_ref.dim(1) - 1);
    this->temp_ref_delta = (this->temp_ref_max - this->temp_ref_min) / (this->temp_ref.dim(1) - 1);

    // Which species are key in one or more bands?
    // this->flavor is an index into this->gas_names
    // if (allocated(this%is_key)) deallocate(this%is_key) ! Shouldn't ever happen...
    Array<int,1> is_key({get_ngas()}); // CvH bool, defaults to 0.?

    for (int j=1; j<=this->flavor.dim(2); ++j)
        for (int i=1; i<=this->flavor.dim(1); ++i)
        {
            if (this->flavor({i, j}) != 0)
                is_key({this->flavor({i, j})}) = true;
        }

    this->is_key = is_key;
}


void Gas_optics_rrtmgp::set_solar_variability(
        const Float mg_index, const Float sb_index)
{
    constexpr Float a_offset = Float(0.1495954);
    constexpr Float b_offset = Float(0.00066696);

    for (int igpt=1; igpt<=this->solar_source_quiet.dim(1); ++igpt)
    {
        this->solar_source({igpt}) = this->solar_source_quiet({igpt})
                + (mg_index - a_offset) * this->solar_source_facular({igpt})
                + (sb_index - b_offset) * this->solar_source_sunspot({igpt});
    }

    // Scale solar source to input TSI value
    // if (present(tsi)) error_msg = this%set_tsi(tsi)
}


// Calculate the molecules of dry air.
void Gas_optics_rrtmgp::get_col_dry(
        Array<Float,2>& col_dry, const Array<Float,2>& vmr_h2o,
        const Array<Float,2>& plev)
{
    // CvH: RRTMGP uses more accurate method based on latitude.
    constexpr Float g0 = 9.80665;

    constexpr Float avogad = 6.02214076e23;
    constexpr Float m_dry = 0.028964;
    constexpr Float m_h2o = 0.018016;

    Array<Float,2> delta_plev({col_dry.dim(1), col_dry.dim(2)});
    Array<Float,2> m_air     ({col_dry.dim(1), col_dry.dim(2)});

    for (int ilay=1; ilay<=col_dry.dim(2); ++ilay)
        for (int icol=1; icol<=col_dry.dim(1); ++icol)
            delta_plev({icol, ilay}) = std::abs(plev({icol, ilay}) - plev({icol, ilay+1}));

    for (int ilay=1; ilay<=col_dry.dim(2); ++ilay)
        for (int icol=1; icol<=col_dry.dim(1); ++icol)
            m_air({icol, ilay}) = (m_dry + m_h2o * vmr_h2o({icol, ilay})) / (1. + vmr_h2o({icol, ilay}));

    for (int ilay=1; ilay<=col_dry.dim(2); ++ilay)
        for (int icol=1; icol<=col_dry.dim(1); ++icol)
        {
            col_dry({icol, ilay}) = Float(10.) * delta_plev({icol, ilay}) * avogad / (Float(1000.)*m_air({icol, ilay})*Float(100.)*g0);
            col_dry({icol, ilay}) /= (Float(1.) + vmr_h2o({icol, ilay}));
        }
}


// Gas optics solver longwave variant.
void Gas_optics_rrtmgp::gas_optics(
        const Array<Float,2>& play,
        const Array<Float,2>& plev,
        const Array<Float,2>& tlay,
        const Array<Float,1>& tsfc,
        const Gas_concs& gas_desc,
        std::unique_ptr<Optical_props_arry>& optical_props,
        Source_func_lw& sources,
        const Array<Float,2>& col_dry,
        const Array<Float,2>& tlev) const
{
    const int ncol = play.dim(1);
    const int nlay = play.dim(2);
    const int ngpt = this->get_ngpt();
    const int nband = this->get_nband();

    // Check if any of the values is out of range.
    if (any_vals_outside(play, this->press_ref_min, this->press_ref_max))
        throw std::range_error("play is out of range");
    if (any_vals_outside(plev, this->press_ref_min, this->press_ref_max))
        throw std::range_error("plev is out of range");

    if (any_vals_outside(tlay, this->temp_ref_min, this->temp_ref_max))
        throw std::range_error("tlay is out of range");
    if (any_vals_outside(tlev, this->temp_ref_min, this->temp_ref_max))
        throw std::range_error("tlev is out of range");
    if (any_vals_outside(tsfc, this->temp_ref_min, this->temp_ref_max))
        throw std::range_error("tsfc is out of range");

    if (any_vals_less_than(col_dry, Float(0.)))
        throw std::range_error("col_dry is out of range");
    // End of checks.

    Array<int,2> jtemp({play.dim(1), play.dim(2)});
    Array<int,2> jpress({play.dim(1), play.dim(2)});
    Array<Bool,2> tropo({play.dim(1), play.dim(2)});
    Array<Float,6> fmajor({2, 2, 2, play.dim(1), play.dim(2), this->get_nflav()});
    Array<int,4> jeta({2, play.dim(1), play.dim(2), this->get_nflav()});

    // Gas optics.
    compute_gas_taus(
            ncol, nlay, ngpt, nband,
            play, plev, tlay, gas_desc,
            optical_props,
            jtemp, jpress, jeta, tropo, fmajor,
            col_dry);

    // External sources.
    source(
            ncol, nlay, nband, ngpt,
            play, plev, tlay, tsfc,
            jtemp, jpress, jeta, tropo, fmajor,
            sources, tlev);
}


// Gas optics solver shortwave variant.
void Gas_optics_rrtmgp::gas_optics(
        const Array<Float,2>& play,
        const Array<Float,2>& plev,
        const Array<Float,2>& tlay,
        const Gas_concs& gas_desc,
        std::unique_ptr<Optical_props_arry>& optical_props,
        Array<Float,2>& toa_src,
        const Array<Float,2>& col_dry) const
{
    const int ncol = play.dim(1);
    const int nlay = play.dim(2);
    const int ngpt = this->get_ngpt();
    const int nband = this->get_nband();

    // Check if any of the values is out of range.
    if (any_vals_outside(play, this->press_ref_min, this->press_ref_max))
        throw std::range_error("play is out of range");
    if (any_vals_outside(plev, this->press_ref_min, this->press_ref_max))
        throw std::range_error("plev is out of range");

    if (any_vals_outside(tlay, this->temp_ref_min, this->temp_ref_max))
        throw std::range_error("tlay is out of range");

    if (any_vals_less_than(col_dry, Float(0.)))
        throw std::range_error("col_dry is out of range");
    // End of checks.

    Array<int,2> jtemp({play.dim(1), play.dim(2)});
    Array<int,2> jpress({play.dim(1), play.dim(2)});
    Array<Bool,2> tropo({play.dim(1), play.dim(2)});
    Array<Float,6> fmajor({2, 2, 2, play.dim(1), play.dim(2), this->get_nflav()});
    Array<int,4> jeta({2, play.dim(1), play.dim(2), this->get_nflav()});

    // Gas optics.
    compute_gas_taus(
            ncol, nlay, ngpt, nband,
            play, plev, tlay, gas_desc,
            optical_props,
            jtemp, jpress, jeta, tropo, fmajor,
            col_dry);

    // External source function is constant.
    for (int igpt=1; igpt<=ngpt; ++igpt)
        for (int icol=1; icol<=ncol; ++icol)
            toa_src({icol, igpt}) = this->solar_source({igpt});
}


namespace rrtmgp_kernel_launcher
{
    template<typename Float> void zero_array(
            int ni, int nj, int nk, Array<Float,3>& array)
    {
        rrtmgp_kernels::zero_array_3D(&ni, &nj, &nk, array.ptr());
    }

    template<typename Float> void zero_array(
            int ni, int nj, int nk, int nl, Array<Float,4>& array)
    {
        rrtmgp_kernels::zero_array_4D(&ni, &nj, &nk, &nl, array.ptr());
    }

    template<typename Float>
    void interpolation(
            int ncol, int nlay,
            int ngas, int nflav, int neta, int npres, int ntemp,
            const Array<int,2>& flavor,
            const Array<Float,1>& press_ref_log,
            const Array<Float,1>& temp_ref,
            Float press_ref_log_delta,
            Float temp_ref_min,
            Float temp_ref_delta,
            Float press_ref_trop_log,
            const Array<Float,3>& vmr_ref,
            const Array<Float,2>& play,
            const Array<Float,2>& tlay,
            Array<Float,3>& col_gas,
            Array<int,2>& jtemp,
            Array<Float,6>& fmajor, Array<Float,5>& fminor,
            Array<Float,4>& col_mix,
            Array<Bool,2>& tropo,
            Array<int,4>& jeta,
            Array<int,2>& jpress)
    {
        rrtmgp_kernels::rrtmgp_interpolation(
                &ncol, &nlay,
                &ngas, &nflav, &neta, &npres, &ntemp,
                const_cast<int*>(flavor.ptr()),
                const_cast<Float*>(press_ref_log.ptr()),
                const_cast<Float*>(temp_ref.ptr()),
                &press_ref_log_delta,
                &temp_ref_min,
                &temp_ref_delta,
                &press_ref_trop_log,
                const_cast<Float*>(vmr_ref.ptr()),
                const_cast<Float*>(play.ptr()),
                const_cast<Float*>(tlay.ptr()),
                col_gas.ptr(),
                jtemp.ptr(),
                fmajor.ptr(), fminor.ptr(),
                col_mix.ptr(),
                tropo.ptr(),
                jeta.ptr(),
                jpress.ptr());
    }

    template<typename Float>
    void compute_tau_absorption(
            int ncol, int nlay, int nband, int ngpt,
            int ngas, int nflav, int neta, int npres, int ntemp,
            int nminorlower, int nminorklower,
            int nminorupper, int nminorkupper,
            int idx_h2o,
            const Array<int,2>& gpoint_flavor,
            const Array<int,2>& band_lims_gpt,
            const Array<Float,4>& kmajor,
            const Array<Float,3>& kminor_lower,
            const Array<Float,3>& kminor_upper,
            const Array<int,2>& minor_limits_gpt_lower,
            const Array<int,2>& minor_limits_gpt_upper,
            const Array<Bool,1>& minor_scales_with_density_lower,
            const Array<Bool,1>& minor_scales_with_density_upper,
            const Array<Bool,1>& scale_by_complement_lower,
            const Array<Bool,1>& scale_by_complement_upper,
            const Array<int,1>& idx_minor_lower,
            const Array<int,1>& idx_minor_upper,
            const Array<int,1>& idx_minor_scaling_lower,
            const Array<int,1>& idx_minor_scaling_upper,
            const Array<int,1>& kminor_start_lower,
            const Array<int,1>& kminor_start_upper,
            const Array<Bool,2>& tropo,
            const Array<Float,4>& col_mix, const Array<Float,6>& fmajor,
            const Array<Float,5>& fminor, const Array<Float,2>& play,
            const Array<Float,2>& tlay, Array<Float,3>& col_gas,
            const Array<int,4>& jeta, const Array<int,2>& jtemp,
            const Array<int,2>& jpress, Array<Float,3>& tau)
    {
        rrtmgp_kernels::rrtmgp_compute_tau_absorption(
            &ncol, &nlay, &nband, &ngpt,
            &ngas, &nflav, &neta, &npres, &ntemp,
            &nminorlower, &nminorklower,
            &nminorupper, &nminorkupper,
            &idx_h2o,
            const_cast<int*>(gpoint_flavor.ptr()),
            const_cast<int*>(band_lims_gpt.ptr()),
            const_cast<Float*>(kmajor.ptr()),
            const_cast<Float*>(kminor_lower.ptr()),
            const_cast<Float*>(kminor_upper.ptr()),
            const_cast<int*>(minor_limits_gpt_lower.ptr()),
            const_cast<int*>(minor_limits_gpt_upper.ptr()),
            const_cast<Bool*>(minor_scales_with_density_lower.ptr()),
            const_cast<Bool*>(minor_scales_with_density_upper.ptr()),
            const_cast<Bool*>(scale_by_complement_lower.ptr()),
            const_cast<Bool*>(scale_by_complement_upper.ptr()),
            const_cast<int*>(idx_minor_lower.ptr()),
            const_cast<int*>(idx_minor_upper.ptr()),
            const_cast<int*>(idx_minor_scaling_lower.ptr()),
            const_cast<int*>(idx_minor_scaling_upper.ptr()),
            const_cast<int*>(kminor_start_lower.ptr()),
            const_cast<int*>(kminor_start_upper.ptr()),
            const_cast<Bool*>(tropo.ptr()),
            const_cast<Float*>(col_mix.ptr()), const_cast<Float*>(fmajor.ptr()), const_cast<Float*>(fminor.ptr()),
            const_cast<Float*>(play.ptr()), const_cast<Float*>(tlay.ptr()), const_cast<Float*>(col_gas.ptr()),
            const_cast<int*>(jeta.ptr()), const_cast<int*>(jtemp.ptr()), const_cast<int*>(jpress.ptr()),
            tau.ptr());
    }

    template<typename Float>
    void compute_tau_rayleigh(
            int ncol, int nlay, int nband, int ngpt,
            int ngas, int nflav, int neta, int npres, int ntemp,
            const Array<int,2>& gpoint_flavor,
            const Array<int,2>& band_lims_gpt,
            const Array<Float,4>& krayl,
            int idx_h2o, const Array<Float,2>& col_dry, const Array<Float,3>& col_gas,
            const Array<Float,5>& fminor, const Array<int,4>& jeta,
            const Array<Bool,2>& tropo, const Array<int,2>& jtemp,
            Array<Float,3>& tau_rayleigh)
    {
        rrtmgp_kernels::rrtmgp_compute_tau_rayleigh(
                &ncol, &nlay, &nband, &ngpt,
                &ngas, &nflav, &neta, &npres, &ntemp,
                const_cast<int*>(gpoint_flavor.ptr()),
                const_cast<int*>(band_lims_gpt.ptr()),
                const_cast<Float*>(krayl.ptr()),
                &idx_h2o,
                const_cast<Float*>(col_dry.ptr()), const_cast<Float*>(col_gas.ptr()),
                const_cast<Float*>(fminor.ptr()), const_cast<int*>(jeta.ptr()),
                const_cast<Bool*>(tropo.ptr()), const_cast<int*>(jtemp.ptr()),
                tau_rayleigh.ptr());
    }

    template<typename Float>
    void reorder123x321(
            const Array<Float,3>& data,
            Array<Float,3>& data_out)
    {
        int dim1 = data.dim(1);
        int dim2 = data.dim(2);
        int dim3 = data.dim(3);
        rrtmgp_kernels::reorder_123x321_kernel(
                &dim1, &dim2, &dim3,
                const_cast<Float*>(data.ptr()),
                data_out.ptr());
    }

    template<typename Float>
    void compute_Planck_source(
            int ncol, int nlay, int nbnd, int ngpt,
            int nflav, int neta, int npres, int ntemp, int nPlanckTemp,
            const Array<Float,2>& tlay, const Array<Float,2>& tlev, const Array<Float,1>& tsfc, int sfc_lay,
            const Array<Float,6>& fmajor, const Array<int,4>& jeta, const Array<Bool,2>& tropo, const Array<int,2>& jtemp, const Array<int,2>& jpress,
            const Array<int,1>& gpoint_bands, const Array<int,2>& band_lims_gpt, const Array<Float,4>& pfracin, Float temp_ref_min,
            Float totplnk_delta, const Array<Float,2>& totplnk, const Array<int,2>& gpoint_flavor,
            Array<Float,2>& sfc_src, Array<Float,3>& lay_src, Array<Float,3>& lev_src_inc, Array<Float,3>& lev_src_dec,
            Array<Float,2>& sfc_src_jac)
    {
        rrtmgp_kernels::rrtmgp_compute_Planck_source(
                &ncol, &nlay, &nbnd, &ngpt,
                &nflav, &neta, &npres, &ntemp, &nPlanckTemp,
                const_cast<Float*>(tlay.ptr()),
                const_cast<Float*>(tlev.ptr()),
                const_cast<Float*>(tsfc.ptr()),
                &sfc_lay,
                const_cast<Float*>(fmajor.ptr()),
                const_cast<int*>(jeta.ptr()),
                const_cast<Bool*>(tropo.ptr()),
                const_cast<int*>(jtemp.ptr()),
                const_cast<int*>(jpress.ptr()),
                const_cast<int*>(gpoint_bands.ptr()), const_cast<int*>(band_lims_gpt.ptr()), const_cast<Float*>(pfracin.ptr()), &temp_ref_min,
                &totplnk_delta, const_cast<Float*>(totplnk.ptr()), const_cast<int*>(gpoint_flavor.ptr()),
                sfc_src.ptr(), lay_src.ptr(), lev_src_inc.ptr(), lev_src_dec.ptr(),
                sfc_src_jac.ptr());
    }
}


void Gas_optics_rrtmgp::compute_gas_taus(
        const int ncol, const int nlay, const int ngpt, const int nband,
        const Array<Float,2>& play,
        const Array<Float,2>& plev,
        const Array<Float,2>& tlay,
        const Gas_concs& gas_desc,
        std::unique_ptr<Optical_props_arry>& optical_props,
        Array<int,2>& jtemp, Array<int,2>& jpress,
        Array<int,4>& jeta,
        Array<Bool,2>& tropo,
        Array<Float,6>& fmajor,
        const Array<Float,2>& col_dry) const
{
    Array<Float,3> vmr({ncol, nlay, this->get_ngas()});
    Array<Float,3> col_gas({ncol, nlay, this->get_ngas()+1});
    col_gas.set_offsets({0, 0, -1});
    Array<Float,4> col_mix({2, ncol, nlay, this->get_nflav()});
    Array<Float,5> fminor({2, 2, ncol, nlay, this->get_nflav()});

    // CvH add all the checking...
    const int ngas = this->get_ngas();
    const int nflav = this->get_nflav();
    const int neta = this->get_neta();
    const int npres = this->get_npres();
    const int ntemp = this->get_ntemp();

    const int nminorlower = this->minor_scales_with_density_lower.dim(1);
    const int nminorklower = this->kminor_lower.dim(3);
    const int nminorupper = this->minor_scales_with_density_upper.dim(1);
    const int nminorkupper = this->kminor_upper.dim(3);

    for (int igas=1; igas<=ngas; ++igas)
    {
        const Array<Float,2>& vmr_2d = gas_desc.get_vmr(this->gas_names({igas}));

        // Fill array with constant value.
        if (vmr_2d.dim(1) == 1 && vmr_2d.dim(2) == 1)
        {
            const Float vmr_c = vmr_2d({1, 1});
            for (int ilay=1; ilay<=nlay; ++ilay)
                for (int icol=1; icol<=ncol; ++icol)
                    vmr({icol, ilay, igas}) = vmr_c;
        }
        // Fill array with constant profile.
        else if (vmr_2d.dim(1) == 1)
        {
            for (int ilay=1; ilay<=nlay; ++ilay)
            {
                const Float vmr_lay = vmr_2d({1, ilay});
                for (int icol=1; icol<=ncol; ++icol)
                    vmr({icol, ilay, igas}) = vmr_lay;
            }
        }
        // Fill array with full 2d data.
        else
        {
            for (int ilay=1; ilay<=nlay; ++ilay)
                for (int icol=1; icol<=ncol; ++icol)
                    vmr({icol, ilay, igas}) = vmr_2d({icol, ilay});
        }
    }

    // CvH: Assume that col_dry is provided.
    for (int ilay=1; ilay<=nlay; ++ilay)
        for (int icol=1; icol<=ncol; ++icol)
            col_gas({icol, ilay, 0}) = col_dry({icol, ilay});

    for (int igas=1; igas<=ngas; ++igas)
        for (int ilay=1; ilay<=nlay; ++ilay)
            for (int icol=1; icol<=ncol; ++icol)
                col_gas({icol, ilay, igas}) = vmr({icol, ilay, igas}) * col_dry({icol, ilay});

    // Call the fortran kernels
    rrtmgp_kernel_launcher::interpolation(
            ncol, nlay,
            ngas, nflav, neta, npres, ntemp,
            this->flavor,
            this->press_ref_log,
            this->temp_ref,
            this->press_ref_log_delta,
            this->temp_ref_min,
            this->temp_ref_delta,
            this->press_ref_trop_log,
            this->vmr_ref,
            play,
            tlay,
            col_gas,
            jtemp,
            fmajor, fminor,
            col_mix,
            tropo,
            jeta, jpress);

    int idx_h2o = -1;
    for (int i=1; i<=this->gas_names.dim(1); ++i)
        if (gas_names({i}) == "h2o")
        {
            idx_h2o = i;
            break;
        }

    if (idx_h2o == -1)
        throw std::runtime_error("idx_h2o cannot be found");

    bool has_rayleigh = (this->krayl.size() > 0);

    if (has_rayleigh)
    {
        Array<Float,3> tau({ncol, nlay, ngpt});
        Array<Float,3> tau_rayleigh({ncol, nlay, ngpt});

        rrtmgp_kernel_launcher::compute_tau_absorption(
                ncol, nlay, nband, ngpt,
                ngas, nflav, neta, npres, ntemp,
                nminorlower, nminorklower,
                nminorupper, nminorkupper,
                idx_h2o,
                this->gpoint_flavor,
                this->get_band_lims_gpoint(),
                this->kmajor,
                this->kminor_lower,
                this->kminor_upper,
                this->minor_limits_gpt_lower,
                this->minor_limits_gpt_upper,
                this->minor_scales_with_density_lower,
                this->minor_scales_with_density_upper,
                this->scale_by_complement_lower,
                this->scale_by_complement_upper,
                this->idx_minor_lower,
                this->idx_minor_upper,
                this->idx_minor_scaling_lower,
                this->idx_minor_scaling_upper,
                this->kminor_start_lower,
                this->kminor_start_upper,
                tropo,
                col_mix, fmajor, fminor,
                play, tlay, col_gas,
                jeta, jtemp, jpress,
                tau);

        rrtmgp_kernel_launcher::compute_tau_rayleigh(
                ncol, nlay, nband, ngpt,
                ngas, nflav, neta, npres, ntemp,
                this->gpoint_flavor,
                this->get_band_lims_gpoint(),
                this->krayl,
                idx_h2o, col_dry, col_gas,
                fminor, jeta, tropo, jtemp,
                tau_rayleigh);

        combine_abs_and_rayleigh(tau, tau_rayleigh, optical_props);
    }
    else
    {
        rrtmgp_kernel_launcher::zero_array(ncol, nlay, ngpt, optical_props->get_tau());

        rrtmgp_kernel_launcher::compute_tau_absorption(
                ncol, nlay, nband, ngpt,
                ngas, nflav, neta, npres, ntemp,
                nminorlower, nminorklower,
                nminorupper, nminorkupper,
                idx_h2o,
                this->gpoint_flavor,
                this->get_band_lims_gpoint(),
                this->kmajor,
                this->kminor_lower,
                this->kminor_upper,
                this->minor_limits_gpt_lower,
                this->minor_limits_gpt_upper,
                this->minor_scales_with_density_lower,
                this->minor_scales_with_density_upper,
                this->scale_by_complement_lower,
                this->scale_by_complement_upper,
                this->idx_minor_lower,
                this->idx_minor_upper,
                this->idx_minor_scaling_lower,
                this->idx_minor_scaling_upper,
                this->kminor_start_lower,
                this->kminor_start_upper,
                tropo,
                col_mix, fmajor, fminor,
                play, tlay, col_gas,
                jeta, jtemp, jpress,
                optical_props->get_tau());
    }
}


void Gas_optics_rrtmgp::combine_abs_and_rayleigh(
        const Array<Float,3>& tau,
        const Array<Float,3>& tau_rayleigh,
        std::unique_ptr<Optical_props_arry>& optical_props) const
{
    int ncol = tau.dim(1);
    int nlay = tau.dim(2);
    int ngpt = tau.dim(3);

    /*
    for (int igpt=1; igpt<=ngpt; ++igpt)
        for (int ilay=1; ilay<=nlay; ++ilay)
            for (int icol=1; icol<=ncol; ++icol)
            {
                const Float t = tau({icol, ilay, igpt}) + tau_rayleigh({icol, ilay, igpt});

                if (t > Float(2.) * std::numeric_limits<Float>::epsilon())
                    optical_props->get_ssa()({icol, ilay, igpt}) = tau_rayleigh({icol, ilay, igpt}) / t;
                else
                    optical_props->get_ssa()({icol, ilay, igpt}) = Float(0.);

                optical_props->get_tau()({icol, ilay, igpt}) = t;
            }
            */

    combine_abs_and_rayleigh_kernel(
            ncol, nlay, ngpt,
            tau.ptr(), tau_rayleigh.ptr(),
            optical_props->get_tau().ptr(), optical_props->get_ssa().ptr());

    rrtmgp_kernel_launcher::zero_array(ncol, nlay, ngpt, optical_props->get_g());
}


void Gas_optics_rrtmgp::source(
        const int ncol, const int nlay, const int nbnd, const int ngpt,
        const Array<Float,2>& play, const Array<Float,2>& plev,
        const Array<Float,2>& tlay, const Array<Float,1>& tsfc,
        const Array<int,2>& jtemp, const Array<int,2>& jpress,
        const Array<int,4>& jeta, const Array<Bool,2>& tropo,
        const Array<Float,6>& fmajor,
        Source_func_lw& sources,
        const Array<Float,2>& tlev) const
{
    // CvH Assume tlev is available.
    // Compute internal (Planck) source functions at layers and levels,
    // which depend on mapping from spectral space that creates k-distribution.
    const int nflav = this->get_nflav();
    const int neta = this->get_neta();
    const int npres = this->get_npres();
    const int ntemp = this->get_ntemp();
    const int nPlanckTemp = this->get_nPlanckTemp();
    auto gpoint_bands = this->get_gpoint_bands();
    auto band_lims_gpoint = this->get_band_lims_gpoint();

    int sfc_lay = play({1, 1}) > play({1, nlay}) ? 1 : nlay;

    rrtmgp_kernel_launcher::compute_Planck_source(
            ncol, nlay, nbnd, ngpt,
            nflav, neta, npres, ntemp, nPlanckTemp,
            tlay, tlev, tsfc, sfc_lay,
            fmajor, jeta, tropo, jtemp, jpress,
            gpoint_bands, band_lims_gpoint, this->planck_frac, this->temp_ref_min,
            this->totplnk_delta, this->totplnk, this->gpoint_flavor,
            sources.get_sfc_source(), sources.get_lay_source(), sources.get_lev_source_inc(), sources.get_lev_source_dec(),
            sources.get_sfc_source_jac());
}


Float Gas_optics_rrtmgp::get_tsi() const
{
    const int n_gpt = this->get_ngpt();

    Float tsi = 0.;
    for (int igpt=1; igpt<=n_gpt; ++igpt)
        tsi += this->solar_source({igpt});

    return tsi;
}
