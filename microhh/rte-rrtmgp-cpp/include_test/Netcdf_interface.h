/*
 * This file is imported from MicroHH (https://github.com/earth-system-radiation/earth-system-radiation)
 * and is adapted for the testing of the C++ interface to the
 * RTE+RRTMGP radiation code.
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef NETCDF_INTERFACE_H
#define NETCDF_INTERFACE_H

#include <vector>
#include <map>
#include <iostream>
#include <numeric>
#include <netcdf.h>

#include "Status.h"

enum class Netcdf_mode { Create, Read, Write };

class Netcdf_handle;
class Netcdf_group;

template<typename T>
class Netcdf_variable
{
    public:
        Netcdf_variable(Netcdf_handle&, const int, const std::vector<int>&);
        Netcdf_variable(const Netcdf_variable&) = default;
        void insert(const std::vector<T>&, const std::vector<int>);
        void insert(const std::vector<T>&, const std::vector<int>, const std::vector<int>);
        void insert(const T, const std::vector<int>);
        const std::vector<int> get_dim_sizes() { return dim_sizes; }
        void add_attribute(const std::string&, const std::string&);
        void add_attribute(const std::string&, const double);
        void add_attribute(const std::string&, const float);

    private:
        Netcdf_handle& nc_file;
        const int var_id;
        const std::vector<int> dim_sizes;
};

class Netcdf_handle
{
    public:
        Netcdf_handle();
        void add_dimension(const std::string&, const int dim_size = NC_UNLIMITED);

        Netcdf_group add_group(const std::string&);
        Netcdf_group get_group(const std::string&) const;

        int get_dimension_size(const std::string&);

        std::map<std::string, int> get_variable_dimensions(const std::string&) const;

        bool variable_exists(const std::string&) const;

        template<typename T>
        Netcdf_variable<T> add_variable(
                const std::string&,
                const std::vector<std::string>&);

        template<typename T>
        Netcdf_variable<T> add_variable(
                const std::string&);

        template<typename T>
        T get_variable(
            const std::string&) const;

        template<typename T>
        std::vector<T> get_variable(
            const std::string&,
            const std::vector<int>&) const;

        template<typename T>
        void get_variable(
                std::vector<T>&,
                const std::string&,
                const std::vector<int>&,
                const std::vector<int>&) const;

        template<typename T>
        void insert(
                const std::vector<T>&,
                const int var_id,
                const std::vector<int>&,
                const std::vector<int>&);

        template<typename T>
        void insert(
                const T,
                const int var_id,
                const std::vector<int>&,
                const std::vector<int>&);

        void add_attribute(
                const std::string&,
                const std::string&,
                const int);

        void add_attribute(
                const std::string&,
                const double,
                const int);

        void add_attribute(
                const std::string&,
                const float,
                const int);

    protected:
        int ncid;
        int root_ncid;
        std::map<std::string, int> dims;
        int record_counter;
};

class Netcdf_file : public Netcdf_handle
{
    public:
        Netcdf_file(const std::string&, Netcdf_mode);
        ~Netcdf_file();

        void sync();
};

class Netcdf_group : public Netcdf_handle
{
    public:
        Netcdf_group(const int, const int);
};

// Implementation.
namespace
{
    void nc_throw(const int return_value)
    {
        std::string error(nc_strerror(return_value));
        throw std::runtime_error(error);
    }

    void nc_check(int return_value)
    {
        if (return_value != NC_NOERR)
            nc_throw(return_value);
    }

    // Get the NetCDF data type based on TF
    template<typename TF> nc_type netcdf_dtype();
    template<> nc_type netcdf_dtype<double>() { return NC_DOUBLE; }
    template<> nc_type netcdf_dtype<float>()  { return NC_FLOAT; }
    template<> nc_type netcdf_dtype<int>()    { return NC_INT; }

    // Wrapper for the `nc_get_vara_TYPE` functions
    template<typename TF>
    int nc_get_vara_wrapper(
            int, int, const std::vector<size_t>&, const std::vector<size_t>&, std::vector<TF>&);

    template<>
    int nc_get_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, std::vector<double>& values)
    {
        return nc_get_vara_double(ncid, var_id, start.data(), count.data(), values.data());
    }

    template<>
    int nc_get_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, std::vector<float>& values)
    {
        return nc_get_vara_float(ncid, var_id, start.data(), count.data(), values.data());
    }

    template<>
    int nc_get_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, std::vector<int>& values)
    {
        return nc_get_vara_int(ncid, var_id, start.data(), count.data(), values.data());
    }

    template<>
    int nc_get_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, std::vector<char>& values)
    {
        return nc_get_vara_text(ncid, var_id, start.data(), count.data(), values.data());
    }

    template<>
    int nc_get_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, std::vector<signed char>& values)
    {
        return nc_get_vara_schar(ncid, var_id, start.data(), count.data(), values.data());
    }

    // Wrapper for the `nc_put_vara_TYPE` functions
    template<typename TF>
    int nc_put_vara_wrapper(
            int, int, const std::vector<size_t>&, const std::vector<size_t>&, const std::vector<TF>&);

    template<>
    int nc_put_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, const std::vector<double>& values)
    {
        return nc_put_vara_double(ncid, var_id, start.data(), count.data(), values.data());
    }

    template<>
    int nc_put_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, const std::vector<float>& values)
    {
        return nc_put_vara_float(ncid, var_id, start.data(), count.data(), values.data());
    }

    template<>
    int nc_put_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, const std::vector<int>& values)
    {
        return nc_put_vara_int(ncid, var_id, start.data(), count.data(), values.data());
    }


    template<typename TF>
    int nc_put_vara_wrapper(
            int, int, const std::vector<size_t>&, const std::vector<size_t>&, const TF);

    template<>
    int nc_put_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, const double value)
    {
        return nc_put_vara_double(ncid, var_id, start.data(), count.data(), &value);
    }

    template<>
    int nc_put_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, const float value)
    {
        return nc_put_vara_float(ncid, var_id, start.data(), count.data(), &value);
    }

    template<>
    int nc_put_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, const int value)
    {
        return nc_put_vara_int(ncid, var_id, start.data(), count.data(), &value);
    }
}

inline Netcdf_file::Netcdf_file(const std::string& name, Netcdf_mode mode) :
        Netcdf_handle()
{
    int nc_check_code = 0;

    if (mode == Netcdf_mode::Create)
        nc_check_code = nc_create(name.c_str(), NC_NOCLOBBER | NC_NETCDF4, &ncid);
    else if (mode == Netcdf_mode::Write)
        nc_check_code = nc_open(name.c_str(), NC_WRITE | NC_NETCDF4, &ncid);
    else if (mode == Netcdf_mode::Read)
    {
        // CvH: I changed this line, as I can no longer open the "classic" coefficient files.
        // nc_check_code = nc_open(name.c_str(), NC_NOWRITE | NC_NETCDF4, &ncid);
        nc_check_code = nc_open(name.c_str(), NC_NOWRITE, &ncid);
    }

    try
    {
        nc_check(nc_check_code);
    }
    catch (std::runtime_error& e)
    {
        std::string error = "Opening of file " + name + " returned: " + e.what();
        throw std::runtime_error(error);
    }

    root_ncid = ncid;

    if (mode == Netcdf_mode::Create)
        nc_check_code =  nc_enddef(root_ncid);

    nc_check(nc_check_code);
}

inline Netcdf_file::~Netcdf_file()
{
    int nc_check_code = 0;

    nc_check_code = nc_close(ncid);
    nc_check(nc_check_code);
}

inline void Netcdf_file::sync()
{
    int nc_check_code = 0;

    nc_check_code = nc_sync(ncid);
    nc_check(nc_check_code);
}

inline void Netcdf_handle::add_dimension(const std::string& dim_name, const int dim_size)
{
    int nc_check_code = 0;

    nc_check_code = nc_redef(root_ncid);
    nc_check(nc_check_code);

    int dim_id;
    int def_out;

    def_out = nc_def_dim(ncid, dim_name.c_str(), dim_size, &dim_id);

    // Dimension is written or already exists.
    if (def_out == NC_NOERR)
        dims.emplace(dim_name, dim_id);
    else if (def_out == NC_ENAMEINUSE)
    {}
        // Error.
    else
        nc_throw(def_out);

    nc_check_code = nc_enddef(root_ncid);
    nc_check(nc_check_code);
}

template<typename T>
inline Netcdf_variable<T> Netcdf_handle::add_variable(
        const std::string& var_name)
{
    int nc_check_code = 0;

    int var_id = -1;

    nc_check_code = nc_redef(root_ncid);
    nc_check(nc_check_code);

    std::vector<int> dims = {1};

    nc_check_code = nc_def_var(ncid, var_name.c_str(), netcdf_dtype<T>(), 0, dims.data(), &var_id);
    nc_check(nc_check_code);

    nc_check_code = nc_enddef(root_ncid);
    nc_check(nc_check_code);

    return Netcdf_variable<T>(*this, var_id, dims);
}

template<typename T>
inline Netcdf_variable<T> Netcdf_handle::add_variable(
        const std::string& var_name,
        const std::vector<std::string>& dim_names)
{
    int nc_check_code = 0;

    int var_id = -1;
    std::vector<int> dim_sizes;

    nc_check_code = nc_redef(root_ncid);
    nc_check(nc_check_code);

    int ndims = dim_names.size();
    std::vector<int> dim_ids;

    for (const std::string& dim_name : dim_names)
        dim_ids.push_back(dims.at(dim_name));

    nc_check_code = nc_def_var(ncid, var_name.c_str(), netcdf_dtype<T>(), ndims, dim_ids.data(), &var_id);
    nc_check(nc_check_code);

    nc_check_code = nc_enddef(root_ncid);
    nc_check(nc_check_code);

    // Broadcast the dim_ids size of the main process to run the for loop on all processes.
    int dim_ids_size = dim_ids.size();

    for (int i=0; i<dim_ids_size; ++i)
    {
        const int dim_id = dim_ids.at(i);

        size_t dim_len = 0;

        nc_check_code = nc_inq_dimlen(ncid, dim_id, &dim_len);
        nc_check(nc_check_code);

        int dim_len_int = static_cast<int>(dim_len);

        if (dim_len_int == NC_UNLIMITED)
        {
            if (i == 0)
                dim_sizes.push_back(1);
            else
                throw std::runtime_error("Only the leftmost dimension is allowed to be an NC_UNLIMITED dimension");
        }
        else
            dim_sizes.push_back(dim_len);
    }

    return Netcdf_variable<T>(*this, var_id, dim_sizes);
}

inline Netcdf_handle::Netcdf_handle() :
        record_counter(0)
{}

template<typename T>
inline void Netcdf_handle::insert(
        const std::vector<T>& values,
        const int var_id,
        const std::vector<int>& i_start,
        const std::vector<int>& i_count)
{
    const std::vector<size_t> i_start_size_t (i_start.begin(), i_start.end());
    const std::vector<size_t> i_count_size_t (i_count.begin(), i_count.end());

    int nc_check_code = 0;

    // CvH: Add proper size checking.
    nc_check_code = nc_put_vara_wrapper<T>(ncid, var_id, i_start_size_t, i_count_size_t, values);
    nc_check(nc_check_code);
}

template<typename T>
inline void Netcdf_handle::insert(
        const T value,
        const int var_id,
        const std::vector<int>& i_start,
        const std::vector<int>& i_count)
{
    const std::vector<size_t> i_start_size_t (i_start.begin(), i_start.end());
    const std::vector<size_t> i_count_size_t (i_count.begin(), i_count.end());

    int nc_check_code = 0;

    // CvH: Add proper size checking.
    nc_check_code = nc_put_vara_wrapper<T>(ncid, var_id, i_start_size_t, i_count_size_t, value);
    nc_check(nc_check_code);
}

inline void Netcdf_handle::add_attribute(
        const std::string& name,
        const std::string& value,
        const int var_id)
{
    int nc_check_code = 0;

    nc_check_code = nc_redef(root_ncid);
    nc_check(nc_check_code);

    // CvH what if string is too long?
    nc_check_code = nc_put_att_text(ncid, var_id, name.c_str(), value.size(), value.c_str());
    nc_check(nc_check_code);

    nc_check_code = nc_enddef(root_ncid);
    nc_check(nc_check_code);
}

inline void Netcdf_handle::add_attribute(
        const std::string& name,
        const double value,
        const int var_id)
{
    int nc_check_code = 0;

    nc_check_code = nc_redef(root_ncid);
    nc_check(nc_check_code);

    // CvH what if string is too long?
    nc_check_code = nc_put_att_double(ncid, var_id, name.c_str(), NC_DOUBLE, 1, &value);
    nc_check(nc_check_code);

    nc_check_code = nc_enddef(root_ncid);
    nc_check(nc_check_code);
}

inline void Netcdf_handle::add_attribute(
        const std::string& name,
        const float value,
        const int var_id)
{
    int nc_check_code = 0;

    nc_check_code = nc_redef(root_ncid);
    nc_check(nc_check_code);

    // CvH what if string is too long?
    nc_check_code = nc_put_att_float(ncid, var_id, name.c_str(), NC_FLOAT, 1, &value);
    nc_check(nc_check_code);

    nc_check_code = nc_enddef(root_ncid);
    nc_check(nc_check_code);
}

inline Netcdf_group Netcdf_handle::add_group(const std::string& name)
{
    int group_ncid = -1;
    int nc_check_code = 0;

    nc_check_code = nc_redef(root_ncid);
    nc_check(nc_check_code);

    nc_check_code = nc_def_grp(ncid, name.c_str(), &group_ncid);
    nc_check(nc_check_code);

    nc_check_code = nc_enddef(root_ncid);
    nc_check(nc_check_code);

    return Netcdf_group(group_ncid, root_ncid);
}

inline Netcdf_group Netcdf_handle::get_group(const std::string& name) const
{
    int group_ncid = -1;
    int nc_check_code = 0;

    nc_check_code = nc_inq_ncid(ncid, name.c_str(), &group_ncid);
    nc_check(nc_check_code);

    return Netcdf_group(group_ncid, root_ncid);
}

inline int Netcdf_handle::get_dimension_size(const std::string& name)
{
    int nc_check_code = 0;
    int dim_id;

    nc_check_code = nc_inq_dimid(ncid, name.c_str(), &dim_id);
    nc_check(nc_check_code);

    size_t dim_len_size_t = 0;
    nc_check_code = nc_inq_dimlen(ncid, dim_id, &dim_len_size_t);
    nc_check(nc_check_code);

    int dim_len = dim_len_size_t;

    return dim_len;
}

inline Netcdf_group::Netcdf_group(const int ncid_in, const int root_ncid_in) :
        Netcdf_handle()
{
    ncid = ncid_in;
    root_ncid = root_ncid_in;
}

inline std::map<std::string, int> Netcdf_handle::get_variable_dimensions(
        const std::string& name) const
{
    int nc_check_code = 0;
    int var_id;

    nc_check_code = nc_inq_varid(ncid, name.c_str(), &var_id);
    nc_check(nc_check_code);

    int ndims;
    int dimids[NC_MAX_VAR_DIMS];

    nc_check_code = nc_inq_var(ncid, var_id, NULL, NULL, &ndims, dimids, NULL);
    nc_check(nc_check_code);

    std::map<std::string, int> dims;

    for (int n=0; n<ndims; ++n)
    {
        char dim_name[NC_MAX_NAME+1];
        size_t dim_length_size_t;

        nc_check_code = nc_inq_dim(ncid, dimids[n], dim_name, &dim_length_size_t);
        nc_check(nc_check_code);

        int dim_length = dim_length_size_t;

        dims.emplace(std::string(dim_name), dim_length);
    }

    return dims;
}

inline bool Netcdf_handle::variable_exists(const std::string& name) const
{
    int nc_check_code = 0;
    int var_id;

    try
    {
        nc_check_code = nc_inq_varid(ncid, name.c_str(), &var_id);
        nc_check(nc_check_code);
    }
    catch (std::runtime_error& e)
    {
        return false;
    }

    return true;
}

template<typename TF>
inline TF Netcdf_handle::get_variable(
        const std::string& name) const
{
    // std::string message = "Retrieving from NetCDF (single value): " + name;
    // Status::print_message(message);

    int nc_check_code = 0;
    int var_id;

    try
    {
        nc_check_code = nc_inq_varid(ncid, name.c_str(), &var_id);
        nc_check(nc_check_code);
    }
    catch (std::runtime_error& e)
    {
        std::string error = "Netcdf variable: " + name + " not found";
        Status::print_error(error);
        throw;
    }

    TF value = 0;
    std::vector<TF> values(1);

    nc_check_code = nc_get_vara_wrapper(ncid, var_id, {0}, {1}, values);
    nc_check(nc_check_code);

    value = values[0];

    return value;
}

template<typename TF>
inline std::vector<TF> Netcdf_handle::get_variable(
        const std::string& name,
        const std::vector<int>& i_count) const
{
    // std::string message = "Retrieving from NetCDF (full array): " + name;
    // Status::print_message(message);

    const std::vector<size_t> i_start_size_t(i_count.size());
    const std::vector<size_t> i_count_size_t(i_count.begin(), i_count.end());

    int nc_check_code = 0;
    int var_id;

    try
    {
        nc_check_code = nc_inq_varid(ncid, name.c_str(), &var_id);
        nc_check(nc_check_code);
    }
    catch (std::runtime_error& e)
    {
        std::string error = "Netcdf variable: " + name + " not found";
        Status::print_error(error);
        throw;
    }

    int total_count = std::accumulate(i_count.begin(), i_count.end(), 1, std::multiplies<>());
    // CvH check needs to be added if total count matches multiplication of all dimensions.

    std::vector<TF> values(total_count);
    nc_check_code = nc_get_vara_wrapper(ncid, var_id, i_start_size_t, i_count_size_t, values);
    nc_check(nc_check_code);

    return values;
}

template<typename TF>
inline void Netcdf_handle::get_variable(
        std::vector<TF>& values,
        const std::string& name,
        const std::vector<int>& i_start,
        const std::vector<int>& i_count) const
{
    // std::string message = "Retrieving from NetCDF: " + name;
    // Status::print_message(message);

    const std::vector<size_t> i_start_size_t (i_start.begin(), i_start.end());
    const std::vector<size_t> i_count_size_t (i_count.begin(), i_count.end());

    int nc_check_code = 0;
    int var_id;

    bool zero_fill = false;
    try
    {
        nc_check_code = nc_inq_varid(ncid, name.c_str(), &var_id);
        nc_check(nc_check_code);
    }
    catch (std::runtime_error& e)
    {
        std::string warning = "Netcdf variable: " + name + " not found, filling with zeros";
        Status::print_warning(warning);
        zero_fill = true;
    }

    // CvH: Add check if the vector is large enough.
    int total_count = std::accumulate(i_count.begin(), i_count.end(), 1, std::multiplies<>());

    if (zero_fill)
    {
        std::fill(values.begin(), values.begin() + total_count, 0);
    }
    else
    {
        nc_check_code = nc_get_vara_wrapper(ncid, var_id, i_start_size_t, i_count_size_t, values);
        nc_check(nc_check_code);
    }
}

// Variable does not communicate with NetCDF library directly.
template<typename T>
inline Netcdf_variable<T>::Netcdf_variable(Netcdf_handle& nc_file, const int var_id, const std::vector<int>& dim_sizes) :
        nc_file(nc_file), var_id(var_id), dim_sizes(dim_sizes)
{}

template<typename T>
inline void Netcdf_variable<T>::insert(const std::vector<T>& values, const std::vector<int> i_start)
{
    nc_file.insert(values, var_id, i_start, dim_sizes);
}

template<typename T>
inline void Netcdf_variable<T>::insert(
        const std::vector<T>& values,
        const std::vector<int> i_start,
        const std::vector<int> i_count)
{
    nc_file.insert(values, var_id, i_start, i_count);
}

template<typename T>
inline void Netcdf_variable<T>::insert(const T value, const std::vector<int> i_start)
{
    nc_file.insert(value, var_id, i_start, dim_sizes);
}

template<typename T>
inline void Netcdf_variable<T>::add_attribute(const std::string& name, const std::string& value)
{
    nc_file.add_attribute(name, value, var_id);
}

template<typename T>
inline void Netcdf_variable<T>::add_attribute(const std::string& name, const double value)
{
    nc_file.add_attribute(name, value, var_id);
}

template<typename T>
inline void Netcdf_variable<T>::add_attribute(const std::string& name, const float value)
{
    nc_file.add_attribute(name, value, var_id);
}
#endif
