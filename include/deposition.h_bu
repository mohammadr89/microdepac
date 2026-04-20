#ifndef DEPOSITION_H
#define DEPOSITION_H

#include <vector>
#include <string>
#include <map>
#include "timeloop.h"
#include "boundary_surface_lsm.h"
#include "boundary.h"
#include "boundary_cyclic.h"
#include "radiation.h" 

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;
template<typename> class Cross;
template<typename> class Boundary_surface_lsm;
template<typename> class Radiation;
template<typename> class Chemistry;
enum class Deposition_type {disabled, enabled, simple, average};
template<typename TF>
struct Deposition_tile
{
    std::string long_name;
    std::vector<TF> vdnh3;
    std::vector<TF> ra;
    std::vector<TF> rb;
    std::vector<TF> T_surface;
    std::vector<TF> rh_surface;
    std::vector<TF> obuk;
    std::vector<TF> ustar;
    std::vector<TF> ccomp_tot;
    std::vector<TF> cw;
    std::vector<TF> cstom;
    std::vector<TF> csoil_eff;
    std::vector<TF> cw_out;
    std::vector<TF> cstom_out;
    std::vector<TF> csoil_out;
    std::vector<TF> rc_tot;
    std::vector<TF> rc_eff;
};

template<typename TF>
using Deposition_tile_map = std::map<std::string, Deposition_tile<TF>>;

template<typename TF>
class Deposition
{
    public:

    Deposition(Master&, Grid<TF>&, Fields<TF>&, Radiation<TF>&, Chemistry<TF>&, Input&);
    ~Deposition();
    void init(Input&);
    void create(Stats<TF>&, Cross<TF>&);
    void update_time_dependent
        (
        Timeloop<TF>&, 
        Boundary<TF>&,
        Thermo<TF>& thermo,
        TF* restrict vdnh3
        );
    const TF get_vd(const std::string&) const;
    void get_tiled_mean(TF*, std::string, TF, const TF*, const TF*, const TF*);
    void update_vd_water(TF*, std::string, const TF*, const TF*, const int*, const TF*, const TF*);
    void exec_cross(Cross<TF>&, unsigned long);
    void spatial_avg_vd(TF*);

    protected:

    std::vector<std::string> cross_list;

    private:

    Master& master;
    Grid<TF>& grid;
    Fields<TF>& fields;
    Radiation<TF>& radiation; 
    Chemistry<TF>& chemistry; 
    bool sw_deposition;
    bool use_depac;
    void sync_depac_fields(Boundary<TF>& boundary);
    std::shared_ptr<Boundary_surface_lsm<TF>> boundary_surface_lsm;
    TF deposition_var;
    TF henry_so2;
    TF rsoil_so2;
    TF rwat_so2;
    TF rws_so2;
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
    TF t0;
    TF td;
    TF max_rad;
    TF glrad;
    TF sinphi;
    TF temperature;
    TF rh;
    TF sai;
    TF lat;
    bool sw_sinphi_prescr;
    TF t_sunrise;
    TF t_sunset;
    int day_of_year;
    int nwet;
    int nwet_veg;
    int nwet_soil;
    int nwet_wet;
    int lu;
    int iratns;
    TF hlaw;
    TF react;
    TF c_ave_prev_nh3;
    TF catm;
    TF c_ug;
    TF pressure;
    std::vector<int> lu_map;
    std::vector<std::string> deposition_tile_names {"veg", "soil", "wet"};
    Deposition_tile_map<TF> deposition_tiles;
    const std::string tend_name = "deposition";
    const std::string tend_longname = "Deposition";
    bool sw_override_ccomp;
    TF ccomp_override_value;
    std::vector<TF> ra_mean;
    std::vector<TF> rb_mean;
    std::vector<TF> obuk_mean;
    std::vector<TF> ustar_mean;
    std::vector<TF> ccomp_mean;
    std::vector<TF> cw_mean;
    std::vector<TF> cstom_mean;
    std::vector<TF> csoil_eff_mean;
    std::vector<TF> cw_out_mean;
    std::vector<TF> cstom_out_mean;
    std::vector<TF> csoil_out_mean;
    std::vector<TF> rc_tot_mean;
    std::vector<TF> rc_eff_mean;
    std::vector<TF> T_surface_mean;
    std::vector<TF> rh_surface_mean;
};
#endif
