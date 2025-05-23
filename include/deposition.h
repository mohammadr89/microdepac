#ifndef DEPOSITION_H
#define DEPOSITION_H

#include <vector>
#include <string>
#include <map>
#include "timeloop.h"
#include "boundary_surface_lsm.h"
#include "boundary.h"
#include "boundary_cyclic.h"
#include "radiation.h" // Include radiation header

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Stats;
template<typename> class Cross;
template<typename> class Boundary_surface_lsm;
template<typename> class Radiation;

enum class Deposition_type {disabled, enabled, simple, average};

template<typename TF>
struct Deposition_tile
{
    std::string long_name;    // Descriptive name of tile
                              // Land surface
    std::vector<TF> vdnh3;    // Deposition velocity of NH3 (m s-1)
    std::vector<TF> ra;
    std::vector<TF> rb;
    std::vector<TF> obuk;
    std::vector<TF> ustar;
    std::vector<TF> ccomp_tot; // Compensation point (ug/m3)
    std::vector<TF> cw;       // External leaf resistance (s/m)
    std::vector<TF> cstom;    // Stomatal resistance (s/m)
    std::vector<TF> csoil_eff; // Soil effective resistance (s/m)

    // New fields for compensation points
    std::vector<TF> cw_out;   // External leaf compensation point (ug/m3)
    std::vector<TF> cstom_out; // Stomatal compensation point (ug/m3)
    std::vector<TF> csoil_out; // Soil compensation point (ug/m3)
    std::vector<TF> rc_tot;   // Total canopy resistance (s/m)
    std::vector<TF> rc_eff;   // Effective canopy resistance (s/m)
};

template<typename TF>
using Deposition_tile_map = std::map<std::string, Deposition_tile<TF>>;

template<typename TF>
class Deposition
{
public:
    Deposition(Master&, Grid<TF>&, Fields<TF>&, Radiation<TF>&, Input&);
    ~Deposition();

    void init(Input&);
    void create(Stats<TF>&, Cross<TF>&);
    void update_time_dependent(
                               Timeloop<TF>&, 
                               Boundary<TF>&,
                               Thermo<TF>& thermo,
                               TF* restrict vdnh3);  // Modified for NH3 only

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
    Radiation<TF>& radiation; // Add radiation reference

    bool sw_deposition;
    bool use_depac;          // Switch to toggle between original and DEPAC models
    void sync_depac_fields(Boundary<TF>& boundary);

    std::shared_ptr<Boundary_surface_lsm<TF>> boundary_surface_lsm;

    // Base parameters
    TF deposition_var;
    TF henry_so2;
    TF rsoil_so2;
    TF rwat_so2;
    TF rws_so2;

    // Resistance and physical parameters
    std::vector<TF> rmes;
    std::vector<TF> rsoil;
    std::vector<TF> rcut;
    std::vector<TF> rws;
    std::vector<TF> rwat;
    std::vector<TF> diff;
    std::vector<TF> diff_scl;
    std::vector<TF> henry;
    std::vector<TF> f0;

    // NH3 specific parameters
    TF vd_nh3;            // NH3 deposition velocity

    // DEPAC configuration - Environmental parameters
    // Time-dependent radiation parameters
    TF t0;               // Start time (s)
    TF td;               // Day length (s)
    TF max_rad;          // Maximum radiation (W/m2)
    TF glrad;            // Global radiation (W/m2)
    TF sinphi;           // Sine of solar elevation
    TF temperature;      // Air temperature (K)
    TF rh;               // Relative humidity (%)
    TF sai;              // Stem area index (m2/m2)
    TF lat;              // Latitude (degrees)


    // Time and surface parameters
    int day_of_year;     // Day of year
    int nwet;            // Surface wetness indicator
    int nwet_veg;
    int nwet_soil;
    int nwet_wet;
    int lu;              // Land use type

    // Chemical parameters
    int iratns;          // NH3 compensation point option
    TF hlaw;             // Henry's law constant for NH3
    TF react;            // Reactivity factor
    TF c_ave_prev_nh3;   // Previous NH3 concentration (μg/m3)
    TF catm;             // Atmospheric NH3 concentration (μg/m3)
    TF c_ug;  // Synchronized conversion factor for NH3
    TF pressure;        // Pressure (Pa)

    // Tile management
    std::vector<std::string> deposition_tile_names {"veg", "soil", "wet"};
    Deposition_tile_map<TF> deposition_tiles;

    const std::string tend_name = "deposition";
    const std::string tend_longname = "Deposition";

    bool sw_override_ccomp;       // Flag to override compensation point
    TF ccomp_override_value;      // Value to use when overriding

    // New arrays for grid-mean values
    std::vector<TF> ra_mean;      // Grid-mean aerodynamic resistance
    std::vector<TF> rb_mean;      
    std::vector<TF> obuk_mean;
    std::vector<TF> ustar_mean;
    std::vector<TF> ccomp_mean;   // Grid-mean compensation point
    std::vector<TF> cw_mean;      // Grid-mean external leaf resistance
    std::vector<TF> cstom_mean;   // Grid-mean stomatal resistance
    std::vector<TF> csoil_eff_mean; // Grid-mean soil effective resistance

    // New grid-mean arrays for compensation points
    std::vector<TF> cw_out_mean;   // Grid-mean external leaf compensation point
    std::vector<TF> cstom_out_mean; // Grid-mean stomatal compensation point
    std::vector<TF> csoil_out_mean; // Grid-mean soil compensation point
    std::vector<TF> rc_tot_mean;   // Grid-mean total canopy resistance
    std::vector<TF> rc_eff_mean;   // Grid-mean effective canopy resistance
};
#endif
