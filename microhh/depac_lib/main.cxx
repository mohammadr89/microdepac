#include <iostream>

// C linkage for Fortran wrapper
extern "C" 
{
   void depac_wrapper(
       const char* compnam,
       int day_of_year,
       float lat,
       float t,
       float ust,
       float glrad,
       float sinphi,
       float rh,
       float lai, 
       float sai,
       int nwet,
       int lu,
       int iratns,
       float *rc_tot,
       float *ccomp_tot,
       float hlaw,
       float react,
       int *status,
       float c_ave_prev_nh3,
       float ra,
       float rb,
       float catm,
       float *rc_eff,
       float *gsoil_eff_out,
       float *rsoil_eff_out
   );
}

// Template functions 
template<typename TF>
TF get_rc_tot()
{
   char compnam[4] = "NH3";
   int day_of_year = 226;
   float lat = 51.0;
   float t = 298.0;
   float ust = 0.2;
   float glrad = 800.0;
   float sinphi = 0.6;
   float rh = 70.0;
   float lai = 5.0;
   float sai = 6.0;
   int nwet = 0;
   int lu = 4;
   int iratns = 2;
   float hlaw = 57.0;
   float react = 0.1;
   float c_ave_prev_nh3 = 6.4;
   float ra = 1.0;
   float rb = 1.0;
   float catm = 10.0;
   float rc_tot;
   float ccomp_tot;
   float rc_eff;
   float gsoil_eff_out;
   float rsoil_eff_out;
   int status;

   depac_wrapper(
       compnam,
       day_of_year,
       lat,
       t,
       ust,
       glrad, 
       sinphi,
       rh,
       lai,
       sai,
       nwet,
       lu,
       iratns,
       &rc_tot,
       &ccomp_tot,
       hlaw,
       react,
       &status,
       c_ave_prev_nh3,
       ra,
       rb,
       catm,
       &rc_eff,
       &gsoil_eff_out,
       &rsoil_eff_out
   );

   std::cout << "Total canopy resistance (rc_tot): " << rc_tot << " s/m" << std::endl;
   
   return static_cast<TF>(rc_tot);
}

template<typename TF>
TF get_rc_eff()
{
   char compnam[4] = "NH3";
   int day_of_year = 226;
   float lat = 51.0;
   float t = 298.0;
   float ust = 0.2;
   float glrad = 800.0;
   float sinphi = 0.6;
   float rh = 70.0;
   float lai = 5.0;
   float sai = 6.0;
   int nwet = 0;
   int lu = 4;
   int iratns = 2;
   float hlaw = 57.0;
   float react = 0.1;
   float c_ave_prev_nh3 = 6.4;
   float ra = 1.0;
   float rb = 1.0;
   float catm = 10.0;
   float rc_tot;
   float ccomp_tot;
   float rc_eff;
   float gsoil_eff_out;
   float rsoil_eff_out;
   int status;

   depac_wrapper(
       compnam,
       day_of_year,
       lat,
       t,
       ust,
       glrad,
       sinphi,
       rh,
       lai,
       sai,
       nwet,
       lu,
       iratns,
       &rc_tot,
       &ccomp_tot,
       hlaw,
       react,
       &status,
       c_ave_prev_nh3,
       ra,
       rb,
       catm,
       &rc_eff,
       &gsoil_eff_out,
       &rsoil_eff_out
   );

   std::cout << "Effective total canopy resistance (rc_eff): " << rc_eff << " s/m" << std::endl;
   
   return static_cast<TF>(rc_eff);
}

template<typename TF>
TF get_rsoil_eff()
{
   char compnam[4] = "NH3";
   int day_of_year = 226;
   float lat = 51.0;
   float t = 298.0;
   float ust = 0.2;
   float glrad = 800.0;
   float sinphi = 0.6;
   float rh = 70.0;
   float lai = 5.0;
   float sai = 6.0;
   int nwet = 0;
   int lu = 4;
   int iratns = 2;
   float hlaw = 57.0;
   float react = 0.1;
   float c_ave_prev_nh3 = 6.4;
   float ra = 1.0;
   float rb = 1.0;
   float catm = 10.0;
   float rc_tot;
   float ccomp_tot;
   float rc_eff;
   float gsoil_eff_out;
   float rsoil_eff_out;
   int status;
   depac_wrapper(
       compnam,
       day_of_year,
       lat,
       t,
       ust,
       glrad,
       sinphi,
       rh,
       lai,
       sai,
       nwet,
       lu,
       iratns,
       &rc_tot,
       &ccomp_tot,
       hlaw,
       react,
       &status,
       c_ave_prev_nh3,
       ra,
       rb,
       catm,
       &rc_eff,
       &gsoil_eff_out,
       &rsoil_eff_out
   );

   std::cout << "Effective soil resistance (rsoil_eff): " << rsoil_eff_out << " s/m" << std::endl;
   
   return static_cast<TF>(rsoil_eff_out);
}

int main()
{
   auto rc_tot = get_rc_tot<float>();
   auto rc_eff = get_rc_eff<float>();
   auto rsoil_eff = get_rsoil_eff<float>();
   
   return 0;
}
