#include <iostream>
using std::cin;
using std::cout;
using std::endl;

extern "C" {
void depac_wrapper(const char* compnam,
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
                   float *rc_eff);

// New function to calculate and return vdnh3
float get_vdnh3() {
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
    int status;

    depac_wrapper(
        compnam, day_of_year, lat, t, ust, glrad, sinphi, rh,
        lai, sai, nwet, lu, iratns,
        &rc_tot, &ccomp_tot, hlaw, react,
        &status, c_ave_prev_nh3, ra, rb, catm, &rc_eff
    );

    float vdnh3 = 1.0 / (rc_eff + ra + rb);
    cout << "Deposition velocity (vdnh3): " << vdnh3 << " m/s" << endl;
    return vdnh3;
}
}

int main()
{
    float vdnh3 = get_vdnh3();
    return 0;
}
