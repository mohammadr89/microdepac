#ifndef DEPAC_WRAPPER_H
#define DEPAC_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif
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
                   float c_ave_prev_nh3);

#ifdef __cplusplus
}
#endif

#endif // DEPAC_WRAPPER_H

