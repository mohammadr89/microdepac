/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Sparse Hessian Data Structures File                              */
/*                                                                  */
/* Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  */
/*       (http://www.cs.vt.edu/~asandu/Software/KPP)                */
/* KPP is distributed under GPL, the general public licence         */
/*       (http://www.gnu.org/copyleft/gpl.html)                     */
/* (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           */
/* (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            */
/*     With important contributions from:                           */
/*        M. Damian, Villanova University, USA                      */
/*        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany */
/*                                                                  */
/* File                 : tm5_ifs_HessianSP.c                       */
/* Time                 : Fri Nov 19 10:37:06 2021                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : tm5_ifs.kpp                               */
/* Output root filename : tm5_ifs                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tm5_ifs_Parameters.h"
#include "tm5_ifs_Global.h"
#include "tm5_ifs_Sparse.h"



/* Hessian Sparse Data                                              */
/*                                                                  */
 /* Index i of Hessian element d^2 f_i/dv_j.dv_k */

  int  IHESS_I[] = {
       0,  0,  1,  2,  2,  3,  3,  3,  4,  4,  5,  5,
       5,  5,  5,  6,  6,  6,  7,  7,  7,  7,  7,  7,
       8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,
       9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 11, 11,
      11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13,
      13, 13, 13, 13, 13, 13, 13, 13, 13 }; 

 /* Index j of Hessian element d^2 f_i/dv_j.dv_k */

  int  IHESS_J[] = {
       8, 11,  9,  2,  8,  3,  5,  6,  4,  7,  4,  5,
       6,  7,  7,  6,  6,  6,  4,  6,  6,  6,  7,  7,
       2,  3,  5,  6,  7,  7,  8,  8,  8,  8,  8, 10,
       6,  8,  9,  9, 10,  6,  8, 10, 10, 10,  6,  7,
       8,  9,  9, 10, 10, 11,  7,  8,  9, 10,  2,  3,
       4,  5,  6,  6,  8,  8,  8, 10, 11 }; 

 /* Index k of Hessian element d^2 f_i/dv_j.dv_k */

  int  IHESS_K[] = {
       9, 13, 11, 13,  8, 13, 13, 10, 13,  8, 13, 13,
      10,  8, 12,  9, 10, 13, 13,  9, 10, 13,  8, 12,
      13, 13, 13, 10,  8, 12,  8,  9, 10, 12, 13, 13,
       9,  9, 11, 12, 11, 10, 10, 11, 12, 13,  9, 12,
      12, 11, 12, 11, 12, 13, 12, 12, 12, 12, 13, 13,
      13, 13, 10, 13, 10, 12, 13, 13, 13 }; 


