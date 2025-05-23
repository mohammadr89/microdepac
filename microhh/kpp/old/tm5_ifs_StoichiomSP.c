/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Sparse Stoichiometric Data Structures File                       */
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
/* File                 : tm5_ifs_StoichiomSP.c                     */
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



/* Row-compressed sparse data for the Jacobian of reaction products JVRP */
 /* Beginning of rows in JVRP */

  int  CROW_JVRP[] = {
       0,  2,  4,  6,  7,  9, 10, 12, 14, 16, 17, 19,
      21, 23, 25, 27, 29, 31, 33, 35, 36, 37, 38, 39,
      40, 41, 43, 45, 47, 48, 50, 51, 53, 54, 55, 55,
      55, 56, 57, 58, 59, 60, 61, 62 }; 

 /* Column indices in JVRP */

  int  ICOL_JVRP[] = {
      10, 13,  8, 10,  8, 13,  8,  2, 13, 13, 10, 12,
       8, 12, 11, 13, 13,  7, 12,  7,  8,  7,  8,  4,
      13,  4, 13,  3, 13,  5, 13,  6, 10,  6, 13, 10,
      11,  5,  5,  4,  2, 10, 11,  9, 12,  9, 11,  1,
       8,  9,  1,  6,  9,  9,  1, 10, 12, 11,  0,  4,
       2,  5 }; 

 /* Row indices in JVRP */

  int  IROW_JVRP[] = {
       0,  0,  1,  1,  2,  2,  3,  4,  4,  5,  6,  6,
       7,  7,  8,  8,  9, 10, 10, 11, 11, 12, 12, 13,
      13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19,
      20, 21, 22, 23, 24, 25, 25, 26, 26, 27, 27, 28,
      29, 29, 30, 31, 31, 32, 33, 36, 37, 38, 39, 40,
      41, 42 }; 




/*  Stoichiometric Matrix in Compressed Column Sparse Format        */

 /* Beginning of columns in STOICM */

  int  CCOL_STOICM[] = {
       0,  3,  6,  8, 10, 13, 15, 18, 22, 25, 27, 32,
      35, 38, 42, 44, 47, 51, 58, 61, 63, 66, 69, 71,
      75, 77, 80, 83, 86, 89, 92, 94, 98,101,103,104,
     105,106,107,108,109,110,111,112 }; 

 /* Row indices in STOICM */

  int  IROW_STOICM[] = {
       8, 10, 13,  8, 10, 13,  8, 13,  2,  8,  2,  8,
      13,  8, 13, 10, 11, 12,  8, 11, 12, 13,  0, 11,
      13,  7, 13,  5,  7,  8, 11, 12,  4,  7,  8,  5,
       7,  8,  4,  5,  7, 13,  4, 13,  3,  8, 13,  3,
       5,  8, 13,  3,  5,  6,  7,  8, 10, 13,  6,  7,
      13, 10, 13, 10, 11, 12,  3,  5,  8,  3,  5,  4,
       5,  8, 13,  2, 13,  9, 10, 11,  9, 11, 12,  1,
       9, 11,  1,  9, 11,  0,  8,  9,  0,  1,  6,  7,
       9, 11,  9, 10, 11,  1, 11,  6, 12, 10, 12, 11,
       0,  4,  2,  5 }; 

 /* Column indices in STOICM */

  int  ICOL_STOICM[] = {
       0,  0,  0,  1,  1,  1,  2,  2,  3,  3,  4,  4,
       4,  5,  5,  6,  6,  6,  7,  7,  7,  7,  8,  8,
       8,  9,  9, 10, 10, 10, 10, 10, 11, 11, 11, 12,
      12, 12, 13, 13, 13, 13, 14, 14, 15, 15, 15, 16,
      16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 18, 18,
      18, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 23,
      23, 23, 23, 24, 24, 25, 25, 25, 26, 26, 26, 27,
      27, 27, 28, 28, 28, 29, 29, 29, 30, 30, 31, 31,
      31, 31, 32, 32, 32, 33, 33, 34, 35, 36, 37, 38,
      39, 40, 41, 42 }; 

 /* Stoichiometric Matrix in compressed column format */

  double  STOICM[] = {
         1,   -1,   -1,   -1,   -1,    1,   -1,   -1,
         1,   -2,   -1,    1,   -1,    1,   -1,   -1,
         1,   -1,   -1,    1,   -1,    1,    1,   -1,
        -1,    1,   -1,    1,   -1,    1,    1,   -1,
         1,   -1,   -1,    1,   -1,   -1,   -1,  0.4,
       0.6, -0.6,   -1,-0.77,   -1,    1,   -1,    1,
        -1,    1,   -1, 0.56, 0.54,   -1, 0.31, 0.19,
        -1, 0.33,   -1,    1,   -1,   -1,    2,    1,
        -1,    1,    1,   -1,    2,    1,   -1,   -1,
         1,    1,    1,   -1,    2,    1,   -1,   -1,
        -1,    2,   -1,    1,   -1,   -1,   -1,    1,
         1,    1,   -1,   -1,    2,   -1,   -1,    1,
        -1,    1,   -1,    1,    1,   -1,    2,    1,
         1,   -1,   -1,   -1,   -1,   -1,   -1,   -1
      }; 


