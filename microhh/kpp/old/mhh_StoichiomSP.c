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
/* File                 : mhh_StoichiomSP.c                         */
/* Time                 : Tue Oct 18 15:45:13 2022                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : mhh.kpp                                   */
/* Output root filename : mhh                                       */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mhh_Parameters.h"
#include "mhh_Global.h"
#include "mhh_Sparse.h"



/* Row-compressed sparse data for the Jacobian of reaction products JVRP */
 /* Beginning of rows in JVRP */

  int  CROW_JVRP[] = {
       0,  2,  4,  6,  7,  9, 10, 12, 14, 16, 18, 19,
      21, 23, 25, 27, 28, 29, 31, 33, 35, 37, 39, 41,
      43, 45, 46, 48, 50, 52, 53, 54, 55, 56, 57, 58,
      59, 60, 60, 60, 61, 62, 63, 64, 65, 66, 67 }; 

 /* Column indices in JVRP */

  int  ICOL_JVRP[] = {
       8,  9,  8, 10,  9, 10, 10,  0,  9,  9,  8, 13,
       7,  8, 11, 13,  7, 11,  1, 10, 13,  7,  9, 10,
      11,  3,  9,  1,  9, 10, 12, 10, 12, 12, 13, 11,
      12,  4,  9,  6,  9,  6, 11,  2,  9, 12,  5,  8,
       5,  9,  5, 11,  8,  7,  1, 11,  4,  6,  6,  0,
       8, 13,  7,  3,  0,  6,  4 }; 

 /* Row indices in JVRP */

  int  IROW_JVRP[] = {
       0,  0,  1,  1,  2,  2,  3,  4,  4,  5,  6,  6,
       7,  7,  8,  8,  9,  9, 10, 11, 11, 12, 12, 13,
      13, 14, 14, 15, 16, 17, 17, 18, 18, 19, 19, 20,
      20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 26, 26,
      27, 27, 28, 28, 29, 30, 31, 32, 33, 34, 35, 36,
      39, 40, 41, 42, 43, 44, 45 }; 




/*  Stoichiometric Matrix in Compressed Column Sparse Format        */

 /* Beginning of columns in STOICM */

  int  CCOL_STOICM[] = {
       0,  3,  6,  8, 10, 13, 15, 18, 21, 24, 27, 30,
      34, 37, 40, 43, 45, 47, 50, 53, 58, 63, 67, 71,
      74, 77, 80, 87, 91, 93, 95, 98,101,104,108,110,
     113,115,116,117,118,119,120,121,122,123,124 }; 

 /* Row indices in STOICM */

  int  IROW_STOICM[] = {
       8,  9, 10,  8,  9, 10,  9, 10,  0, 10,  0,  9,
      10,  9, 10,  7,  8, 13,  7,  8, 11,  7, 11, 13,
       1,  7, 11,  1,  7, 11,  7,  9, 10, 13,  3,  7,
       9,  3, 10, 11,  3,  9, 11,  1,  3,  9, 12,  4,
      10, 12,  6, 10, 12,  6,  7, 10, 12, 13,  6,  7,
      10, 11, 12,  4,  6,  9, 12,  2,  6,  9, 10,  3,
       6, 11,  2,  9, 10,  6, 10, 12,  2,  5,  6,  8,
       9, 10, 12,  5,  6,  9, 12,  5, 11,  8,  9,  7,
       8, 13,  1,  7, 11,  7,  8, 11,  4,  6,  9, 10,
       2,  6,  2,  6, 10,  0,  9, 13,  5,  8, 13,  7,
       3,  0,  6,  4 }; 

 /* Column indices in STOICM */

  int  ICOL_STOICM[] = {
       0,  0,  0,  1,  1,  1,  2,  2,  3,  3,  4,  4,
       4,  5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,
       9,  9,  9, 10, 10, 10, 11, 11, 11, 11, 12, 12,
      12, 13, 13, 13, 14, 14, 14, 15, 15, 16, 16, 17,
      17, 17, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20,
      20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 23,
      23, 23, 24, 24, 24, 25, 25, 25, 26, 26, 26, 26,
      26, 26, 26, 27, 27, 27, 27, 28, 28, 29, 29, 30,
      30, 30, 31, 31, 31, 32, 32, 32, 33, 33, 33, 33,
      34, 34, 35, 35, 35, 36, 36, 37, 38, 39, 40, 41,
      42, 43, 44, 45 }; 

 /* Stoichiometric Matrix in compressed column format */

  double  STOICM[] = {
        -1,   -1,    1,   -1,    1,   -1,   -1,   -1,
         1,   -2,   -1,   -1,    1,   -1,    1,    1,
        -1,   -1,   -1,   -1,    1,    2,   -1,   -1,
         1,   -1,   -1,   -1,    1,    1,    1,    1,
        -1,   -1,    1,   -1,   -1,    1,   -1,   -1,
        -1,   -1,    1,   -1,    2,   -1,    1,    1,
        -1,   -1,    1,   -1,   -1,    1,    1,    1,
        -1,   -1,    1,    1,    1,   -1,   -1,   -1,
       0.4, -0.6,  0.6,    1,   -1,   -1,    1,    1,
        -1,   -1,   -1,   -1,    1, 1.37, 0.74,   -2,
      0.56,   -1, 1.04,   -1, 0.33, 0.19, 0.31,   -1,
       1.5,   -1,    1,   -1,   -1,   -1,    2,   -1,
         1,    1,   -1,    1,    1,    1,    1,   -1,
        -1,    1,    1,    1,    1,   -1,    1,   -1,
         2,   -1,    2,    1,    1,   -1,   -1,   -1,
        -1,   -1,   -1,   -1 }; 


