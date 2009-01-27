#define _CS_RCINVWISHART_H
#include "cs.h"

#ifndef _CS_RINVWISHART_H
#include "cs_rinvwishart.h"
#endif

#ifndef _CS_KRONECKERA_H
#include "cs_kroneckerA.h"
#endif

#ifndef _CS_INV_H
#include "cs_inv.h"
#endif

#include "R.h" 
#include "Rmath.h" 

cs *cs_rCinvwishart(const cs *A, double nu, int split);

/* This Function samples from the conditional inverse wishart where A is the inverse sum of squares and split defines the matrix partion. For example: */
/*                                                                                                                                                     */
/*       1  0  0                                                                                                                                       */
/*   A = 0  1 0.5    and split =1  then A11 = 1, A12=t(A21) = [0 0], and A22 =  1 0.5   and A11 and A12 are sampled conditioning on A22                */
/*       0 0.5 1                                                                0.5 1                                                                  */
/*                             alternatively if                                                                                                        */
/*       1  0  0                                                                                                                                       */
/*   A = 0  1 0.5    and split =2  then A11 = 1 0, A12=t(A21) = [0 0], and A22 =  1     and A11 and A12 are sampled conditioning on A22                */
/*       0 0.5 1                              0 1                                                                                                      */
/*                                                                                                                                                     */
/* The rinvwishart samples directly from A where As is the symbolic cholesky factoristion of A.                                                        */
/* For rCinvwishart  As needs to be the symbolic cholesky factoristion of A11 and T2s is the symbolic cholesky factoristaion of var(T2) from G&S       */
/* Note that for sampling a covariance matrix C with C22=F then A22 needs to be (F*nu)^{2}                                                             */
/*                                                                                                                                                     */
/*******************************************************************************************************************************************************/
