#define _CS_DCMVNORM_H
#include "cs.h"

#include "R.h" 
#include "Rmath.h" 

#ifndef _CS_INV_H
#include "cs_inv.h"
#endif

double cs_dcmvnorm(const cs *beta,  const cs *mu, double ldet, const cs *Minv, const cs *M, int *keep, int nkeep, int *cond, int ncond);

