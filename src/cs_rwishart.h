#define _CS_RWISHART_H
#include "cs.h"

#ifndef _CS_INV_H
#include "cs_inv.h"
#endif

#include "R.h" 
#include "Rmath.h" 

cs *cs_rwishart(const cs *A, double nu, const css *As);

