#define _CS_RR_H
#include "cs.h"

#ifndef _CS_INV_H
#include "cs_inv.h"
#endif

#ifndef _CS_COV2COR_H
#include "cs_cov2cor.h"
#endif

#ifndef _CS_RINVWISHART_H
#include "cs_rinvwishart.h"
#endif

#include "R.h" 
#include "Rmath.h" 

cs *cs_rR(const cs *A, double nu, double nuR, const css *As, const cs *Roldinv, double Roldldet, const cs *pG);

