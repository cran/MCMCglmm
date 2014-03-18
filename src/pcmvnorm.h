#define _PCMVNORM_H
#include "cs.h"

#include "R.h" 
#include "Rmath.h" 

#ifndef _CS_INV_H
#include "cs_inv.h"
#endif

double pcmvnorm(const cs *predi, const cs *linki, const cs *G, int keep, double lower, double upper);

