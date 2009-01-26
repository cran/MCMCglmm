#define _CS_DMVNORM_H
#include "cs.h"

#include "R.h" 
#include "Rmath.h" 

double cs_dmvnorm(const cs *beta,  const cs *mu, double ldet, const cs *Minv);

