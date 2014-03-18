#define _RTCMVNORM_H
#include "cs.h"

#include "R.h" 
#include "Rmath.h" 

#ifndef _RTNORM_H
#include "rtnorm.h"
#endif

double rtcmvnorm(const cs *predi, const cs *linki, const cs *G, int keep, double lower, double upper);
