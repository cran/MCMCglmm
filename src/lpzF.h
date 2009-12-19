#define _LPZF_H
#include "cs.h"

#include "R.h" 
#include "Rmath.h" 


double lpzF(double x, int lowertail);

/* log cdf of pnorm with mean=0 and sd=1, using Didemenko's feller approximation for extreme values */


