#define _DCUTPOINTS_H
#include "cs.h"

#include "R.h" 
#include "Rmath.h" 


double dcutpoints(const cs *liab, double *yP, int *observed, int start,int finish, double *oldcutopints, double *newcutopints, int stcutpoints, int ncutpoints, double sdcp, double sdl);
