#define _CS_KRONECKERIMII_H
#include "cs.h"

#ifndef _CS_INV_H
#include "cs_inv.h"
#endif

cs *cs_kroneckerImii(const cs *A, const cs *B, int nI, double *mii);
void cs_kroneckerImiiupdate(const cs *A, const cs *B, int nI, double *mii, const cs*C);

