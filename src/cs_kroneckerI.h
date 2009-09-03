#define _CS_KRONECKERI_H
#include "cs.h"

cs *cs_kroneckerI(const cs *A, int nI);
void cs_kroneckerIupdate(const cs *A, int nI, const cs*C);

