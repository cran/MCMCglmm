#define _CS_OMEGA_H
#include "cs.h"

cs *cs_omega(cs **KGinv, int nG, cs *pvB);
void cs_omegaupdate(cs **KGinv, int nG, cs *pvB, const cs *C);
