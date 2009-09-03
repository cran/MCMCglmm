#define _CS_DIRECTSUM_H
#include "cs.h"

cs *cs_directsum(cs **KGinv, int nG, int nGR);
void cs_directsumupdate(cs **KGinv, int nG, int nGR, const cs *C);
