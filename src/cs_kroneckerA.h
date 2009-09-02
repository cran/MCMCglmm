#define _CS_KRONECKERA_H
#include "cs.h"

cs *cs_kroneckerA(const cs *G, const cs *A);
void cs_kroneckerAupdate(const cs *G, const cs *A, const cs *C);

