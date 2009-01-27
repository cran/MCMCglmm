#include "cs_kroneckerH.h"

cs *cs_kroneckerH(const cs *A, int *rl){

    int i, j, n, nI, cnt;

    n = A->n;
    nI = 0;
    cnt = 0;

    for (j = 0 ; j < n ; j++){
      nI+=rl[j];
    }

    cs *C;
    if (!CS_CSC (A)) return (NULL);                         
    C = cs_spalloc (nI, nI, nI, 1, 0) ;	 /* allocate result */
    if (!C ) return (cs_done (C, NULL, NULL, 0));   
//    if (!C ) error("cs_kroneckerH out of memory");   

    for (j = 0 ; j < n ; j++){
      for (i = 0 ; i < rl[j] ; i++){
        C->i[cnt] = cnt;
        C->p[cnt] = cnt;
        C->x[cnt] = A->x[j*n+j];
        cnt++;
      }
    }

    C->p[nI] = nI;

    cs_sprealloc (C, 0) ;		// remove extra space from C 
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}



