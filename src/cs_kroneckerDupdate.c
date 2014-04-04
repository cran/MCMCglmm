#include "MCMCglmm.h"

void cs_kroneckerDupdate(const cs *A, int nI, double *diag, const cs *C, int reciprocal){

    int i, j, k, cnt, an, am;
    double *Ax;

    an = A->n;
    am = A->m;
    Ax = A->x;
	
    cnt = 0;
 
    if(reciprocal){
      for(i = 0; i < an; i++){
     	for(j = 0 ; j < nI ; j++){
          for(k = 0; k < am; k++){
            C->x[cnt] = Ax[i*an+k]/diag[j];
            cnt++;
          }
        }
      }
    }else{
      for(i = 0; i < an; i++){
     	for(j = 0 ; j < nI ; j++){
          for(k = 0; k < am; k++){
            C->x[cnt] = Ax[i*an+k]*diag[j];
            cnt++;
          }
        }
      }
    }
}



