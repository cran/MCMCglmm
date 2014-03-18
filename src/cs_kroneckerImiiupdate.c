#include "cs_kroneckerImii.h"

void cs_kroneckerImiiupdate(const cs *A, const cs *B, int nI, double *mii, const cs *C){

    int i, j, k, cnt, anz, am, an;
    cs *AB, *ABinv;
                  
    am = A->m ; an = A->n ; anz = A->nzmax;

    if(am>1){
      AB = cs_spalloc (am, an, anz, 1, 0) ;	 /* allocate result */
      ABinv = cs_spalloc (am, an, anz, 1, 0);

      cnt = 0;

      for (i = 0 ; i < an; i++){
        AB->p[i] = i*an;
        ABinv->p[i] = i*an;
        for (j = 0 ; j <an; j++){
          AB->i[cnt] = j;
          ABinv->i[cnt] = j;
          AB->x[cnt] = 0.0;
          ABinv->x[cnt] = 0.0;
          cnt++;
        }
      }
      AB->p[an] = an*an;
      ABinv->p[an] = an*an;
    }


    if(am>1){

       for(j = 0 ; j < nI ; j++){
        cnt = 0;
        for(i = 0; i < an; i++){
          for(k = 0; k < an; k++){
            AB->x[cnt] = A->x[cnt]+B->x[cnt]*mii[j];
            cnt ++;
          }
        }

        cs_invR(AB, ABinv);


        for(i = 0; i < an; i++){
          for(k = 0; k < am; k++){
            C->x[(i*nI+j)*an+k] = ABinv->x[i*an+k];
          } 
        }
      }

      cs_spfree(AB);
      cs_spfree(ABinv);

    }else{
      for(j = 0 ; j < nI ; j++){
        C->x[j] = 1.0/(A->x[0]+B->x[0]*mii[j]);
      }
    }
}



