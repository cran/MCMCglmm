#include "cs_kroneckerImii.h"

cs *cs_kroneckerImii(const cs *A, const cs *B, int nI, double *mii){

    int i, j, k, cnt, anz, cnz, *Cp,  *Ci, am, an, cm, cn;
    double *Cx;
    cs *C, *AB, *ABinv;

    if (!CS_CSC (A)) return (NULL);                         
    am = A->m ; an = A->n ; anz = A->nzmax;
    cm = am*nI; cn = an*nI; cnz = anz*nI;
    C = cs_spalloc (cm, cn, cnz, 1, 0) ;	 /* allocate result */
    if (!C ) return (cs_done (C, NULL, NULL, 0));   

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


    Cp = C->p ; Ci = C->i ; Cx = C->x;   
    cnt = 0;	
    for (j = 0 ; j < cn ; j++){
         for (i = 0 ; i < am ; i++){
            Ci[cnt] = i*nI+j%nI;
            cnt++;
         }      
    }
    cnt = 0;
    Cp[0] = 0;
    for(i = 0; i<an; i++){
     	for (j = 0 ; j < nI ; j++){
          cnt++;
          Cp[cnt] = Cp[cnt-1]+am;
	}
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
            Cx[(i*nI+j)*an+k] = ABinv->x[i*an+k];
          } 
        }
      }

      cs_spfree(AB);
      cs_spfree(ABinv);

    }else{
      for(j = 0 ; j < nI ; j++){
        Cx[j] = 1.0/(A->x[0]+B->x[0]*mii[j]);
      }
    }

    cs_sprealloc (C, 0) ;		// remove extra space from C 
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}



