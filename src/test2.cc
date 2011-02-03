#include "test.h"


extern "C"{  

/**************************/
/* Meuwissen and Luo 1992 */
/**************************/

void test2(
        double *rbv,     
 	double *xGP,
        int *nGP,
	int *nzmaxGP,
        int *n,
        int *nG
){         


  int dimG = nGP[0];
  int cnt, i,j;

  cs *Ginv, *KGinv, *Grv;
  csn *GinvL;
  css *GinvS;


  Ginv = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
  Grv = cs_spalloc(1, dimG*nG[0], dimG*nG[0], true, false);

  cnt = 0;
  for(i = 0 ; i < dimG ; i++){
    Ginv->p[i] = i*dimG;    
    for(j = 0 ; j < dimG; j++){
      Ginv->i[cnt] = j;
      Ginv->x[cnt] = xGP[cnt];
      cnt++;
    }
  }
  Ginv->p[dimG] = dimG*dimG;

  cnt = 0;
  for(i = 0 ; i < dimG*nG[0]; i++){
    Grv->p[i] = i;
    Grv->i[i] = 0;
    Grv->x[i] = 1.0;
  }
  Grv->p[dimG*nG[0]] = dimG*nG[0];

  double me = 1.001;
  double *me_ptr = &me;

  KGinv = cs_kroneckerI(Ginv,nG[0]); 

  for(i = 0 ; i < (dimG*dimG*nG[0]); i++){
      KGinv->x[i] = *me_ptr;
  }

  cs_print(KGinv, true);

  *me_ptr = 2.0;

  cs_print(KGinv, true);


  GinvS = cs_schol(0, KGinv);                 // Symbolic factorisation of G
  GinvL = cs_chol(KGinv, GinvS);              // cholesky factorisation of G^{-1} for forming N(0, G \otimes I)

  GetRNGstate();
  
  cnt = 0;
  for(j = 0 ; j < n[0]; j++){
    for(i = 0 ; i < dimG*nG[0] ; i++){
      Grv->x[i] = rnorm(0.0,1.0);
    }
    cs_ltsolve(GinvL->L, Grv->x);
    for(i = 0 ; i < dimG*nG[0] ; i++){
      rbv[cnt] = Grv->x[i];
      cnt++;
    }
  }

  PutRNGstate();

  cs_spfree(Ginv);
  cs_spfree(KGinv);
  cs_spfree(Grv);
  cs_sfree(GinvS);                
  cs_nfree(GinvL);        

}
}
