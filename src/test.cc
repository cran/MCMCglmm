#include "test.h"


extern "C"{  

/**************************/
/* Meuwissen and Luo 1992 */
/**************************/

void test(
        double *rbv,     
 	double *xGP,
        int *nGP,
	int *nzmaxGP,
        int *n,
        int *nG
){         


  int dimG = nGP[0];
  int cnt, i,j;

  cs *Ginv, *Grv;
  csn *GinvL;
  css *GinvS;


  Ginv = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
  Grv = cs_spalloc(1, dimG, dimG, true, false);

  cnt = 0;
  for(i = 0 ; i < dimG ; i++){
    Ginv->p[i] = i*dimG;    
    Grv->p[i] = i;
    Grv->i[i] = 0;
    Grv->x[i] = 1.0;
    for(j = 0 ; j < dimG; j++){
      Ginv->i[cnt] = j;
      Ginv->x[cnt] = xGP[cnt];
      cnt++;
    }
  }
  Ginv->p[dimG] = dimG*dimG;
  Grv->p[dimG] = dimG;

  GinvS = cs_schol(0, Ginv);                 // Symbolic factorisation of G
  GinvL = cs_chol(Ginv, GinvS);              // cholesky factorisation of G^{-1} for forming N(0, G \otimes I)

  GetRNGstate();
  
  cnt = 0;
  for(j = 0 ; j < nG[0]*n[0]; j++){
    for(i = 0 ; i < dimG ; i++){
      Grv->x[i] = rnorm(0.0,1.0);
    }
    cs_ltsolve(GinvL->L, Grv->x);
    for(i = 0 ; i < dimG ; i++){
      rbv[cnt] = Grv->x[i];
      cnt++;
    }
  }

  PutRNGstate();

  cs_spfree(Ginv);
  cs_spfree(Grv);
  cs_sfree(GinvS);                
  cs_nfree(GinvL);        

}
}
