#include "cs_rR.h"

cs *cs_rR(const cs *A, double nu, double nuR, const css *As, const cs *Roldinv, double Roldldet){
    
	cs *Rnew, *Rnewinv, *Ainv;
	double Rnewldet, MH;
    int dimG = A->n;
	int cnt = 0;
	int i, j;
	
	Rnewinv = cs_spalloc (dimG, dimG, dimG*dimG, 1, 0);
	
	for (i = 0 ; i < dimG; i++){
	  Rnewinv->p[i] = i*dimG;
	  for (j = 0 ; j < dimG; j++){
		Rnewinv->i[cnt] = j;
		Rnewinv->x[cnt] = 0.0;
		cnt++;
	  }
	}
	Rnewinv->p[dimG] = dimG*dimG;
		
	cs_cov2cor(A);
	Ainv = cs_inv(A);
	
	Rnew = cs_rinvwishart(Ainv, nu, As);	
	cs_cov2cor(Rnew);
		
	Rnewldet = log(cs_invR(Rnew, Rnewinv));
		
	MH = Roldldet-Rnewldet;

	for (i = 0 ; i < dimG; i++){
		MH += log(Roldinv->x[i*dimG+i]);	
		MH -= log(Rnewinv->x[i*dimG+i]);
	}
	
	MH *= 0.5*nuR;
	
	if(MH<log(runif(0.0,1.0))){
		cs_invR(Roldinv, Rnew);	// save old R
		for (i = 0 ; i < dimG; i++){
			Rnew->x[i*dimG+i] = 1.0;
        }
	}	
	
    cs_spfree(Rnewinv);
    cs_spfree(Ainv);

    return (cs_done (Rnew, NULL, NULL, 1)) ;	/* success; free workspace, return C */

}


