#include "cs_tmultiplyupdate.h"
/************************************************************************************************/
/*  A'BW where all matrices must be sorted row within column and all rows/columns have entries  */
/************************************************************************************************/

/*       0 1 0            p = {0, 1, 3, 4}                                                         */
/*   X = 1 0 2    then,   i = {1, 0, 2, 1}                                                         */
/*       0 1 0            x = {1.0, 1.0, 1.0, 2.0}                                                 */
/* row i A' col j B  = col i A col j C  = A'Bij                                                    */

void cs_tmultiplyupdate(const cs *A, const cs *B, const cs *C)
{
	
	
	int i, j, k, l, m, n, cnt=0, cnt2, *Ci, *Cp, *Ap, *Ai, *Bp, *Bi;
    double  *Ax, *Bx, *Cx, sv;
	
	n = C->n;

	Ap = A->p;
	Ai = A->i;
	Ax = A->x;
	Bp = B->p;
	Bi = B->i;
	Bx = B->x;
	Cp = C->p; 
	Ci = C->i;
	
	for(j = 0; j<n; j++){
		for(m=Cp[j]; m<Cp[j+1]; m++){
			i = Ci[m];
			sv = 0.0;
			cnt2=Ap[i];
			for(k=Bp[j]; k<Bp[j+1]; k++){
				for(l=cnt2; l<Ap[i+1]; l++){
				  if(Bi[k]<=Ai[l]){						  
				    if(Bi[k]==Ai[l]){
				      sv += Bx[k]*Ax[l];
				    }
                  cnt2=l;
				  break;
				  }	 
				}
			}	

 C->x[cnt] = sv;
			cnt++;
		}
	}
		
}
