#include "cs_dcmvnorm.h"

#define LPIx2 1.837877066409345339082

double cs_dcmvnorm(const cs *beta,  const cs *mu, double ldet, const cs *Minv, const cs *M, int *keep, int nkeep, int *cond, int ncond){

  double llik = 0.0;

  int i,
      j, 
      cnt;
  double tmp;

 
  cs *invSc, *S22, *invS22, *S12, *muC, *dev, *cdev;

  invSc = cs_spalloc (nkeep, nkeep, nkeep*nkeep, 1, 0);
  S22 = cs_spalloc (ncond, ncond, ncond*ncond, 1, 0);
  invS22 = cs_spalloc (ncond, ncond, ncond*ncond, 1, 0);
  S12 = cs_spalloc (nkeep, ncond, nkeep*ncond, 1, 0);
  muC = cs_spalloc (ncond, nkeep, nkeep*ncond, 1, 0);
  dev = cs_spalloc (nkeep, 1, nkeep, 1, 0);
  cdev = cs_spalloc (ncond, 1, ncond, 1, 0);

  cnt = 0;
  for(i = 0; i<nkeep; i++){
    invSc->p[i] = i*nkeep;
    for(j = 0; j<nkeep; j++){
      invSc->i[cnt] = j;
      invSc->x[cnt] = Minv->x[Minv->p[keep[i]]+keep[j]];
      cnt ++;
    }
  }
  invSc->p[nkeep] = nkeep*nkeep;

  cnt = 0;
  for(i = 0; i<ncond; i++){ 
    S22->p[i] = i*ncond;
    for(j = 0; j<ncond; j++){
      S22->i[cnt] = j;
      S22->x[cnt] = M->x[M->p[cond[i]]+cond[j]];
      cnt++;
    }
  }
  S22->p[ncond] = ncond*ncond;

  cnt = 0;
  for(i = 0; i<ncond; i++){
    invS22->p[i] = i*ncond;
    cdev->i[i] = i;
    cdev->x[i] = beta->x[cond[i]]-mu->x[cond[i]];
    for(j = 0; j<ncond; j++){
      invS22->i[cnt] = j;
      invS22->x[cnt] = 1.0;
      cnt++;
    }
  }
  cdev->p[0] = 0;
  cdev->p[1] = ncond;
  invS22->p[ncond] = ncond*ncond;

  ldet -= log(cs_invR(S22, invS22));

  cnt = 0;
  for(i = 0; i<ncond; i++){
    S12->p[i] = i*nkeep;
    for(j = 0; j<nkeep; j++){
      S12->i[cnt] = j;
      S12->x[cnt] = M->x[M->p[cond[i]]+keep[j]]; 
      cnt ++;
    }
  }
  S12->p[ncond] = ncond*nkeep;

  muC = cs_multiply(S12, invS22);

  dev = cs_multiply(muC, cdev);
 
  for(i = 0; i<nkeep; i++){
    dev->x[i]  = beta->x[keep[i]]-dev->x[i]-mu->x[keep[i]];
  }

  for(i=0; i<nkeep; i++){
    tmp= 0.0;
    for(j=0; j<nkeep; j++){
      tmp += dev->x[j]*(invSc->x[j+i*nkeep]);
    }    
    llik-= tmp*dev->x[i];
  }

  llik -= LPIx2*(nkeep);
  llik -= ldet;
  llik /= 2.0;

  cs_spfree(invSc);
  cs_spfree(S22);
  cs_spfree(invS22);
  cs_spfree(S12);
  cs_spfree(muC);
  cs_spfree(dev);
  cs_spfree(cdev);

  return (llik);
}                

