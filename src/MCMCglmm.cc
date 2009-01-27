#include "MCMCglmm.h"

extern "C"{  

/***************************************************************************************************/
/* Matrices are stored in compressed-column format where i indicates the row indicies of non-zero  */
/* values, p indexes the first i of each column, and x the actual non-zero values.  For example,if */
/*                                                                                                 */
/*       0 1 0            p = {0, 1, 3, 4}                                                         */
/*   X = 1 0 2    then,   i = {1, 0, 2, 1}                                                         */
/*       0 1 0            x = {1.0, 1.0, 1.0, 2.0}                                                 */
/*                                                                                                 */
/* where the final element of p is always length(i).                                               */
/* dim is a vector with the number of rows and columns, and nzmax the max number of non-zero values*/
/***************************************************************************************************/

void MCMCglmm(
        double *yP,      // read in data
        double *y2P,     // additional data variables (censoring points, truncation points, binomial n etc)
        double *liabP,   // starting liabilities (= y if gaussian, and no missing data)
        int *mvtype,     // missing data types 
        int *nyP,        // number of records
        int *iXP,        // X       
	int *pXP,	         
	double *xXP,	      
	int *dimXP,	         
	int *nzmaxXP,
        int *iZP,        // Z       
	int *pZP,	         
	double *xZP,	         
	int *dimZP,	         
	int *nzmaxZP,
        int *iAP,        // A       
	int *pAP,	         
	double *xAP,	      
	int *dimAP,	         
	int *nzmaxAP,
        double *MSsdP,       // Mendelian sampling standard deviation for each individual (*1.5 if an individual doesn't have both parents)
        int *idP,            // idenitifier note this is the position in the pedigree of analysed indiviuals
        int *damP,           // dam idenitifier     
        int *sireP,          // sire identifier
        int *PedDimP,        // dimensions of pedigree (#individuals/#columns #columns=2 for normal pedigree or 1 for phylogeny/clones)
        int *AtermP,         // boolean: 1 if kth G is associated with A
        int *GRdimP,         // GR structure dimensions
        int *levelsRP,       // number of GR levels 
        int *updateP,
        int *splitP, 
        int *nGP,            // number of G structures and number of R structures  
        double *GRinvP,      // GR starting values 
        double *GRvpP,       // GR prior V 
        double *GRnpP,       // GR prior df:
        int *nittP,          // Number of MCMC iterations	
        int *thinP,          // Thinning interval
        int *burninP,        // Number of burnin iterations	
        int *prP,            // Should posterior distribution of random effects and liabilities be stored	
        double *LocP,        // Posterior distribution of MME solutions (fixed and random)
        double *VarP,        // Posterior distribution of (co)variances
        double *PLiabP,      // Posterior distribution of liabilities
        int *familyP,        // distribution of response variables
        double *propCP,      // proposal distribition for liabilities
        bool *verboseP,      // print iterations and MH acceptance ratio to screen
        double *BvpP,        // inverse prior covariance matrix for beta        
        double *BmupP,       // prior mean vector for beta     
        int *mfacP,          // vector of J-1 levels for each multinomial response 
	int *observedP,	     // vector of 1 (observed) and 0 (missing)
	int *diagP,          // is a us matrix in fact diagonal?
        int *AMtuneP          // should adaptive Metroplis algorithm be used
){         

int     i, j, k,l,cnt,cnt2,itt,record,dimG,

        nG = nGP[0],          // number of G structures
        nR = nGP[1],          // number of R structures 
        nGR = nG+nR,          // number of variance structures
        ny = nyP[0],          // number of records
        *nlGR = levelsRP,     // number of random levels 
        *GRdim = GRdimP,      // dimensions of variance structures

        nitt = nittP[0], 
        thin = thinP[0], 
        burnin = burninP[0], 

        post_cnt = 0,
        tvc =0; // total number of (co)variance components 

// Multinomial counters 
 
int     nthmn = 0,        // positionin mfac of  the nth multinomial
        nthmnl = 1;       // counter for the lth level of nth multinomial

double   mndenom1 = 1.0,  // sum( exp(l_{j}) for all j) and l_old
         mndenom2 = 1.0;  // sum( exp(l_{j}) for all j) and l_new

int     nrowX =  dimXP[0], ncolX =  dimXP[1],  nzmaxX = nzmaxXP[0]; 
int     nrowZ =  dimZP[0], ncolZ =  dimZP[1],  nzmaxZ = nzmaxZP[0];
int     nrowA =  dimAP[0], ncolA =  dimAP[1],  nzmaxA = nzmaxAP[0];

int 	dimAS =  ncolX+ncolZ;
int     nMH = 0;
       
bool    pr = prP[0];
bool    pl = prP[1];
bool    Aexists = FALSE;
bool    missing = FALSE;

        for(i=0; i<nG; i++){if(AtermP[i]==1){Aexists=TRUE;}}

        cnt2=0;
        for(k=nG; k<nGR; k++){
          for(i=0; i<nlGR[k]; i++){
            if(mvtype[i+cnt2]!=2){
              missing=TRUE;
              if(mvtype[i+cnt2]!=1){
                nMH++; 
              }
            }
          }
          cnt2+=nlGR[k];
        }

double  densityl1,
        densityl2,
        Eaccl= 0.0,
        alpha_star = (-double(GRdimP[nG])/(1.0-2.75*double(GRdimP[nG])))-0.133,  // optimal acceptance ratio
        rACCEPT = 0.9,
        qACCEPT = 2.0,
        t[nGR],
        sd[nGR],
        wn[nGR],
        zn[nGR],
        ldet[nR+nG],
        interval,
        remainder,
        dm = 1.0/PedDimP[1];  // 0.5 for pedigrees 1 for phylogenies for taking averaging of potential bv.

        for(k=0; k<nGR; k++){
           t[k]=0.0;
          sd[k]=0.0;
          wn[k]=0.0;
          zn[k]=0.0;
        }

cs      *X, *Z, *W,  *Wt, *KRinv, *WtmKRinv, *M, *Omega, *MME, *zstar, *zstar_tmp, *zstar_tmp2, *astar, *astar_tmp, *location, *location_tmp, *linky,  *pred, *dev, *linki, *linki_tmp, *predi, *A, *bv, *bv_tmp, *bvA, *bvAbv, *tbv, *pvB, *pmuB;

csn	*L;
css     *S;

cs*     Ginv[nGR];
cs*     propC[nGR];
cs*     propCinv[nGR];
cs*     muC[nGR];
cs*     G[nGR];
cs*     pG[nGR];
cs*     Gtmp[nGR];
cs*     Grv[nGR];
css*    GinvS[nGR];
csn*    GinvL[nGR];
css*    propCinvS[nR];
csn*    propCinvL[nR];
cs*     KGinv[nGR];

// read in fixed-effects design matrix X 

        X = cs_spalloc(nrowX, ncolX, nzmaxX, true, false); /* X =  0  X_2  0  */
                                                           /*      0   0  X_3 */
        for (i = 0 ; i < nzmaxX ; i++){
          X->i[i] = iXP[i];
          X->x[i] = xXP[i];
        }
        for (i = 0 ; i <= ncolX ; i++){
          X->p[i] = pXP[i];
        }

// read in prior covaraince matrix and mean vector for fixed effects.

       pvB = cs_spalloc(ncolX, ncolX, ncolX*ncolX, true, false);
       pmuB = cs_spalloc(ncolX, 1, ncolX, true, false);

        cnt = 0;
        for (i = 0 ; i < ncolX ; i++){
               pvB->p[i] = i*ncolX; 
               pmuB->i[i] = i;
               pmuB->x[i] = BmupP[i];
            for (j = 0 ; j < ncolX ; j++){
               pvB->i[cnt] = j;
               pvB->x[cnt] = BvpP[cnt];
               cnt++;
            }
        }
        pvB->p[ncolX] = ncolX*ncolX;
        pmuB->p[0] = 0; 
        pmuB->p[1] = ncolX;
 
// read in random-effects design matrix Z 
                                                            /*     Za_1 0   0   Zb_1 0   0  */
        Z = cs_spalloc(nrowZ, ncolZ, nzmaxZ, true, false);  /* Z =  0  Za_2 0    0  Zb_2 0  */
                                                            /*     0   0  Za_3  0   0  Zb_3 */
        for (i = 0 ; i < nzmaxZ ; i++){
          Z->i[i] = iZP[i];
          Z->x[i] = xZP[i];
        }
        for (i = 0 ; i <= ncolZ ; i++){
          Z->p[i] = pZP[i];
        }

/*************************************************************/
/* read in inverse numerator matrix A and pedigree/phylogeny */
/*************************************************************/

        for(k=0; k<nG; k++){
          if(AtermP[k]==1){
            dimG = GRdim[k];
            A = cs_spalloc(nrowA, ncolA, nzmaxA, true, false);  
            bv = cs_spalloc(nrowA, dimG, nrowA*dimG, true, false);
            bv_tmp = cs_spalloc(PedDimP[0], dimG, PedDimP[0]*dimG, true, false);

            for (i = 0 ; i < nzmaxA ; i++){
              A->i[i] = iAP[i];
              A->x[i] = xAP[i];
            }
            for (i = 0 ; i <= ncolA ; i++){
              A->p[i] = pAP[i];
            }
       
// create matrix for breeding values when sampling kronecker(G,A) of dimesion ntXpedigree
        
            cnt=0;

            for (i = 0 ; i < dimG; i++){
              bv_tmp->p[i] = cnt;
              for (j = 0 ; j < PedDimP[0] ; j++){
                bv_tmp->i[cnt] = j;
                bv_tmp->x[cnt] = 0.0;
                cnt++;
              }
            }

            bv_tmp->p[dimG] = dimG*PedDimP[0];

// create matrix for breeding values when sampling G (bv'Abv) of dimesion ntXindividuals anlaysed

            cnt=0;

            for (i = 0 ; i < dimG; i++){
              bv->p[i] = cnt;
              for (j = 0 ; j < ncolA ; j++){
                bv->i[cnt] = j;
                bv->x[cnt] = 0.0;
                cnt++;
              }
            }

            bv->p[dimG] = dimG*ncolA;
          }
        }

 
// read in G/R matrices 
 
        for (k = 0 ; k < nGR; k++){
          dimG = GRdim[k];
          Ginv[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          G[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          pG[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          Grv[k] = cs_spalloc(1, dimG, dimG, true, false);
          cnt=0;
           for (i = 0 ; i < dimG; i++){
              Ginv[k]->p[i] = i*dimG;
              G[k]->p[i] = i*dimG;
              pG[k]->p[i] = i*dimG;
              Grv[k]->p[i] = i;
              Grv[k]->i[i] = 0;
              for (j = 0 ; j < dimG; j++){
                 Ginv[k]->i[cnt] = j;
                 Ginv[k]->x[cnt] = GRinvP[cnt+tvc];
                 G[k]->i[cnt] = j;
                 G[k]->x[cnt] = GRinvP[cnt+tvc];
                 pG[k]->i[cnt] = j;
                 pG[k]->x[cnt] = GRvpP[cnt+tvc];
                 cnt++;
              }
          }
          Ginv[k]->p[dimG] = dimG*dimG;
          G[k]->p[dimG] = dimG*dimG;
          pG[k]->p[dimG] = dimG*dimG;
          Grv[k]->p[dimG] = dimG;
          ldet[k] = log(cs_invR(Ginv[k], G[k]));
          Gtmp[k] = cs_inv(Ginv[k]);  //delete this eventually  
          tvc += dimG*dimG; 
        }

/*********************************/
/* Read in proposal distribution */
/*********************************/
        cnt2=0;
        for (k = nG ; k < nGR; k++){
          dimG = GRdim[k];
          propC[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          muC[k] = cs_spalloc(dimG, 1, dimG, true, false);
          cnt=0;
          for (i = 0 ; i < dimG; i++){
            propC[k]->p[i] = i*dimG;
            for (j = 0 ; j < dimG; j++){
              propC[k]->i[cnt] = j;
              propC[k]->x[cnt] = propCP[cnt2];
              cnt++;
              cnt2++;
            }
          }
          propC[k]->p[dimG] = dimG*dimG;
          muC[k]->p[0] = 0;
          muC[k]->p[1] = dimG;
          for (i = 0 ; i < dimG; i++){
            muC[k]->i[i] = i;
            muC[k]->x[i] = 0.0;
          }
          propCinv[k] = cs_inv(propC[k]);
          propCinvS[k] = cs_schol(1, propCinv[k]);   
          propCinvL[k] = cs_chol(propCinv[k], propCinvS[k]);       
        }

// allocate vecotors for pseudo-random effects z* and [0, *a] 

        zstar = cs_spalloc(ny, 1, ny, true, false);
        linky = cs_spalloc(ny, 1, ny, true, false);
        pred = cs_spalloc(ny, 1, ny, true, false);
        dev = cs_spalloc(ny, 1, ny, true, false);

        astar = cs_spalloc(dimAS, 1, dimAS, true, false);
        location = cs_spalloc(dimAS, 1, dimAS, true, false);
        location_tmp = cs_spalloc(dimAS, 1, dimAS, true, false);

        linki = cs_spalloc(dimG, 1, dimG, true, false);
        linki_tmp = cs_spalloc(dimG, 1, dimG, true, false);
        predi = cs_spalloc(dimG, 1, dimG, true, false);                 // This is dimG because this is the diemsnion of the R structures (all of which are the same)

        for (i = 0 ; i < ny; i++){   
           zstar->i[i] = i;
           linky->i[i] = i;
           pred->i[i] = i;
           dev->i[i] = i;
           linky->x[i] = liabP[i];                         /* this needs to be changed for random regression */
        }
        zstar->p[0] = 0; 
        zstar->p[1] = ny;
        linky->p[0] = 0; 
        linky->p[1] = ny;
        pred->p[0] = 0; 
        pred->p[1] = ny;
        dev->p[0] = 0; 
        dev->p[1] = ny;

        for (i = 0 ; i < (dimAS-ncolX); i++){
           astar->i[i] = i+ncolX;    
        }
        astar->p[0] = 0;
        astar->p[1] = dimAS-ncolX;

        for (i = 0 ; i < dimAS; i++){
           location->i[i] = i;    
           location_tmp->i[i] = i;  
           location->x[i] = 0.0;    
           location_tmp->x[i] = 0.0;      
        }
        location->p[0] = 0;
        location->p[1] = dimAS;
        location_tmp->p[0] = 0;
        location_tmp->p[1] = dimAS;

        for (i = 0 ; i < dimG; i++){
           linki->i[i] = i;    
           linki_tmp->i[i] = i;    
           predi->i[i] = i;    
        }
        linki->p[0] = 0;
        linki->p[1] = dimG;
        linki_tmp->p[0] = 0;
        linki_tmp->p[1] = dimG;
        predi->p[0] = 0;
        predi->p[1] = dimG;

// form W = [X, Z] 

        if(nG>0){
          W = cs_cbind(X,Z);
        }else{
          W = cs_spalloc(nrowX, ncolX, nzmaxX, true, false);                                                           
          for (i = 0 ; i < nzmaxX ; i++){
            W->i[i] = iXP[i];
            W->x[i] = xXP[i];
          }
          for (i = 0 ; i <= ncolX ; i++){
            W->p[i] = pXP[i];
          }
        }
        Wt = cs_transpose(W, true);

        for (k = 0 ; k < nGR ; k++){
           GinvS[k] = cs_schol(1, Ginv[k]);                    // Symbolic factorisation of G
           GinvL[k] = cs_chol(Ginv[k], GinvS[k]);              // cholesky factorisation of G^{-1} for forming N(0, G \otimes I)
        }
            

// form KGinv = G^{-1} \otimes I    

        for (k = 0 ; k < nGR ; k++){
           if(AtermP[k]==1){
             KGinv[k] = cs_kroneckerA(Ginv[k],A);            //  form kronecker(G^{-1}, A^{-1}) structure
           }else{
             KGinv[k] = cs_kroneckerI(Ginv[k],nlGR[k]);      //  form G^{-1} structure
           }
        }

        KRinv = cs_directproduct(KGinv, nG, nGR);

// form WtmKRinv = W^{t}%*%KRinv    

        WtmKRinv = cs_multiply(Wt, KRinv); 

// form M = WtmKRinv%*%W  
   
        M = cs_multiply(WtmKRinv, W);

// form Omega = bdiag(0, KGinv) 

        Omega = cs_omega(KGinv, nG, pvB);

// form MME = M + Omega; mixed model equations 

        MME = cs_add(M, Omega, 1.0, 1.0); 
        S = cs_schol(1, MME);                            // Symbolic factorisation - only has to be done once
	
  	GetRNGstate();                                   // get seed for random number generation

/**********************************************************************************************************/
/***** MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC *****/
/**********************************************************************************************************/

        for (itt = 0; itt < nitt; itt++){

/***************************/
/* sample location effects */
/***************************/

          if(itt>0){
            for (i = 0 ; i < nGR; i++){
              if(updateP[i]==1){
                cs_spfree(G[i]); 
              }
            }
            cs_spfree(astar_tmp);
	    cs_spfree(zstar_tmp);
            cs_spfree(zstar_tmp2);
            cs_spfree(pred);
            cs_spfree(dev);
	    cs_nfree(L);
            if(Aexists){
              cs_spfree(tbv);               
              cs_spfree(bvA);               
	      cs_spfree(bvAbv);
	    }                                
          }

          for (i = 0 ; i < nGR; i++){
            if(updateP[i]==1){
	      cs_nfree(GinvL[i]);
	      cs_spfree(KGinv[i]);   
            }
          }
          if(missing){	
            for (i = nG ; i < nGR; i++){
	      cs_nfree(propCinvL[i]);
	      cs_spfree(propCinv[i]);   
            }
          }

	  cs_spfree(KRinv);  
	  cs_spfree(WtmKRinv);     
	  cs_spfree(M);                   
	  cs_spfree(Omega);               
	  cs_spfree(MME);  

          for (i = 0 ; i < nGR ; i++){
            if(updateP[i]==1){
              GinvL[i] = cs_chol(Ginv[i], GinvS[i]);                 // cholesky factorisation of G^{-1} for forming N(0, G \otimes I)
            }
          } 

          cnt = 0;

			for(k=0; k<nG; k++){
            dimG = GRdim[k];
            if(AtermP[k]==1){          

/***************/
/*  pedigrees  */
/* phylogenies */
/***************/
 
              if(PedDimP[0]==nlGR[k]){            // All individuals are present - write directly to astar
                for(i=0; i<nlGR[k]; i++){
                  for(j=0; j<dimG; j++){
                    Grv[k]->x[j] = rnorm(0.0,MSsdP[i]);
                  }
                  
                  cs_ltsolve(GinvL[k]->L, Grv[k]->x);

                  for(j=0; j<dimG; j++){
                    astar->x[j*nlGR[k]+i+cnt] = Grv[k]->x[j];     
                  }
                  if(sireP[i]!=-999){  
                    for(j=0; j<dimG; j++){
                      astar->x[j*nlGR[k]+i+cnt] += 0.5*(astar->x[j*nlGR[k]+sireP[i]+cnt]);
                    }  
                  }
                  if(damP[i]!=-999){  
                    for(j=0; j<dimG; j++){
                      astar->x[j*nlGR[k]+i+cnt] += dm*(astar->x[j*nlGR[k]+damP[i]+cnt]);
                    }  
                  }
                }
              }else{    // For situations where the number of individuals < the dimension of the pedigree/phylogeny

/***************/
/*   partial   */
/*  pedigrees  */
/* phylogenies */
/***************/

                for(i=0; i<PedDimP[0]; i++){
                  for(j=0; j<dimG; j++){
                    Grv[k]->x[j] = rnorm(0.0,MSsdP[i]);
                  }
                  
                  cs_ltsolve(GinvL[k]->L, Grv[k]->x);

                  for(j=0; j<dimG; j++){
                    bv_tmp->x[j*PedDimP[0]+i] = Grv[k]->x[j];     
                  }
                  if(sireP[i]!=-999){  
                    for(j=0; j<dimG; j++){
                      bv_tmp->x[j*PedDimP[0]+i] += 0.5*(bv_tmp->x[j*PedDimP[0]+sireP[i]]);
                    }  
                  }
                  if(damP[i]!=-999){  
                    for(j=0; j<dimG; j++){
                      bv_tmp->x[j*PedDimP[0]+i] += dm*(bv_tmp->x[j*PedDimP[0]+damP[i]]);
                    }  
                  }
                }
                for(i=0; i<nlGR[k]; i++){  // need to read into astar
                  for(j=0; j<dimG; j++){
                    astar->x[j*nlGR[k]+i+cnt] = bv_tmp->x[j*PedDimP[0]+idP[i]];
                  }
                }
              }
            }else{

/***********/
/* blocked */
/***********/ 
              for(i=0; i<nlGR[k]; i++){
                for(j=0; j<dimG; j++){
                  Grv[k]->x[j] = rnorm(0.0,1.0);
                }
                cs_ltsolve(GinvL[k]->L, Grv[k]->x);
                for(j=0; j<dimG; j++){
                  astar->x[j*nlGR[k]+i+cnt] = Grv[k]->x[j];
                }
              }
            }   
            cnt += dimG*nlGR[k];                                                           // form [0, a*]
          }

/************/
/* residual */
/************/

          cnt=0;

          for(k=nG; k<nGR; k++){
            dimG = GRdim[k];
            for(i=0; i<nlGR[k]; i++){
              for(j=0; j<dimG; j++){
                Grv[k]->x[j] = rnorm(0.0,1.0);
              }
              cs_ltsolve(GinvL[k]->L, Grv[k]->x);
              for(j=0; j<dimG; j++){
                zstar->x[j*nlGR[k]+i+cnt] = Grv[k]->x[j];
              }
            }
            cnt += dimG*nlGR[k];  
          }    

/******************/
/* form equations */
/******************/

        for (k = 0 ; k < nGR; k++){                             
          if(updateP[k]==1){
            if(AtermP[k]==1){
              KGinv[k] = cs_kroneckerA(Ginv[k],A);              //  form kronecker(G^{-1}, A^{-1}) structure
            }else{
              KGinv[k] = cs_kroneckerI(Ginv[k],nlGR[k]);        //  form G^{-1} structure
            }
          }  
        }

         KRinv = cs_directproduct(KGinv, nG, nGR);

         WtmKRinv = cs_multiply(Wt, KRinv);     

         M = cs_multiply(WtmKRinv, W);                          // form M = WtmKRinv%*%W  

         Omega = cs_omega(KGinv, nG, pvB);                      // form Omega = bdiag(0, KGinv) 

         MME = cs_add(M, Omega, 1.0, 1.0);                      // form MME = M + Omega; mixed model equations 

/************************/
/* sample pseudo vector */
/************************/

         zstar_tmp = cs_multiply(W, astar);                     // Za*
         zstar_tmp2 = cs_add(zstar,zstar_tmp, 1.0 ,1.0);        // z* = N(Za*, R \otimes I)

         for (i = 0 ; i < ny ; i++){
            zstar->x[i] = linky->x[i] - zstar_tmp2->x[i];       // form y - z*
         }

         astar_tmp = cs_multiply(WtmKRinv, zstar);              // WtmKRinv(y - z*)

         for (i = 0 ; i < astar_tmp->nzmax ; i++){
            location->x[astar_tmp->i[i]] = astar_tmp->x[i];
         }



/*************/
/* solve MME */
/*************/

         L = cs_chol(MME, S); 
 
         for (i = 0 ; i < dimAS; i++){
            location_tmp->x[i] = 0.0;
         }

         cs_ipvec (S->pinv, location->x, location_tmp->x, MME->n);	 // x = P*b 
         cs_lsolve(L->L,location_tmp->x);                                // x = L\x 
	 cs_ltsolve (L->L, location_tmp->x);		                 // x = L'\x 
	 cs_pvec (S->pinv, location_tmp->x, location->x, MME->n);        // b = P'*x 

          for (i = 0 ; i < ncolZ ; i++){
            location->x[i+ncolX]+= astar->x[i];
         }

/***********************/
/* sample VCV matrices */
/***********************/

         pred = cs_multiply(W, location);
         dev = cs_add(linky, pred, 1.0, -1.0);    
         cnt2 = ncolX;

	 for(i=0; i< nG; i++){        
           if(updateP[i]==1){            
             dimG = GRdim[i];
             if(AtermP[i]==1){        
               for(j=0; j<(dimG*ncolA); j++){
                 bv->x[j]=location->x[cnt2+j];
               }
               tbv = cs_transpose(bv, true);
               bvA = cs_multiply(tbv, A);
               bvAbv = cs_multiply(bvA, bv);
               for(k=0; k<(dimG*dimG); k++){
                 Gtmp[i]->x[k] = bvAbv->x[k] + pG[i]->x[k];
               }
             }else{
               for(j=0; j<dimG; j++){
                 for(k=j; k<dimG; k++){
                   cnt = j*dimG+k;
                   Gtmp[i]->x[cnt] = 0.0;
                   for(l=0; l<nlGR[i]; l++){
                     Gtmp[i]->x[cnt] += location->x[cnt2+nlGR[i]*j+l]*location->x[cnt2+nlGR[i]*k+l];
                   }
                   Gtmp[i]->x[cnt] += pG[i]->x[cnt];
                 }
               }
               for(j=1; j<dimG; j++){
                 for(k=j; k<dimG; k++){
                   Gtmp[i]->x[j+dimG*k-1] = Gtmp[i]->x[(j-1)*dimG+k];
                 }
               }
             }
             cnt2 += nlGR[i]*dimG;
             if(splitP[i]==-999){  
               cs_invR(Gtmp[i], Ginv[i]);
               G[i] = cs_rinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], GinvS[i]);
             }else{
               for(j=splitP[i]; j<dimG; j++){
                 for(k=splitP[i]; k<dimG; k++){
                   Gtmp[i]->x[j*dimG+k] = pG[i]->x[j*dimG+k];
                 }
               }
               cs_invR(Gtmp[i], Ginv[i]);
               G[i] = cs_rCinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], splitP[i]);
             }	
	     ldet[i] = log(cs_invR(G[i], Ginv[i]));  
           }
         }

/**********************/
/* Sample R Structure */
/**********************/
         cnt2=0;
         for(i=nG; i<nGR; i++){    
            if(updateP[i]==1){ 
             dimG = GRdim[i];
             for(j=0; j<dimG; j++){
               for(k=j; k<dimG; k++){
                 cnt = j*dimG+k;
                 Gtmp[i]->x[cnt] = 0.0;
                 for(l=0; l<nlGR[i]; l++){
                   Gtmp[i]->x[cnt] += dev->x[cnt2+nlGR[i]*j+l]*dev->x[cnt2+nlGR[i]*k+l];
                 }
                  Gtmp[i]->x[cnt] += pG[i]->x[cnt];
               }
             }
             for(j=1; j<dimG; j++){
               for(k=j; k<dimG; k++){
                 Gtmp[i]->x[j+dimG*k-1] = Gtmp[i]->x[(j-1)*dimG+k];
               }
             }
             cnt2 += nlGR[i]*dimG;
	     if(diagP[i]==1){
	       cnt=0;  
	       for(j=0; j<dimG; j++){
		 for(k=0; k<(dimG-1); k++){
		   if(j==k){cnt++;}
		   Gtmp[i]->x[cnt] = 0.0;
		   cnt++;
		 }
	       }	
	     }
             if(splitP[i]==-999){  
               cs_invR(Gtmp[i], Ginv[i]);
               G[i] = cs_rinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], GinvS[i]);
             }else{
               for(j=splitP[i]; j<dimG; j++){
                 for(k=splitP[i]; k<dimG; k++){
                   Gtmp[i]->x[j*dimG+k] = pG[i]->x[j*dimG+k];
                 }
               }
               cs_invR(Gtmp[i], Ginv[i]);
               G[i] = cs_rCinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], splitP[i]);	 
             }	
             if(diagP[i]==1){
	       cnt=0;  
	       for(j=0; j<dimG; j++){
		 for(k=0; k<(dimG-1); k++){
		   if(j==k){cnt++;}
		   G[i]->x[cnt] = 0.0;
		   cnt++;
		 }
	       }	
	     }		
             ldet[i] = log(cs_invR(G[i], Ginv[i]));		 
           }
         }

/***********************/
/* sample liabilities  */   
/***********************/

     if(missing){

       cnt2=0;
       cnt=0;
       for(k=nG; k<nGR; k++){      // Iterate through R-structures
         dimG = GRdim[k];
		   propCinv[k] = cs_inv(propC[k]);
		   propCinvL[k] = cs_chol(propCinv[k], propCinvS[k]);                 // cholesky factorisation of G^{-1} for forming N(0, G \otimes I)

         for(j=0; j<nlGR[k]; j++){     // Iterate through levels
           densityl1 = 0.0;
           densityl2 = 0.0;
 
           nthmn = 0;           //variables for the multinomial
           nthmnl = 1;
           mndenom1 = 1.0;
           mndenom2 = 1.0;

           for(i=0; i<dimG; i++){
             linki_tmp->x[i] = rnorm(0.0,1.0);
           }

           if(mvtype[cnt+j]==1){         // can be Gibbsed

             cs_ltsolve(GinvL[k]->L, linki_tmp->x);  
             for(i=0; i<dimG; i++){
               linky->x[nlGR[k]*i+j+cnt2] = linki_tmp->x[i]+pred->x[nlGR[k]*i+j+cnt2];
             } 
           }
           if(mvtype[cnt+j]<1){         // has to be MHed
             cs_ltsolve(propCinvL[k]->L, linki_tmp->x);
           }

           for(i=0; i<dimG; i++){       // Iterate through fixed levels

             record=cnt2+nlGR[k]*i+j;

             if(mvtype[cnt+j]>0){break;} // has been Gibbsed, or is fully observed and Guassian and therefore knwon

             linki->x[i] = linky->x[record];
             predi->x[i] =  pred->x[record];
             linki_tmp->x[i] += linki->x[i];

             if(familyP[record]==1 && observedP[record]==1){;
               linki_tmp->x[i] = linki->x[i];

             }else{

               if(observedP[record]==1){

                 switch(familyP[record]){

                   case 1:  /* Normal */
                   break;
  
                   case 2:  /* Posisson */
                     densityl1 += dpois(yP[record], exp(linki->x[i]), true);
                     densityl2 += dpois(yP[record], exp(linki_tmp->x[i]), true);
                   break;

                   case 3:  /* Multinomial */
                     mndenom1 += exp(linki->x[i]);
                     mndenom2 += exp(linki_tmp->x[i]);

                     densityl1 += yP[record]*linki->x[i];
                     densityl2 += yP[record]*linki_tmp->x[i];

                     if(mfacP[nthmn]==nthmnl){ 
                       densityl1 -= y2P[record]*log(mndenom1);
                       densityl2 -= y2P[record]*log(mndenom2);
                       nthmnl = 1;
                       mndenom1 = 1.0;
                       mndenom2 = 1.0;
                     }else{
                       nthmn ++;
                       nthmnl++;
                     }

                   break;
     
                   case 4: /* Weibull */
                     densityl1 += dweibull(yP[record], 1.0, 1.0/exp(linki->x[i]), true);
                     densityl2 += dweibull(yP[record], 1.0, 1.0/exp(linki_tmp->x[i]), true);
                   break;
    
                   case 5: /* Exponential */
                     densityl1 += dexp(yP[record], 1.0/exp(linki->x[i]), true);
                     densityl2 += dexp(yP[record], 1.0/exp(linki_tmp->x[i]), true);
                   break;
    
                   case 6: /* Censored Gaussian */
                     if(linki_tmp->x[i]<yP[record] || linki_tmp->x[i]>y2P[record]){
                       if(linki_tmp->x[i]<yP[record]){
                         if(int(interval)%2==1){
                           linki_tmp->x[i] = y2P[record] - remainder;
                         }else{
                           linki_tmp->x[i] = yP[record] + remainder;
                        }
                      }else{
                        interval = (linki_tmp->x[i]-y2P[record])/(y2P[record]-yP[record]);
                        remainder = (linki_tmp->x[i]-y2P[record])-double(int(interval))*(y2P[record]-yP[record]);
                        if(int(interval)%2==1){
                          linki_tmp->x[i] = yP[record] + remainder;
                        }else{
                          linki_tmp->x[i] = y2P[record] - remainder;
                       }
                     }
                   }
		   break;
 
		   case 7: /* Censored Poisson */
	             densityl1 += log(ppois(y2P[record], exp(linki->x[i]), true, false)-ppois(yP[record], exp(linki->x[i]), true, false));
		     densityl2 += log(ppois(y2P[record], exp(linki_tmp->x[i]), true, false)-ppois(yP[record], exp(linki_tmp->x[i]), true, false));
		   break;
						 
                   case 8: /* Censored Weibull */
                     densityl1 += log(pweibull(y2P[record], 1.0, 1.0/exp(linki->x[i]), true, false)-pweibull(yP[record],1.0, 1.0/exp(linki->x[i]), true, false));
                     densityl2 += log(pweibull(y2P[record], 1.0, 1.0/exp(linki_tmp->x[i]), true, false)-pweibull(yP[record],1.0, 1.0/exp(linki_tmp->x[i]), true, false));
                   break;                                            

                   case 9: /* Censored Exponential */
                     densityl1 += log(pexp(y2P[record], 1.0/exp(linki->x[i]), true, false)-pexp(yP[record], 1.0/exp(linki->x[i]), true, false));
                     densityl2 += log(pexp(y2P[record], 1.0/exp(linki_tmp->x[i]), true, false)-pexp(yP[record], 1.0/exp(linki_tmp->x[i]), true, false));
                   break;  
						 
		   case 10: /* Zero-inflated Gaussian */
		   break;
						 
		   case 11: /* Zero-inflated Poisson */
                     if(nthmnl==1){
                       mndenom1 = dpois(yP[record], exp(linki->x[i]), false);  
                       mndenom2 = dpois(yP[record], exp(linki_tmp->x[i]), false);  
		       nthmn ++;
                       nthmnl++;
                     }else{
			mndenom1 *= (1.0-exp(linki->x[i])/(1.0+exp(linki->x[i])));
			mndenom2 *= (1.0-exp(linki_tmp->x[i])/(1.0+exp(linki_tmp->x[i])));
                        if(yP[record]>0.5){
			  mndenom1 += exp(linki->x[i])/(1.0+exp(linki->x[i]));
			  mndenom2 += exp(linki_tmp->x[i])/(1.0+exp(linki_tmp->x[i]));    
			}
                        densityl1 += log(mndenom1);
                        densityl2 += log(mndenom2);
                        nthmnl = 1;
                        mndenom1 = 1.0;
                        mndenom2 = 1.0;
                     }
		   break;
						 
		   case 12: /* Zero-inflated Weibull */ 
                   break;                                            
						 
                   case 13: /* Zero-inflated Exponential */
		   break;   
                 }
               }
             }
           }

           if(mvtype[cnt+j]<1){
			       
             densityl1 += cs_dmvnorm(linki, predi, ldet[k], Ginv[k]);
             densityl2 += cs_dmvnorm(linki_tmp, predi, ldet[k], Ginv[k]);

             zn[k] *= rACCEPT; 
             wn[k] *= rACCEPT;
             wn[k] ++;

             if((densityl2-densityl1)>log(runif(0.0,1.0))){
               for(i=0; i<dimG; i++){
                 linky->x[nlGR[k]*i+j+cnt2] = linki_tmp->x[i];
               }
               zn[k]++;
               Eaccl++;  
             }

/***************/
/* Adaptive MH */
/***************/

             if(AMtuneP[k]==1 && itt<burnin){
               t[k] ++;
               for(i=0; i<dimG; i++){
                 for(l=0; l<dimG; l++){
                  propC[k]->x[i*dimG+l] *= (t[k]-1.0)/t[k];
                  propC[k]->x[i*dimG+l] += muC[k]->x[i]*muC[k]->x[l];
                 }
               }
			   for(i=0; i<dimG; i++){
				  muC[k]->x[i] *= (t[k]-1.0);
					 muC[k]->x[i] += linky->x[nlGR[k]*i+j+cnt2];
					 muC[k]->x[i] /= t[k];
			   }				 
               for(i=0; i<dimG; i++){
                 for(l=0; l<dimG; l++){
                   propC[k]->x[i*dimG+l] -= ((t[k]+1.0)/t[k])*muC[k]->x[i]*muC[k]->x[l];
                   propC[k]->x[i*dimG+l] += linky->x[nlGR[k]*i+j+cnt2]*linky->x[nlGR[k]*l+j+cnt2]/t[k];
                 }
                 propC[k]->x[i*(dimG+1)] += 0.001/t[k];
               }
             }
           }
         }

         if(AMtuneP[k]==1 && itt<burnin){
           for(i=0; i<dimG; i++){
             for(l=0; l<dimG; l++){
			   propC[k]->x[i*dimG+l] += (muC[k]->x[i]*muC[k]->x[l])/t[k];
               propC[k]->x[i*dimG+l] *= sd[k];
             }
           }

		   if(wn[k]>0.0){
             sd[k] = pow(qACCEPT, ((zn[k]/wn[k])-alpha_star));
		   }		 
           zn[k] = 0.0; 
           wn[k] = 0.0;
         }

         cnt2+=nlGR[k]*dimG;
         cnt+=nlGR[k];
       }
     }

/***********************/
/*   store posterior   */
/***********************/

     if(itt%1000 == 0 && verboseP[0]){
       Rprintf("\nMCMC iteration = %i\n",itt);
       if(nMH>0){
         Rprintf("\nE[MH acceptance ratio] = %f\n", Eaccl/(nMH*1000.0));
         Eaccl = 0.0;
       }
     }
     if(itt>=burnin && itt%thin == 0){
       if(pr){
         for (i = 0 ; i < dimAS ; i++){
           LocP[i+post_cnt*dimAS] = location->x[i];
         }
       }else{
         for (i = 0 ; i < ncolX ; i++){
           LocP[i+post_cnt*ncolX] = location->x[i];
         }
       } 
       if(pl){
         for (i = 0 ; i < ny; i++){
           PLiabP[i+post_cnt*ny]  = linky->x[i];
         }
       }
       cnt=0;
       for(i=0; i<nGR; i++){
         dimG = GRdim[i];
         for(j=0; j<(dimG*dimG); j++){
           VarP[cnt+post_cnt*tvc] = G[i]->x[j];
           cnt++;
         }
       }
       post_cnt++;
     }
   }

	PutRNGstate();

        cs_spfree(X);
        cs_spfree(Z);
        cs_spfree(W);
        cs_spfree(Wt);
        cs_spfree(KRinv);
        cs_spfree(WtmKRinv);
        cs_spfree(M);
        cs_spfree(Omega);
        cs_spfree(MME);
        cs_spfree(zstar);
        cs_spfree(zstar_tmp);
        cs_spfree(zstar_tmp2);
        cs_spfree(astar);
        cs_spfree(astar_tmp);
        cs_spfree(location);
        cs_spfree(location_tmp);
        cs_spfree(linky);
        cs_spfree(linki);
        cs_spfree(linki_tmp);
        cs_spfree(pred);
        cs_spfree(predi);
        cs_spfree(dev);                              
        cs_spfree(pvB);
        cs_spfree(pmuB);                                

                                
        if(Aexists){
          cs_spfree(A);               
          cs_spfree(bv);                                                              
          cs_spfree(bv_tmp);                                
          cs_spfree(bvA);               
          cs_spfree(bvAbv);                                                              
          cs_spfree(tbv);                                
        }


        cs_nfree(L);
        cs_sfree(S);

        for(i=0; i<nGR; i++){
	    cs_spfree(Ginv[i]);
	    cs_spfree(G[i]);
	    cs_spfree(Gtmp[i]);
 	    cs_spfree(Grv[i]);
	    cs_spfree(pG[i]);
	    cs_sfree(GinvS[i]);
	    cs_nfree(GinvL[i]);
            cs_spfree(KGinv[i]);
        }


        for(i=nG; i<nGR; i++){
	    cs_spfree(propC[i]);
	    cs_spfree(propCinv[i]);
	    cs_sfree(propCinvS[i]);
	    cs_nfree(propCinvL[i]);
        }
}
}

