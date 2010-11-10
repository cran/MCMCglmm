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
	int *dimP,	 // dimensions of X,Z,A and L 
        int *nzmaxP,     // number of non-zero's in X,Z,A and L 
        int *iXP,        // X       
	int *pXP,	         
	double *xXP,	               
        int *iZP,        // Z       
	int *pZP,	         
	double *xZP,	                  
        int *iAP,        // A       
	int *pAP,	         
	double *xAP,	               
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
        double *BvpP,        // inverse prior covariance matrix for beta  with final number indicating diagonal (1) or non diagonal (0)      
        double *BmupP,       // prior mean vector for beta     
        int *mfacP,          // vector of J-1 levels for each multinomial response 
	int *observedP,	     // vector of 1 (observed) and 0 (missing)
	int *diagP,          // is a us matrix in fact diagonal?
        int *AMtuneP,        // should adaptive Metroplis algorithm be used
	int *DICP,	     // should DIC be computed
        double *dbarP,
        int *proposal,
        int *ncutpointsP,    // number of cutpoints with -Inf, 0, C1,.... Cncutpoints, Inf
        int *nordinalP,      // number of ordinal traits  
        double *stcutpointsP,// starting vector of cutpoints    
        double *CPP,
        double *AmupP,
        int *iAVpP,          // priors for redundant parameters in parameter expanded models.      
        int *pAVpP,        
        double *xAVpP,     
        int *nzmaxAVpP,        
        int *PXtermsP,
        int *iLXP,              // Design matrix to form Gianola & Sorensen's Lambda         
        int *pLXP,        
        double *xLXP,
        double *lambdaP,        
        double *LvpP,          // prior for structural parameters
        double *LmupP  
  
){         

int     i, j, k,l,p,cnt,cnt2,rterm,itt,record,dimG,nthordinal,

        nG = nGP[0],          // number of G structures
        nR = nGP[1],          // number of R structures 
        nGR = nG+nR,          // number of variance structures
        ny = nyP[0],          // number of records
        *nlGR = levelsRP,     // number of random levels 
        *GRdim = GRdimP,      // dimensions of variance structures
        ncond=0,
        nkeep=0,
        nitt = nittP[0], 
        thin = thinP[0], 
        burnin = burninP[0], 
        post_cnt = 0,
        tvc =0,          // total number of (co)variance components 
        nordinal = nordinalP[0],
        diagR = diagP[0];  // is the us R-structure diagionalised 0-no, 1-yes, 2-yes with identical variances

	int *cond = new int[GRdim[nG]];
	int *keep = new int[GRdim[nG]];

// Multinomial counters 
 
int     nthmnl = 0;  // counter for the lth level of nth multinomial

double   mndenom1 = 1.0,  // sum(exp(l_{j}) for all j) and l_old
         mndenom2 = 1.0,  // sum(exp(l_{j}) for all j) and l_new
         u;               // slice sampling runifs

int     nrowX =  dimP[0], ncolX =  dimP[1],  nzmaxX = nzmaxP[0]; 
int     nrowZ =  dimP[2], ncolZ =  dimP[3],  nzmaxZ = nzmaxP[1];
int     nrowA =  dimP[4], ncolA =  dimP[5],  nzmaxA = nzmaxP[2];
int     nrowLX =  dimP[6], ncolLX =  dimP[7],  nzmaxLX = nzmaxP[3];
int     nL = ncolLX/nrowLX;  // number of structural parameters

int 	dimAS =  ncolX+ncolZ;
int     nMH = 0;
int     nalpha = 0; 
 
bool    pr = prP[0];
bool    pl = prP[1];
bool    cp = nordinal>0;
bool    Aexists = FALSE;
bool    missing = FALSE;

double detLambda[2];   // holds old and proposed fabs(det(Lambda))
int sign_detLambda[2];    // holds sign of fabs(det(Lambda))
int lambda_old = 0;    // indexes freed Lambda
int lambda_new = 1;    // indexes freed Lambda
int *diagLambdaL = new int[ny]; // indexes diagonal elements of L and U from LU decomposition of Lambda
int *diagLambdaU = new int[ny];


double log_alphaL;     // MH ratio for lambda



        for(i=0; i<nG; i++){if(AtermP[i]==1){Aexists=TRUE;}}

        cnt2=0;
        for(k=nG; k<nGR; k++){
          for(i=0; i<nlGR[k]; i++){
            if(mvtype[i+cnt2]!=2){
              missing=TRUE;
              if(mvtype[i+cnt2]<0){
                nMH++; 
              }
            }
          }
          cnt2+=nlGR[k];
        }

        for(k=0; k<nGR; k++){           
           nalpha += PXtermsP[k]*GRdim[k];
        }

double  densityl1,
        densityl2,
        dbar = 0.0,
        mdbar = 0.0,
        Eaccl= 0.0,
        alpha_star = (-double(GRdimP[nG])/(1.0-2.75*double(GRdimP[nG])))-0.133,  // optimal acceptance ratio
        rACCEPT = 0.9,
        qACCEPT = 2.0,
        *t = new double[nGR*2],
        *sd = new double[nGR*2],
        *wn = new double[nGR*2],
        *zn = new double[nGR*2],
        *ldet = new double[nR+nG],
        interval,
        remainder,
        dm = 1.0/PedDimP[1];  // 0.5 for pedigrees 1 for phylogenies for taking averaging of potential bv.

	
//double inf = std::numeric_limits<double>::max();

        int *cumsum_ncutpoints = new int[nordinal+1];
            cumsum_ncutpoints[0] = 0;
        int ncutpoints_store = 0;

        if(cp){
          for(i=0; i<nordinal; i++){          
            cumsum_ncutpoints[i+1] = ncutpointsP[i]+cumsum_ncutpoints[i];  
            if(ncutpointsP[i]>3){
              ncutpoints_store += ncutpointsP[i]-3;
            }
          }
        } 

        double  *oldcutpoints = new double[cumsum_ncutpoints[nordinal]],
                *newcutpoints = new double[cumsum_ncutpoints[nordinal]],
                *sdcp = new double[nordinal],
                *wncp = new double[nordinal],
                *zncp = new double[nordinal],
                *accp = new double[nordinal],
                *cutpointMHR = new double[nordinal];

        if(cp){
          cnt = 0;
          for(i=0; i<nordinal; i++){
            sdcp[i] =1.0;
            wncp[i] =1.0;
            zncp[i] =1.0;
            accp[i] =0.0;
            for(j=0; j<ncutpointsP[i]; j++){
              oldcutpoints[cnt] = stcutpointsP[cnt];
              newcutpoints[cnt] = stcutpointsP[cnt];
              cnt ++;
            }
          }
        }

        if(ncutpoints_store==0){cp=FALSE;}  // even if cutpoints are present they do not need to be updated if there are only 2 categories

        for(k=0; k<(nGR*2); k++){
          t[k]=0.0;
          sd[k]=1.0;
          wn[k]=0.0;
          zn[k]=0.0;
        }

cs      *X, *Z, *W,  *Wt, *KRinv, *WtmKRinv, *WtmKRinvtmp, *M, *Omega, *MME, *zstar, *astar, *astar_tmp, *location, *location_tmp, *linky, *mulinky,  *pred, *mupred, *dev, *linki, *linki_tmp, *predi, *A, *bv, *bv_tmp, *bvA, *bvAbv, *tbv, *pvB, *pmuB, *Brv, *Xalpha, *tXalpha, *Alphainv, *muAlpha, *XtmKRinv, *XtmKRinvtmp, *alphaM, *alphaMME, *alphaastar, *alphapred, *alphazstar, *alphaastar_tmp, *alphalocation, *alphalocation_tmp, *Worig, *LambdaX, *pvL, *pmuL, *Lrv, *I, *linky_orig, *Y, *tY, *ILY, *w, *tYKrinv, *tYKrinvY, *tYKrinvw, *lambda_dev, *tl, *tlV, *tlVl;

csn	*L, *pvBL, *alphaL, *AlphainvL, *pvLL, *LambdaLU, *tYKrinvYL;
css     *S, *pvBS, *alphaS, *AlphainvS, *pvLS, *LambdaS, *tYKrinvYS;

cs*     *Ginv = new cs*[nGR];
cs*     *muG = new cs*[nGR];	
cs*     *propC = new cs*[nGR*2];
cs*     *propCinv = new cs*[nGR*2];
cs*     *muC = new cs*[nGR*2];
cs*     *G = new cs*[nGR];
cs*     *pG = new cs*[nGR];
cs*     *CM = new cs*[nGR];	
cs*     *Gtmp = new cs*[nGR];
cs*     *Grv = new cs*[nGR];
css*    *GinvS = new css*[nGR];
csn*    *GinvL = new csn*[nGR];
css*    *propCinvS = new css*[nGR*2];
csn*    *propCinvL = new csn*[nGR*2];
cs*     *KGinv = new cs*[nGR];
cs*     *lambda = new cs*[2];
cs*     *lambdaI = new cs*[2];
cs*     *Lambda = new cs*[2];
cs*     *Lambda_tmp = new cs*[2];

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

       pvB = cs_spalloc(ncolX, ncolX, pow(ncolX, 2.0-BvpP[ncolX*ncolX]), true, false);
       pmuB = cs_spalloc(ncolX, 1, ncolX, true, false);
       Brv = cs_spalloc(ncolX, 1, ncolX, true, false);

       cnt = 0;
       for (i = 0 ; i < ncolX ; i++){
         pmuB->i[i] = i;
         Brv->i[i] = i;
         pmuB->x[i] = BmupP[i];
       }
       pmuB->p[0] = 0; 
       pmuB->p[1] = ncolX;
       Brv->p[0] = 0; 
       Brv->p[1] = ncolX;

       if(BvpP[ncolX*ncolX]>0.5){   // prior for beta is diagonal
         for (i = 0 ; i < ncolX ; i++){
           pvB->p[i] = i; 
           pvB->i[i] = i;
           pvB->x[i] = BvpP[ncolX*i+i];
         }
         pvB->p[ncolX] = ncolX;
       }else{
         cnt = 0;
         for (i = 0 ; i < ncolX ; i++){
           pvB->p[i] = i*ncolX; 
           for (j = 0 ; j < ncolX ; j++){
             pvB->i[cnt] = j;
             pvB->x[cnt] = BvpP[cnt];
             cnt++;
           }
         }
         pvB->p[ncolX] = ncolX*ncolX;
       }

       pvBS = cs_schol(1, pvB);                    // Symbolic factorisation of B
       pvBL = cs_chol(pvB, pvBS);                  // cholesky factorisation of B^{-1} for forming N(0, B)


if(nL>0){

// read in parth analytic design matrix 

        LambdaX = cs_spalloc(nrowLX, ncolLX, nzmaxLX, true, false); 
                                                                   
        for (i = 0 ; i < nzmaxLX ; i++){
          LambdaX->i[i] = iLXP[i];
          LambdaX->x[i] = xLXP[i];
        }
        for (i = 0 ; i <= ncolLX ; i++){
          LambdaX->p[i] = pLXP[i];
        }                                                         

// lambda parameter vector and form kronecker(lamba, I)

        lambda[0] = cs_spalloc(nL, 1, nL, true, false); 
        lambda[1] = cs_spalloc(nL, 1, nL, true, false); 

         for (i = 0 ; i < nL ; i++){
          lambda[0]->i[i] = i;
          lambda[0]->x[i] = i*0.01+0.01;
        }
        lambda[0]->p[0] = 0;
        lambda[0]->p[1] = nL;
        
// read in prior covaraince matrix and mean vector for structural parameters.

       pvL = cs_spalloc(nL, nL, pow(nL, 2.0-LvpP[nL*nL]), true, false);
       pmuL = cs_spalloc(nL, 1, nL, true, false);
       Lrv = cs_spalloc(nL, 1, nL, true, false);

       cnt = 0;
       for (i = 0 ; i < nL ; i++){
         pmuL->i[i] = i;
         Lrv->i[i] = i;
         pmuL->x[i] = LmupP[i];
       }
       pmuL->p[0] = 0; 
       pmuL->p[1] = nL;
       Lrv->p[0] = 0; 
       Lrv->p[1] = nL;

       if(LvpP[nL*nL]>0.5){   // prior for structural parameters is diagonal
         for (i = 0 ; i < nL ; i++){
           pvL->p[i] = i; 
           pvL->i[i] = i;
           pvL->x[i] = LvpP[nL*i+i];
         }
         pvL->p[nL] = nL;
       }else{
         cnt = 0;
         for (i = 0 ; i < nL; i++){
           pvB->p[i] = i*nL; 
           for (j = 0 ; j < nL; j++){
             pvL->i[cnt] = j;
             pvL->x[cnt] = LvpP[cnt];
             cnt++;
           }
         }
         pvL->p[nL] = nL*nL;
       }

       pvLS = cs_schol(1, pvL);                    // Symbolic factorisation of L
       pvLL = cs_chol(pvL, pvLS);                  // cholesky factorisation of L^{-1} for forming N(0, L)
}
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

/**********************************/ 
/* read in G/R and prior matrices */
/**********************************/
 
        for (k = 0 ; k < nGR; k++){
          dimG = GRdim[k];
          Ginv[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          G[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          muG[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          pG[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          Grv[k] = cs_spalloc(1, dimG, dimG, true, false);
          cnt=0;
           for (i = 0 ; i < dimG; i++){
              Ginv[k]->p[i] = i*dimG;
              G[k]->p[i] = i*dimG;
              muG[k]->p[i] = i*dimG;
              pG[k]->p[i] = i*dimG;
              Grv[k]->p[i] = i;
              Grv[k]->i[i] = 0;
              for (j = 0 ; j < dimG; j++){
                 Ginv[k]->i[cnt] = j;
                 Ginv[k]->x[cnt] = GRinvP[cnt+tvc];
                 G[k]->i[cnt] = j;
                 G[k]->x[cnt] = GRinvP[cnt+tvc];
		 muG[k]->i[cnt] = j;
		 muG[k]->x[cnt] = 0.0;
                 pG[k]->i[cnt] = j;
                 pG[k]->x[cnt] = GRvpP[cnt+tvc];
                 cnt++;
              }
          }
          Ginv[k]->p[dimG] = dimG*dimG;
          G[k]->p[dimG] = dimG*dimG;
          muG[k]->p[dimG] = dimG*dimG;
          pG[k]->p[dimG] = dimG*dimG;
          Grv[k]->p[dimG] = dimG;
          ldet[k] = log(cs_invR(Ginv[k], G[k]));
          Gtmp[k] = cs_inv(Ginv[k]);  //delete this eventually  
          tvc += dimG*dimG; 
        }

/**************************************/	
/* Read in any condtional submatrices */
/**************************************/
	
	for (k = 0 ; k < nGR; k++){	
	   if(updateP[k]!=2){
	     CM[k] = cs_spalloc(1,1,1,true, false);
           }else{	
	     cnt = 0;	
	     CM[k] = cs_spalloc(dimG-splitP[k], dimG-splitP[k], (dimG-splitP[k])*(dimG-splitP[k]), true, false);
	     for (i = splitP[k] ; i < dimG; i++){
	        CM[k]->p[i-splitP[k]] = (i-splitP[k])*(dimG-splitP[k]);
		for (j = splitP[k] ; j < dimG; j++){
		   CM[k]->i[cnt] = j-splitP[k];
		   CM[k]->x[cnt] = G[k]->x[i*dimG+j];
		   cnt++;
		}
	      }
	      CM[k]->p[dimG-splitP[k]] = (dimG-splitP[k])*(dimG-splitP[k]);
	   }
	 }

/**********************************/ 
/* read in alpha  prior matrices */
/**********************************/
 
      if(nalpha>0){

        Alphainv = cs_spalloc(nalpha, nalpha, nzmaxAVpP[0], true, false); 
        muAlpha = cs_spalloc(nalpha, 1, nalpha, true, false); 
                                            
        for (i = 0 ; i < nzmaxAVpP[0]; i++){
          Alphainv->i[i] = iAVpP[i];
          Alphainv->x[i] = xAVpP[i];
        }
        for (i = 0 ; i < nalpha ; i++){
          Alphainv->p[i] = pAVpP[i];
          muAlpha->i[i] = i;
          muAlpha->x[i] = AmupP[i];
        }

        Alphainv->p[nalpha] = pAVpP[nalpha];
        muAlpha->p[0] = 0;
        muAlpha->p[1] = nalpha;
        AlphainvS = cs_schol(1, Alphainv);                    
        AlphainvL = cs_chol(Alphainv, AlphainvS);                

        Xalpha = cs_spalloc(nrowZ, nalpha, nrowZ*nalpha, true, false); 
                                                                   
        cnt=0;
        for (i = 0; i < nalpha; i++){
          Xalpha->p[i] = cnt;
          for (j = 0; j < nrowZ; j++){
            Xalpha->i[cnt] = j;
            Xalpha->x[cnt] = 0.0;
            cnt++;
          }
        }

        Xalpha->p[nalpha] = nrowZ*nalpha;
        tXalpha = cs_transpose(Xalpha, true);

        alphaastar = cs_spalloc(nalpha, 1, nalpha, true, false);
        alphalocation = cs_spalloc(nalpha, 1, nalpha, true, false);
        alphalocation_tmp = cs_spalloc(nalpha, 1, nalpha, true, false);

        for (i = 0 ; i < nalpha; i++){
           alphaastar->i[i] = i;    
           alphalocation->i[i] = i;  
           alphalocation->x[i] = 1.0;
           alphalocation_tmp->i[i] = i;    
        }
        alphaastar->p[0] = 0;
        alphaastar->p[1] = nalpha;
        alphalocation->p[0] = 0;
        alphalocation->p[1] = nalpha;
        alphalocation_tmp->p[0] = 0;
        alphalocation_tmp->p[1] = nalpha;

        alphazstar = cs_spalloc(ny, 1, ny, true, false);
        alphapred = cs_spalloc(ny, 1, ny, true, false);

        for (i = 0 ; i < ny; i++){   
           alphazstar->i[i] = i;
           alphapred->i[i] = i;
        }
        alphazstar->p[0] = 0; 
        alphazstar->p[1] = ny;
        alphapred->p[0] = 0; 
        alphapred->p[1] = ny;
      }

/*********************************/
/* Read in proposal distribution */
/*********************************/
        cnt2=0;
        for (k = nG ; k < nGR; k++){
          dimG = GRdim[k];
          propC[k] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          muC[k] = cs_spalloc(dimG, 1, dimG, true, false);
          propC[k+nGR] = cs_spalloc(dimG, dimG, dimG*dimG, true, false);
          muC[k+nGR] = cs_spalloc(dimG, 1, dimG, true, false);
          cnt=0;
          for (i = 0 ; i < dimG; i++){
            propC[k]->p[i] = i*dimG;
            propC[k+nGR]->p[i] = i*dimG;
            for (j = 0 ; j < dimG; j++){
              propC[k]->i[cnt] = j;
              propC[k]->x[cnt] = propCP[cnt2];
              propC[k+nGR]->i[cnt] = j;
              propC[k+nGR]->x[cnt] = propCP[cnt2];
              cnt++;
              cnt2++;
            }
          }
          propC[k]->p[dimG] = dimG*dimG;
          propC[k+nGR]->p[dimG] = dimG*dimG;
          muC[k]->p[0] = 0;
          muC[k]->p[1] = dimG;
          muC[k+nGR]->p[0] = 0;
          muC[k+nGR]->p[1] = dimG;
          for (i = 0 ; i < dimG; i++){
            muC[k]->i[i] = i;
            muC[k]->x[i] = 0.0;
            muC[k+nGR]->i[i] = i;
            muC[k+nGR]->x[i] = 0.0;
          }
          propCinv[k] = cs_inv(propC[k]);
          propCinvS[k] = cs_schol(1, propCinv[k]);   
          propCinvL[k] = cs_chol(propCinv[k], propCinvS[k]);       
          propCinv[k+nGR] = cs_inv(propC[k+nGR]);
          propCinvS[k+nGR] = cs_schol(1, propCinv[k+nGR]);   
          propCinvL[k+nGR] = cs_chol(propCinv[k+nGR], propCinvS[k+nGR]);
        }

// allocate vecotors for pseudo-random effects z* and [0, *a] 

        zstar = cs_spalloc(ny, 1, ny, true, false);
        linky = cs_spalloc(ny, 1, ny, true, false);
        pred = cs_spalloc(ny, 1, ny, true, false);
	mupred = cs_spalloc(ny, 1, ny, true, false);
	mulinky = cs_spalloc(ny, 1, ny, true, false);

        astar = cs_spalloc(dimAS, 1, dimAS, true, false);
        location = cs_spalloc(dimAS, 1, dimAS, true, false);
        location_tmp = cs_spalloc(dimAS, 1, dimAS, true, false);

        linki = cs_spalloc(dimG, 1, dimG, true, false);
        linki_tmp = cs_spalloc(dimG, 1, dimG, true, false);
        predi = cs_spalloc(dimG, 1, dimG, true, false);                 // This is dimG because this is the dimesnion of the R structures (all of which are the same)

        for (i = 0 ; i < ny; i++){   
           zstar->i[i] = i;
           linky->i[i] = i;
	   mupred->i[i] = i;
	   mupred->x[i] = 0.0;
	   pred->i[i] = i;
	   pred->x[i] = 0.0;
	   mulinky->i[i] = i;
	   mulinky->x[i] = 0.0;
           linky->x[i] = liabP[i];                         /* this needs to be changed for random regression */
        }

        zstar->p[0] = 0; 
        zstar->p[1] = ny;
        linky->p[0] = 0; 
        linky->p[1] = ny;
    	mupred->p[0] = 0; 
        mupred->p[1] = ny;
    	pred->p[0] = 0; 
        pred->p[1] = ny;
        mulinky->p[0] = 0; 
	mulinky->p[1] = ny;

        if(nL==0){       
          linky = cs_spalloc(ny, 1, ny, true, false);
          for (i = 0 ; i < ny; i++){   
            linky->i[i] = i;
            linky->x[i] = liabP[i];                        
          }
          linky->p[0] = 0; 
          linky->p[1] = ny; 
        }else{                                    /* for rec/sim models need to presevere non-transformed latent variable i.e. solve(Lambda, linky) - This =y (in gaussian case) */
          linky_orig = cs_spalloc(ny, 1, ny, true, false);
          for (i = 0 ; i < ny; i++){   
            linky_orig->i[i] = i;
            linky_orig->x[i] = liabP[i];                        
          }
          linky_orig->p[0] = 0; 
          linky_orig->p[1] = ny; 
        }

        for (i = 0 ; i < dimAS; i++){
           astar->i[i] = i;    
        }
        astar->p[0] = 0;
        astar->p[1] = dimAS;

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
        
        if(nalpha>0){

          Worig = cs_transpose(Wt, true);
          cnt = ncolX;
          cnt2 = 0;
          for (k = 0; k < nGR; k++){
            dimG = GRdim[k];
            for (i = 0; i < dimG; i++){
              if(PXtermsP[k]==1){    // parameter expanded 
                for (j = W->p[cnt]; j < W->p[cnt+nlGR[k]]; j++){ 
                  W->x[j] = Worig->x[j]*alphalocation->x[cnt2];     
                }
                cnt2 ++; 
              }
              cnt += nlGR[k];
            }
          } 
        }
  
       


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

        KRinv = cs_directsum(KGinv, nG, nGR);

// form WtmKRinv = W^{t}%*%KRinv  and t(t(WtmKRinv)) so its ordered correctly
        
        WtmKRinv = cs_multiply(Wt, KRinv); 
	WtmKRinvtmp = cs_transpose(WtmKRinv, TRUE); 
	cs_spfree(WtmKRinv);
	WtmKRinv = cs_transpose(WtmKRinvtmp, TRUE); 	
	cs_spfree(WtmKRinvtmp);

// form M = WtmKRinv%*%W  
   
        M = cs_multiply(WtmKRinv, W);

// form Omega = bdiag(0, KGinv) 

        Omega = cs_omega(KGinv, nG, pvB);

// form MME = M + Omega; mixed model equations 

        MME = cs_add(M, Omega, 1.0, 1.0); 
        S = cs_schol(1, MME);                            // Symbolic factorisation - only has to be done once

// If sim/rec model then:

       if(nL>0){
     
         lambdaI[0] = cs_kroneckerI(lambda[0], LambdaX->m);
         lambdaI[1] = cs_kroneckerI(lambda[0], LambdaX->m);

         Lambda_tmp[0] = cs_multiply(LambdaX, lambdaI[0]);

         I = cs_spalloc(ny, ny, ny, true, false); // ny by ny I matrix

         for (i = 0 ; i<ny ; i++){
           I->p[i] = i;                                                               
           I->i[i] = i;
           I->x[i] = 1.0;
         }
         I->p[ny] = ny;   

         Lambda[0] = cs_add(I, Lambda_tmp[0], 1.0, -1.0);

         LambdaS = cs_sqr(1, Lambda[0], true);                            // Symbolic LU factorisation of Lambda for determinant calculation

         LambdaLU = cs_lu(Lambda[0], LambdaS, DBL_EPSILON);      

         for (i = 0 ; i<ny ; i++){           
           for (j = LambdaLU->L->p[i] ; j<(LambdaLU->L->p[i+1]); j++){
              if(LambdaLU->L->i[j]==i){
                diagLambdaL[i]=j;
              }
           }
         }

         for (i = 0 ; i<ny ; i++){           
           for (j = LambdaLU->U->p[i] ; j<(LambdaLU->U->p[i+1]); j++){
              if(LambdaLU->U->i[j]==i){
                diagLambdaU[i]=j;
              }
           }
         }

         detLambda[0] = 0.0;
         sign_detLambda[0] = 1;

         for (i = 0 ; i<ny ; i++){ 
            if(LambdaLU->L->x[diagLambdaL[i]]<0.0){
              sign_detLambda[0] *= -1;
            }
            if(LambdaLU->U->x[diagLambdaU[i]]<0.0){
              sign_detLambda[0] *= -1;
            }
            detLambda[0] += log(fabs(LambdaLU->L->x[diagLambdaL[i]]))+log(fabs(LambdaLU->U->x[diagLambdaU[i]]));
         }

 
         linky = cs_multiply(Lambda[0], linky_orig);   

         cs_sortdv(linky);

         ILY = cs_spalloc(ny*nL, nL, ny*nL, true, false); // form kronecker(I,y) where I has dimension nL
       
         for (i = 0 ; i<nL ; i++){
           ILY->p[i] = ny*i;                                                               
           for (j = 0 ; j < ny ; j++){
             ILY->i[j+ny*i] = j+ny*i;
             ILY->x[j+ny*i] = linky_orig->x[j];
           }
         }
         ILY->p[nL] = nL*ny;   
         Y = cs_multiply(LambdaX, ILY);

         tY = cs_transpose(Y, true);

         tYKrinv = cs_multiply(tY, KRinv);
         tYKrinvY = cs_multiply(tYKrinv,Y);
         
         tYKrinvYS = cs_schol(1, tYKrinvY);
         tYKrinvYL = cs_chol(tYKrinvY, tYKrinvYS);
       }

        if(nalpha>0){

          XtmKRinv = cs_multiply(tXalpha, KRinv); 
	  XtmKRinvtmp = cs_transpose(XtmKRinv, TRUE); 
	  cs_spfree(XtmKRinv);
	  XtmKRinv = cs_transpose(XtmKRinvtmp, TRUE); 	// double transposes so there row ordered within columns!
	  cs_spfree(XtmKRinvtmp);

// form alphaM = t(alphaX)%*%solve(R)%*%alphaX%*%
   
          alphaM = cs_multiply(XtmKRinv, Xalpha);

// form MME = M + Omega; mixed model equations 

          alphaMME = cs_add(alphaM, Alphainv, 1.0, 1.0); 

          alphaS = cs_schol(1, alphaMME);                            // Symbolic factorisation - only has to be done once
        }
	
  	GetRNGstate();                                   // get seed for random number generation

/**********************************************************************************************************/
/***** MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC MCMC *****/
/**********************************************************************************************************/

        for (itt = 0; itt < (nitt+DICP[0]); itt++){
			
          if(itt>0){
            for (i = 0 ; i < nGR; i++){
              if(updateP[i]>0){
                cs_spfree(G[i]); 
              }
            }
            cs_spfree(astar_tmp);
            cs_spfree(pred);
            cs_spfree(dev);
  	    cs_nfree(L);

            if(nalpha>0){
              cs_nfree(alphaL);
              cs_spfree(alphaastar_tmp); 
            }
          }

          if(missing){	
            for (i = nG ; i < nGR; i++){
	      cs_nfree(propCinvL[i]);
	      cs_spfree(propCinv[i]);   
	      cs_nfree(propCinvL[i+nGR]);
	      cs_spfree(propCinv[i+nGR]);   
            }
          }
	  cs_spfree(WtmKRinv);     
	  cs_spfree(M);                               
	  cs_spfree(MME);  

          if(nalpha>0){
            cs_spfree(tXalpha);
            cs_spfree(XtmKRinv);
            cs_spfree(alphaM);
            cs_spfree(alphaMME);
          }

          if(nL>0){
  	     cs_nfree(LambdaLU);
             cs_nfree(tYKrinvYL);     
             cs_spfree(tYKrinv);
             cs_spfree(tYKrinvY);                
           //  cs_spfree(Y);     // for non-gaussian or missing data Y and tY need updating as linky_orig changes
           //  cs_spfree(tY); 
          }

/******************/
/* form equations */
/******************/
						
			for (k = 0 ; k < nGR; k++){                             
				if(updateP[k]>0){
					if(AtermP[k]==1){
						cs_kroneckerAupdate(Ginv[k],A,KGinv[k]);              //  form kronecker(G^{-1}, A^{-1}) structure
					}else{
						cs_kroneckerIupdate(Ginv[k],nlGR[k],KGinv[k]);        //  form G^{-1} structure
					}
				}  
			}
			
			cs_directsumupdate(KGinv, nG, nGR, KRinv);
			
			WtmKRinv = cs_multiply(Wt, KRinv);   
			
			M = cs_multiply(WtmKRinv, W);                          // form M = WtmKRinv%*%W  
			
			cs_omegaupdate(KGinv, nG, pvB, Omega);                 // update Omega = bdiag(0, KGinv) 
			
			MME = cs_add(M, Omega, 1.0, 1.0);                      // form MME = M + Omega; mixed model equations 


/*************************/
/*  sample vectors from  */
/*      their prior      */
/*************************/

/* beta */			
          cnt = 0;

          for(i=0; i<ncolX; i++){
              Brv->x[i] = rnorm(0.0,1.0);
          }                   
          cs_ltsolve(pvBL->L, Brv->x);
          for(i=0; i<ncolX; i++){
              astar->x[i]  = Brv->x[i]+pmuB->x[i];
          }                          // sample fixed effects from their prior
          cnt = ncolX;
			
	  for(k=0; k<nG; k++){
            dimG = GRdim[k];
            if(AtermP[k]==1){          

/* pedigree and phylogeny effects - to become non iid effects generally  */

		
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

/* blocked random effects */
				
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

/* residuals */
			
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

/************************/
/* sample pseudo vector */
/************************/

         cs_gaxpy(W, astar->x, zstar->x);

         for (i = 0 ; i < ny ; i++){
            zstar->x[i] *= -1.0;
            zstar->x[i] += linky->x[i];       // form y - z*
         }

         astar_tmp = cs_multiply(WtmKRinv, zstar);              // WtmKRinv(y - z*)

         for (i = 0 ; i < astar_tmp->nzmax ; i++){
            location->x[astar_tmp->i[i]] = astar_tmp->x[i];
         }

/*************/
/* solve MME */
/*************/

         L = cs_chol(MME, S); 

         if(L==NULL){
	   PutRNGstate();
           error("Mixed model equations singular: use a (stronger) prior\n");
         }

         for (i = 0 ; i < dimAS; i++){
            location_tmp->x[i] = 0.0;
         }

         cs_ipvec (S->pinv, location->x, location_tmp->x, MME->n);	 // x = P*b 
         cs_lsolve(L->L,location_tmp->x);                                // x = L\x 
         cs_ltsolve (L->L, location_tmp->x);		                 // x = L'\x 
         cs_pvec (S->pinv, location_tmp->x, location->x, MME->n);        // b = P'*x 

          for (i = 0 ; i < dimAS ; i++){
            location->x[i] += astar->x[i];
         }

/***********************/
/* sample VCV matrices */
/***********************/

         pred = cs_multiply(W, location);
         cs_sortdv(pred);
         dev = cs_add(linky, pred, 1.0, -1.0);    
         cnt2 = ncolX;

	 for(i=0; i< nG; i++){        
           if(updateP[i]>0){            
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
               cs_spfree(tbv);               
               cs_spfree(bvA);               
	       cs_spfree(bvAbv);
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
	     switch(updateP[i]){					   
		case 1:
 		  cs_invR(Gtmp[i], Ginv[i]);
		  G[i] = cs_rinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], GinvS[i]);	 
		break;
					   
		case 2:    
		  G[i] = cs_rCinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], splitP[i], CM[i]); 
		break;	   
					   
		case 3:  				 
		  G[i] = cs_rR(Gtmp[i], double(nlGR[i]), GRnpP[i], GinvS[i], Ginv[i], ldet[i]);
		break;
	     }		
	     ldet[i] = log(cs_invR(G[i], Ginv[i]));  
	     cs_nfree(GinvL[i]);
	     GinvL[i] = cs_chol(Ginv[i], GinvS[i]);                 // cholesky factorisation of G^{-1} for Gibbs sampling PX working parameters. 
           }
         }

/**********************/
/* Sample R Structure */
/**********************/
         cnt2=0;
			

         for(i=nG; i<nGR; i++){    
            if(updateP[i]>0){ 
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
	          if(diagR>0){
	            cnt=0;  
	            for(j=0; j<dimG; j++){
		      for(k=0; k<dimG; k++){
		        if(j==k){
                          if(diagR==2 && cnt!=0){
                            Gtmp[i]->x[0] += Gtmp[i]->x[cnt] - pG[i]->x[cnt];
                          }
                        }else{
			  Gtmp[i]->x[cnt] = 0.0;
                        }
		        cnt++;
		      }
	            }	
	          }
                  switch(updateP[i]){

	            case 1: 
		      cs_invR(Gtmp[i], Ginv[i]);
		      G[i] = cs_rinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], GinvS[i]);
		    break;
						
	            case 2: 
		      cs_invR(Gtmp[i], Ginv[i]);
		      G[i] = cs_rCinvwishart(Ginv[i], double(nlGR[i])+GRnpP[i], splitP[i], CM[i]);
		    break;	   
			
	            case 3:  
		      cs_cov2cor(Gtmp[i]);
		      G[i] = cs_rR(Gtmp[i], double(nlGR[i]), GRnpP[i], GinvS[i], Ginv[i], ldet[i]);  
	            break;
	          }					
	          if(diagR>0){
	            cnt=0;  
	            for(j=0; j<dimG; j++){
		      for(k=0; k<dimG; k++){
		        if(j==k){
                          if(diagR==2){
                            G[i]->x[cnt] = G[i]->x[0];
                          }                          
                        }else{
		          G[i]->x[cnt] = 0.0;
                        }
			cnt++;
		      }
		    }	
		  }
		  ldet[i] = log(cs_invR(G[i], Ginv[i]));
	          cs_nfree(GinvL[i]);	
	          GinvL[i] = cs_chol(Ginv[i], GinvS[i]);                 // cholesky factorisation of R^{-1} for Gibbs sampling fully missing data. 
		}
	    }

/**********************/
/* calculate deviance */   
/**********************/
     dbar =0.0;

     if(itt>=burnin && DICP[0]==1){
       if(itt==nitt){
          for(k=nG; k<nGR; k++){
             dimG = GRdim[k];
             for(i=0; i<(dimG*dimG); i++){          
               G[k]->x[i] = muG[k]->x[i];
             }
             ldet[k] = log(cs_invR(G[k], Ginv[k]));  
          }
          for(i=0; i<ny; i++){          
            linky->x[i] = mulinky->x[i];
            pred->x[i] = mupred->x[i];
          }
       }         
       cnt2=0;
       for(k=nG; k<nGR; k++){          // Iterate through R-structures
         dimG = GRdim[k];
         for(j=0; j<nlGR[k]; j++){     // Iterate through levels
           nkeep=0;
           ncond=0;
           for(i=0; i<dimG; i++){       // Iterate through fixed levels
             record=cnt2+nlGR[k]*i+j;
             linki->x[i] = linky->x[record];              
             predi->x[i] =  pred->x[record];
             if(familyP[record]==1 && observedP[record]==1){
               keep[nkeep] = i;
               nkeep ++;
             }else{
               cond[ncond] = i;
               ncond ++;
             }
           }
           if(nkeep>0){      // some gaussian observed traits
             if(ncond>0){   // some non-gaussian or non-observed traits
               dbar += cs_dcmvnorm(linki, predi, ldet[k], Ginv[k], G[k], keep, nkeep, cond, ncond);    // some gaussian observed
             }else{
               dbar += cs_dmvnorm(linki, predi, ldet[k], Ginv[k]);                                     // all gaussian observed
             }
           }
         }
         cnt2+=nlGR[k]*dimG;
       }
     }

/********************/
/* update cutpoints */
/********************/

     if(cp){  
       for(i=0; i<nordinal; i++){ 
         for(j=2; j<(ncutpointsP[i]-1); j++){ 
            newcutpoints[cumsum_ncutpoints[i]+j] = rtnorm(oldcutpoints[cumsum_ncutpoints[i]+j], sdcp[i], newcutpoints[cumsum_ncutpoints[i]+j-1], oldcutpoints[cumsum_ncutpoints[i]+j+1]);
         } 
       }
       cnt2=0;
       rterm=0;

       for(k=nG; k<nGR; k++){      // Iterate through R-structures
         for(i=0; i<dimG; i++){    // Iterate through first indiviual to find any ordinal variables
           dimG = GRdim[k];
           record=cnt2+nlGR[k]*i;
           if(familyP[record]==14){
             nthordinal = mfacP[rterm+i];
             cutpointMHR[nthordinal] = dcutpoints(linky, yP, observedP, record,record+nlGR[k], oldcutpoints, newcutpoints, cumsum_ncutpoints[nthordinal], ncutpointsP[i], sdcp[nthordinal]);
             wncp[nthordinal] *= rACCEPT;
             zncp[nthordinal] *= rACCEPT;
             wncp[nthordinal] ++;

             if(cutpointMHR[nthordinal]>log(runif(0.0,1.0))){
               zncp[nthordinal] ++;
               accp[nthordinal] ++;
               for(j=2; j<(ncutpointsP[nthordinal]-1); j++){ 
                 oldcutpoints[cumsum_ncutpoints[nthordinal]+j] = newcutpoints[cumsum_ncutpoints[nthordinal]+j];
               }
             } 
             if(itt<burnin){          
               sdcp[nthordinal] *= pow(qACCEPT, ((zncp[nthordinal]/wncp[nthordinal])-0.44));
             }
	   }	 
         }
         rterm += dimG;
         cnt2+=nlGR[k]*dimG;
       }
     }

/********************************************************/
/* update recusrive-simultaneous structural parameters  */   
/********************************************************/


    if(nL>0){    
                
      for (i = 0 ; i<ny ; i++){            // updates missing values for records i when Lambda[i,i] = 1 and Lambda[i,\i] = 0
         if(observedP[i]==0){
           linky_orig->x[i] = linky->x[i];
         }
      }

              
      w = cs_add(linky_orig, pred, 1.0, -1.0);
      tYKrinv = cs_multiply(tY, KRinv);
      tYKrinvY = cs_multiply(tYKrinv,Y);
      tYKrinvw = cs_multiply(tYKrinv,w);

      for(i=0; i<nL; i++){
        Lrv->x[i] = rnorm(0.0,1.0);
      }


      tYKrinvYL = cs_chol(tYKrinvY, tYKrinvYS); 

      cs_ltsolve(tYKrinvYL->L, Lrv->x);       // ~ N(0, V_lambda)

      cs_ipvec (tYKrinvYS->pinv, tYKrinvw->x, w->x, nL);	 
      cs_lsolve(tYKrinvYL->L,w->x);                                
      cs_ltsolve (tYKrinvYL->L, w->x);		                 
      cs_pvec (tYKrinvYS->pinv, w->x, tYKrinvw->x, nL);       

      lambda[lambda_new] = cs_add(Lrv, lambda[lambda_old], 1.0, 1.0);

      cs_kroneckerIupdate(lambda[lambda_new], LambdaX->m, lambdaI[lambda_new]);    
      Lambda_tmp[lambda_new] = cs_multiply(LambdaX, lambdaI[lambda_new]);    
      
      Lambda[lambda_new] = cs_add(I, Lambda_tmp[lambda_new], 1.0, -1.0);

      LambdaLU = cs_lu(Lambda[lambda_new], LambdaS, DBL_EPSILON);     

      detLambda[lambda_new] = 0.0;
      sign_detLambda[lambda_new] = 1;

      for (i = 0 ; i<ny ; i++){ 
         if(LambdaLU->L->x[diagLambdaL[i]]<0.0){
           sign_detLambda[lambda_new] *= -1;
         }
         if(LambdaLU->U->x[diagLambdaU[i]]<0.0){
           sign_detLambda[lambda_new] *= -1;
         }

         detLambda[lambda_new] += log(fabs(LambdaLU->L->x[diagLambdaL[i]]))+log(fabs(LambdaLU->U->x[diagLambdaU[i]]));
      }

        log_alphaL = detLambda[lambda_new] - detLambda[lambda_old];

        lambda_dev = cs_add(lambda[lambda_new], tYKrinvw, 1.0, -1.0);

        tl = cs_transpose(lambda_dev, true);
        tlV = cs_multiply(tl, tYKrinvY);
        tlVl = cs_multiply(tlV, lambda_dev);

        log_alphaL -= 0.5*tlVl->x[0];

        cs_spfree(lambda_dev);
        cs_spfree(tl);
        cs_spfree(tlV);
        cs_spfree(tlVl);

        lambda_dev = cs_add(lambda[lambda_old], tYKrinvw, 1.0, -1.0);

        tl = cs_transpose(lambda_dev, true);
        tlV = cs_multiply(tl, tYKrinvY);
        tlVl = cs_multiply(tlV, lambda_dev);

        log_alphaL += 0.5*tlVl->x[0];

        cs_spfree(lambda_dev);
        cs_spfree(tl);
        cs_spfree(tlV);
        cs_spfree(tlVl);

      if(sign_detLambda[lambda_old]==sign_detLambda[lambda_new] && log_alphaL>log(runif(0.0,1.0))){

        cs_spfree(linky);  
        linky = cs_multiply(Lambda[lambda_new], linky_orig);
        cs_sortdv(linky);                                    // re-sort latent variable (data in Gaussian case)
 
      // only needs upadting if latent variables exist - need to free memory for Y and tY when implemented

      /*
        for (i = 0 ; i<nL ; i++){
          for (j = 0 ; j < ny ; j++){
            ILY->x[j+ny*i] = linky_orig->x[j];
          }
        }
        Y = cs_multiply(LambdaX, ILY);
        tY = cs_transpose(Y, true);
      */      

        detLambda[0] = detLambda[lambda_new];
        detLambda[1] = detLambda[lambda_new];
        sign_detLambda[0] = sign_detLambda[lambda_new];
        sign_detLambda[1] = sign_detLambda[lambda_new];
        lambda_old =  lambda_new;
        lambda_new -= 1;
        lambda_new = abs(lambda_new);

      }
      cs_spfree(Lambda_tmp[lambda_new]);                                     
      cs_spfree(Lambda[lambda_new]);   
      cs_spfree(lambda[lambda_new]);   
      cs_spfree(w);
      cs_spfree(tYKrinvw);
    }

/***********************/
/* sample liabilities  */   
/***********************/
			
     if(missing){

       cnt2=0;
       cnt=0;
       rterm=0;  // indexes the individual R-level so with the terms in the R-structure (1 1-dimensional, and 1 2 dimensional) {1} + {2, 3}

       for(k=nG; k<nGR; k++){      // Iterate through R-structures

         dimG = GRdim[k];
 	 propCinv[k] = cs_inv(propC[k]);
	 propCinvL[k] = cs_chol(propCinv[k], propCinvS[k]); 
	 propCinv[k+nGR] = cs_inv(propC[k+nGR]);
	 propCinvL[k+nGR] = cs_chol(propCinv[k+nGR], propCinvS[k+nGR]);
   
         for(j=0; j<nlGR[k]; j++){     // Iterate through levels

           densityl1 = 0.0;
           densityl2 = 0.0;
 
           nthmnl = 0;          //variables for the multinomial
           mndenom1 = 1.0;
           mndenom2 = 1.0;

           p = k+proposal[cnt+j]*nGR;  // indexes proposal distributions

           for(i=0; i<dimG; i++){
             linki_tmp->x[i] = rnorm(0.0,1.0);
           }

           if(mvtype[cnt+j]==1){         // can be Gibbsed  
             cs_ltsolve(GinvL[k]->L, linki_tmp->x);  
             for(i=0; i<dimG; i++){
               linky->x[nlGR[k]*i+j+cnt2] = linki_tmp->x[i]+pred->x[nlGR[k]*i+j+cnt2];
             }
           }
			 
           if(mvtype[cnt+j]<0){         // has to be MHed
             cs_ltsolve(propCinvL[p]->L, linki_tmp->x);
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

                   case 3:  /* Nominal Multinomial Logit */

                     if(mvtype[cnt+j]==0){   // univraiate binary models can be slice sampled
                       if(yP[record]>0.5){
                         u = linky->x[record]-log1p(exp(linky->x[record]));  // needed for the deviance calculation
                         densityl1 += u;
                         u -= rexp(1.0);
                         linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]),u-log1p(-exp(u)), 1e+35);
                       }else{
                         u = -log1p(exp(linky->x[record])); // needed for the deviance calculation
                         densityl1 += u;
                         u -= rexp(1.0);
                         linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]), -1e+35, log1p(-exp(u))-u);				
                       }

                     }else{

                       mndenom1 += exp(linki->x[i]);
                       mndenom2 += exp(linki_tmp->x[i]);
                       densityl1 += yP[record]*linki->x[i];
                       densityl2 += yP[record]*linki_tmp->x[i];

                       if(mfacP[rterm+i]==nthmnl){ 
                         densityl1 -= y2P[record]*log(mndenom1);
                         densityl2 -= y2P[record]*log(mndenom2);
                         nthmnl = 0;
                         mndenom1 = 1.0;
                         mndenom2 = 1.0;
                       }else{
                         nthmnl++;
                       }
                     }
                   break;
     
                   case 4: /* Weibull */
                     densityl1 += dweibull(yP[record], 1.0, exp(-linki->x[i]), true);
                     densityl2 += dweibull(yP[record], 1.0, exp(-linki_tmp->x[i]), true);
                   break;
    
                   case 5: /* Exponential */
                     densityl1 += dexp(yP[record], exp(-linki->x[i]), true);
                     densityl2 += dexp(yP[record], exp(-linki_tmp->x[i]), true);
                   break;
    
                   case 6: /* Censored Gaussian */

                     if(linki_tmp->x[i]<yP[record] || linki_tmp->x[i]>y2P[record]){
                        interval = (linki_tmp->x[i]-y2P[record])/(y2P[record]-yP[record]);
                        remainder = (linki_tmp->x[i]-y2P[record])-double(int(interval))*(y2P[record]-yP[record]);
                       if(linki_tmp->x[i]<yP[record]){
                         if(int(interval)%2==1){
                           linki_tmp->x[i] = yP[record] - remainder;
                         }else{
                           linki_tmp->x[i] = y2P[record] + remainder;
                        }
                      }else{
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
                     densityl1 += log(pweibull(y2P[record], 1.0, exp(-linki->x[i]), true, false)-pweibull(yP[record],1.0, exp(-linki->x[i]), true, false));
                     densityl2 += log(pweibull(y2P[record], 1.0, exp(-linki_tmp->x[i]), true, false)-pweibull(yP[record],1.0, exp(-linki_tmp->x[i]), true, false));
                   break;                                            

                   case 9: /* Censored Exponential */
                     densityl1 += log(pexp(y2P[record], exp(-linki->x[i]), true, false)-pexp(yP[record], exp(-linki->x[i]), true, false));
                     densityl2 += log(pexp(y2P[record], exp(-linki_tmp->x[i]), true, false)-pexp(yP[record], exp(-linki_tmp->x[i]), true, false));
                   break;  
						 
		   case 10: /* Zero-inflated Gaussian */
		   break;
						 
		   case 11: /* Zero-inflated Poisson */

                     if(mfacP[rterm+i]==0){
                       mndenom1 = dpois(yP[record], exp(linki->x[i]), true);  
                       mndenom2 = dpois(yP[record], exp(linki_tmp->x[i]), true);  
                     }else{
			mndenom1 += log1p(-exp(linki->x[i])/(1.0+exp(linki->x[i])));
			mndenom2 += log1p(-exp(linki_tmp->x[i])/(1.0+exp(linki_tmp->x[i])));
                        if(yP[record]>0.5){
			  mndenom1 = log(exp(mndenom1)+exp(linki->x[i])/(1.0+exp(linki->x[i])));
			  mndenom2 = log(exp(mndenom2)+exp(linki_tmp->x[i])/(1.0+exp(linki_tmp->x[i])));  
			}
                        densityl1 += mndenom1;
                        densityl2 += mndenom2;
                        mndenom1 = 1.0;
                        mndenom2 = 1.0;
                     }

		   break;
						 
		   case 12: /* Zero-inflated Weibull */ 
                   break;                                            
						 
                   case 13: /* Zero-inflated Exponential */
		   break;  

                   case 14: /* Ordered Mulinomial Probit */

                     if(mvtype[cnt+j]==0){   // univraiate binary models can be slice sampled

                       
                       if(yP[record]>1.5){
	                 u = pnorm(linky->x[record], 0.0, 1.0, true, true);
                         densityl1 += u;    // needed for the deviance calculation   
                         u -= rexp(1.0);
                         linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]), qnorm(u, 0.0, 1.0, true, true), 1e+35);
                       }else{
                         u = pnorm(linky->x[record], 0.0, 1.0, false, true);
                         densityl1 += u;   // needed for the deviance calculation
                         u -= rexp(1.0);
                         linky->x[record] = rtnorm(pred->x[record], sqrt(G[k]->x[0]), -1e+35, qnorm(u, 0.0, 1.0, false, true));		
                       }

                     }else{
                       nthordinal = mfacP[rterm+i];

                       if(int(yP[record])==1 || int(yP[record])==(ncutpointsP[nthordinal]-1)){
                         if(int(yP[record])==1){
                           densityl1 += pnorm(linki->x[i], 0.0, 1.0, false,true);
                           densityl2 += pnorm(linki_tmp->x[i], 0.0, 1.0, false,true);
                         }else{
                           densityl1 += pnorm(oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]]-linki->x[i], 0.0, 1.0, false,true);
                           densityl2 += pnorm(oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]]-linki_tmp->x[i], 0.0, 1.0, false,true);
                         }
                       }else{
                         densityl1 += log(pnorm(oldcutpoints[int(yP[record])+cumsum_ncutpoints[nthordinal]]-linki->x[i], 0.0, 1.0, true,false)-pnorm(oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]]-linki->x[i], 0.0, 1.0, true,false));
                         densityl2 += log(pnorm(oldcutpoints[int(yP[record])+cumsum_ncutpoints[nthordinal]]-linki_tmp->x[i], 0.0, 1.0, true,false)-pnorm(oldcutpoints[int(yP[record])-1+cumsum_ncutpoints[nthordinal]]-linki_tmp->x[i], 0.0, 1.0, true,false));
                       }
                     }
                   break;

                   case 15:  /* Hurdle Poisson */

                    if(mfacP[rterm+i]==0){
                       if(yP[record]>0.5){
                         mndenom1 = dpois(yP[record], exp(linki->x[i]), true)-log1p(-exp(-exp(linki->x[i])));
                         mndenom2 = dpois(yP[record], exp(linki_tmp->x[i]), true)-log1p(-exp(-exp(linki_tmp->x[i]))); 
                       } 
                     }else{
	               if(yP[record]>0.5){
                         mndenom1 = (linki->x[i] - log1p(exp(linki->x[i])));
			 mndenom2 = (linki_tmp->x[i] - log1p(exp(linki_tmp->x[i]))); 		       		 
		       }else{
                         mndenom1 += -linki->x[i] - log1p(exp(-linki->x[i]));
			 mndenom2 += -linki_tmp->x[i] - log1p(exp(-linki_tmp->x[i]));     
                       }
                        densityl1 += mndenom1;
                        densityl2 += mndenom2;
                        mndenom1 = 1.0;
                        mndenom2 = 1.0;
                     }
                   break;

                   case 16:  /* Zero-truncated Poisson */

                     densityl1 += dpois(yP[record], exp(linki->x[i]), true)-log1p(-exp(-exp(linki->x[i])));
                     densityl2 += dpois(yP[record], exp(linki_tmp->x[i]), true)-log1p(-exp(-exp(linki_tmp->x[i])));

                   break;

                   case 17:  /* Geometric */

                     densityl1 += linki->x[i]-log1p(exp(linki->x[i]))-yP[record]*(linki->x[i]+log1p(exp(-linki->x[i]))); 
                     densityl2 += linki_tmp->x[i]-log1p(exp(linki_tmp->x[i]))-yP[record]*(linki_tmp->x[i]+log1p(exp(-linki_tmp->x[i])));  

                   break;

                   case 18:  /* Zero-altered Poisson */

                    if(mfacP[rterm+i]==0){
                       if(yP[record]>0.5){
                         mndenom1 = dpois(yP[record], exp(linki->x[i]), true)-log1p(-exp(-exp(linki->x[i])));
                         mndenom2 = dpois(yP[record], exp(linki_tmp->x[i]), true)-log1p(-exp(-exp(linki_tmp->x[i]))); 
                       } 
                     }else{
	               if(yP[record]>0.5){
                         mndenom1 = log1p(-pexp(exp(linki->x[i]), 1.0,true, false));
        		   mndenom2 = log1p(-pexp(exp(linki_tmp->x[i]),1.0, true, false));	
		       		 	       		 
		       }else{
                         mndenom1 += pexp(exp(linki->x[i]),1.0,true, true);
   			 mndenom2 += pexp(exp(linki_tmp->x[i]),1.0,true, true);     
                       }
                        densityl1 += mndenom1;
                        densityl2 += mndenom2;
                        mndenom1 = 1.0;
                        mndenom2 = 1.0;
                     }

                   break;

                   case 19:  /* Zero-inflated Binomial */

                     if(mfacP[rterm+i]==0){  // non-zero bit
                       mndenom1 = dbinom(yP[record], y2P[record], exp(linki->x[i])/(1.0+exp(linki->x[i])), true);  
                       mndenom2 = dbinom(yP[record], y2P[record], exp(linki_tmp->x[i])/(1.0+exp(linki_tmp->x[i])), true);  
                     }else{
			mndenom1 += log1p(-exp(linki->x[i])/(1.0+exp(linki->x[i])));
			mndenom2 += log1p(-exp(linki_tmp->x[i])/(1.0+exp(linki_tmp->x[i])));
                        if(yP[record]>0.5){
			  mndenom1 = log(exp(mndenom1)+exp(linki->x[i])/(1.0+exp(linki->x[i])));
			  mndenom2 = log(exp(mndenom2)+exp(linki_tmp->x[i])/(1.0+exp(linki_tmp->x[i])));  
			}
                        densityl1 += mndenom1;
                        densityl2 += mndenom2;
                        mndenom1 = 1.0;
                        mndenom2 = 1.0;
                     }
                   break;

                 }
               }
             }
           }
           dbar += densityl1;

           if(mvtype[cnt+j]<0){

             densityl1 += cs_dmvnorm(linki, predi, ldet[k], Ginv[k]);
             densityl2 += cs_dmvnorm(linki_tmp, predi, ldet[k], Ginv[k]);
             zn[p] *= rACCEPT; 
             wn[p] *= rACCEPT;
             wn[p] ++;

             if((densityl2-densityl1)>log(runif(0.0,1.0))){
               for(i=0; i<dimG; i++){
                 linky->x[nlGR[k]*i+j+cnt2] = linki_tmp->x[i];
               }
               zn[p]++;
               Eaccl++;  
             }
/***************/
/* Adaptive MH */
/***************/

             if(itt<burnin && AMtuneP[k]==1){
               t[p] ++;
               for(i=0; i<dimG; i++){
                 for(l=0; l<dimG; l++){
                  propC[p]->x[i*dimG+l] *= (t[p]-1.0)/t[p];
                  propC[p]->x[i*dimG+l] += muC[p]->x[i]*muC[p]->x[l];
                 }
               }
	       for(i=0; i<dimG; i++){
		  muC[p]->x[i] *= (t[p]-1.0);
		  muC[p]->x[i] += linky->x[nlGR[k]*i+j+cnt2];
		  muC[p]->x[i] /= t[p];
	       }				 
               for(i=0; i<dimG; i++){
                 for(l=0; l<dimG; l++){
                   propC[p]->x[i*dimG+l] -= ((t[p]+1.0)/t[p])*muC[p]->x[i]*muC[p]->x[l];
                   propC[p]->x[i*dimG+l] += linky->x[nlGR[k]*i+j+cnt2]*linky->x[nlGR[k]*l+j+cnt2]/t[p];
                 }
                 propC[p]->x[i*(dimG+1)] += 0.001/t[p];
               }
             }
           }
         }
         if(itt<burnin && AMtuneP[k]==1){
           for(i=0; i<dimG; i++){
             for(l=0; l<dimG; l++){
	       propC[k]->x[i*dimG+l] += (muC[k]->x[i]*muC[k]->x[l])/(t[k]);
               propC[k]->x[i*dimG+l] *= sd[k];
               if(t[k+nGR]>0.0){
	         propC[k+nGR]->x[i*dimG+l] += (muC[k+nGR]->x[i]*muC[k+nGR]->x[l])/(t[k+nGR]);
                 propC[k+nGR]->x[i*dimG+l] *= sd[k+nGR];
               }
             }
           }
	   if(wn[k]>0.0){
             sd[k] = pow(qACCEPT, ((zn[k]/wn[k])-alpha_star));
	   }	
	   if(wn[k+nGR]>0.0){
             sd[k+nGR] = pow(qACCEPT, ((zn[k+nGR]/wn[k+nGR])-alpha_star));
	   }	 
           zn[k] = 0.0; 
           wn[k] = 0.0;
           zn[k+nGR] = 0.0; 
           wn[k+nGR] = 0.0;
         }
         rterm += dimG;
         cnt2+=nlGR[k]*dimG;
         cnt+=nlGR[k];
       }
     }

/***********************/
/* update alpha for PX */
/***********************/

     if(nalpha>0){

       for (i = 0; i < nalpha*nrowZ; i++){
         Xalpha->x[i] = 0.0;
       }

       cnt = ncolX;
       cnt2 = 0;
       for (k = 0; k < nG; k++){
         dimG = GRdim[k];
         for (i = 0; i < dimG; i++){
           if(PXtermsP[k]==1){   // parameter expanded 
             for (l = 0; l < nlGR[k]; l++){ 
               for (j = Worig->p[cnt+l]; j < Worig->p[cnt+l+1]; j++){ 
                 Xalpha->x[cnt2*nrowZ+Worig->i[j]] += Worig->x[j]*location->x[cnt+l];    
               }
             }
             cnt2 ++; 
           }
           cnt += nlGR[k];
         }
       }

       tXalpha = cs_transpose(Xalpha, true);

       XtmKRinv = cs_multiply(tXalpha, KRinv);   
	
       alphaM = cs_multiply(XtmKRinv, Xalpha);          

       alphaMME = cs_add(alphaM, Alphainv, 1.0, 1.0);                      

       alphaL = cs_chol(alphaMME, alphaS); 

       if(alphaL==NULL){
	 PutRNGstate();
         error("alpha equations singular: use a (stronger) prior for the alphas\n");
       }

       for(i=0; i<nalpha; i++){
          alphaastar->x[i] = rnorm(0.0,1.0);
       }               

       cs_ltsolve(AlphainvL->L, alphaastar->x);
       for(i=0; i<nalpha; i++){
          alphaastar->x[i]  += muAlpha->x[i];
       }      
                    
       cnt=0;
       for(k=nG; k<nGR; k++){
         dimG = GRdim[k];
         for(i=0; i<nlGR[k]; i++){
           for(j=0; j<dimG; j++){
             Grv[k]->x[j] = rnorm(0.0,1.0);
           }
           cs_ltsolve(GinvL[k]->L, Grv[k]->x);
           for(j=0; j<dimG; j++){
             alphazstar->x[j*nlGR[k]+i+cnt] = Grv[k]->x[j];     
           }
         }
         cnt += dimG*nlGR[k];  
       }  
  
// form W%*%theta for non expanded parameter

       for (i = 0; i < ny; i++){
           alphapred->x[i] = linky->x[i];
       }
       for (l = 0; l < ncolX; l++){
          for (j = W->p[l]; j < W->p[l+1]; j++){ 
             alphapred->x[W->i[j]] -= W->x[j]*location->x[l]; 
          }
       }

       cnt = ncolX;

       for (k = 0; k < nG; k++){
         dimG = GRdim[k];
         for (i = 0; i < dimG; i++){
           if(PXtermsP[k]==0){   // non-parameter expanded 
             for (l = 0; l < nlGR[k]; l++){ 
               for (j = W->p[cnt+l]; j < W->p[cnt+l+1]; j++){ 
                 alphapred->x[W->i[j]] -= W->x[j]*location->x[cnt+l];    
               }
             }
           }
           cnt += nlGR[k];
         }
       }

       cs_gaxpy(Xalpha, alphaastar->x, alphazstar->x);

       for (i = 0 ; i < ny ; i++){
         alphapred->x[i] -= alphazstar->x[i];       // form y - z*
       }

       alphaastar_tmp = cs_multiply(XtmKRinv, alphapred);                 // XtmKRinv(y - z*)

       for (i = 0 ; i < alphaastar_tmp->nzmax ; i++){
         alphalocation->x[alphaastar_tmp->i[i]] = alphaastar_tmp->x[i];
       }

       for (i = 0 ; i < nalpha; i++){
         alphalocation_tmp->x[i] = 0.0;
       }

       cs_ipvec (alphaS->pinv, alphalocation->x, alphalocation_tmp->x, alphaMME->n);	 
       cs_lsolve(alphaL->L,alphalocation_tmp->x);                                // x = L\x 
       cs_ltsolve (alphaL->L, alphalocation_tmp->x);		                 // x = L'\x 
       cs_pvec (alphaS->pinv, alphalocation_tmp->x, alphalocation->x, alphaMME->n);        // b = P'*x 

       for (i = 0 ; i < nalpha; i++){
         alphalocation->x[i] += alphaastar->x[i];
       }  

       cnt = ncolX;
       cnt2 = 0;
       for (k = 0; k < nG; k++){
         dimG = GRdim[k];
         for (i = 0; i < dimG; i++){
           if(PXtermsP[k]==1){    // parameter expanded 
             for (j = Worig->p[cnt]; j < Worig->p[cnt+nlGR[k]]; j++){ 
               W->x[j] = Worig->x[j]*alphalocation->x[cnt2];     
             }
             cnt2 ++; 
           }
           cnt += nlGR[k];
         }
       } 
       cs_spfree(Wt);
       Wt = cs_transpose(W, true);                               
     }

/***********************/
/*   store posterior   */
/***********************/

     if(itt%1000 == 0 && verboseP[0]){
         Rprintf("\n                      MCMC iteration = %i\n",itt);
       if(nMH>0){
         Rprintf("\n  Acceptance ratio for latent scores = %f\n", Eaccl/(nMH*1000.0));
         Eaccl = 0.0;
       }
       if(cp){
         for(i=0; i<nordinal; i++){
           Rprintf(" Acceptance ratio for cutpoint set %i = %f\n", i+1, accp[i]/1000.0);
           accp[i] = 0.0;
         }
       }
     }

     if(itt>=burnin && DICP[0]==1){
       mdbar *= (itt-burnin);
       mdbar += dbar;
       mdbar /= (itt-burnin+1.0);
       for(i=0; i< ny; i++){
	 mupred->x[i] *= (itt-burnin);
	 mupred->x[i] += pred->x[i];
	 mupred->x[i] /= (itt-burnin+1.0);
	 mulinky->x[i] *= (itt-burnin);
	 mulinky->x[i] += linky->x[i];
	 mulinky->x[i] /= (itt-burnin+1.0);
       }
       for(i=0; i<nGR; i++){    
	 dimG = GRdim[i];
	 for(j=0; j<(dimG*dimG); j++){
	   muG[i]->x[j] *= (itt-burnin);
     	   muG[i]->x[j] += G[i]->x[j];
	   muG[i]->x[j] /= (itt-burnin+1.0);
         }				 
       }
     }

     if((itt-burnin)%thin == 0 && itt>=burnin && itt!=nitt){
       if(DICP[0]==1){
         dbarP[post_cnt] = dbar;
       }
       if(cp){
         cnt =0;
         for(i=0; i<nordinal; i++){ 
           for(j=2; j<(ncutpointsP[i]-1); j++){ 
             CPP[cnt+post_cnt*ncutpoints_store] = oldcutpoints[cumsum_ncutpoints[i]+j];
             cnt++;
           }
         }
       }
       for (i = 0 ; i < ncolX ; i++){
         LocP[i+post_cnt*(ncolX+pr*ncolZ)] = location->x[i];
       }
       if(nL>0){
         for (i = 0 ; i < nL; i++){
           lambdaP[i+post_cnt*nL] = lambda[lambda_old]->x[i];        
         }
       }
       if(pr){
         if(nalpha>0){
           cnt=0;
           cnt2=0;
           for(k=0; k<nG; k++){
             dimG = GRdim[k];
             if(PXtermsP[k]==1){
               for(j=0; j<dimG; j++){
                 for(i=0; i<nlGR[k]; i++){              
                   LocP[ncolX+post_cnt*dimAS+cnt] = location->x[ncolX+cnt]*alphalocation->x[cnt2];
                   cnt++;
                 }
                 cnt2++;
               }               
             }else{
               for(j=0; j<dimG; j++){
                 for(i=0; i<nlGR[k]; i++){              
                   LocP[ncolX+post_cnt*dimAS+cnt] = location->x[ncolX+cnt];
                   cnt++;
                 }
               }
             }
           }
         }else{
           for (i = ncolX ; i < dimAS ; i++){
             LocP[i+post_cnt*dimAS] = location->x[i];
           }
         }
       }
       if(pl==1){
         for (i = 0 ; i < ny; i++){
           PLiabP[i+post_cnt*ny]  = linky->x[i];
         }
       }

       cnt=0;
       cnt2=0;
       for(i=0; i<nGR; i++){
         dimG = GRdim[i];
         if(PXtermsP[i]==0){
           for(j=0; j<(dimG*dimG); j++){
             VarP[cnt+post_cnt*tvc] = G[i]->x[j];
             cnt++;
           }
         }else{
           for(j=0; j<dimG; j++){              
             for(k=0; k<dimG; k++){
               VarP[cnt+post_cnt*tvc] = G[i]->x[j*dimG+k]*alphalocation->x[cnt2+j]*alphalocation->x[cnt2+k];
               cnt++;
             }
           }
           cnt2 += dimG;
         }
       }
       post_cnt++;
     }
   }

   if(DICP[0]==1){
     mdbar *= (itt-burnin+1.0);
     mdbar -= dbar;
     mdbar /= (itt-burnin);
     dbarP[post_cnt] = mdbar;
     dbarP[post_cnt+1] = dbar;
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
        cs_spfree(astar);
        cs_spfree(astar_tmp);
        cs_spfree(location);
        cs_spfree(location_tmp);
        cs_spfree(linky);
        cs_spfree(linki);
        cs_spfree(linki_tmp);
        cs_spfree(pred);
   	cs_spfree(mupred);
   	cs_spfree(mulinky);
        cs_spfree(predi);
        cs_spfree(dev);                              
        cs_spfree(pvB);
        cs_spfree(pmuB);  
        cs_spfree(Brv); 

        if(nL>0){                               
          cs_spfree(LambdaX);
          cs_spfree(pvL);
          cs_spfree(pmuL);
          cs_spfree(Lrv);
          cs_spfree(Lambda_tmp[lambda_old]);                                     
          cs_spfree(Lambda[lambda_old]);   
          cs_spfree(lambda[lambda_old]);   
          cs_spfree(lambdaI[0]);   
          cs_spfree(lambdaI[1]);   
          cs_spfree(I);
          cs_spfree(linky_orig);
          cs_spfree(ILY);
          cs_spfree(Y);
          cs_spfree(tY);
          cs_nfree(pvLL);
          cs_sfree(pvLS);
	  cs_sfree(LambdaS);
	  cs_nfree(LambdaLU);
          cs_sfree(tYKrinvYS);
          cs_nfree(tYKrinvYL);
        }
                                               
        if(Aexists){
          cs_spfree(A);               
          cs_spfree(bv);                                                              
          cs_spfree(bv_tmp);                                                            
        }


        cs_nfree(L);
        cs_sfree(S);
        cs_nfree(pvBL);
        cs_sfree(pvBS);

        for(i=0; i<nGR; i++){
	    cs_spfree(Ginv[i]);
	    cs_spfree(muG[i]);
	    cs_spfree(G[i]);
	    cs_spfree(Gtmp[i]);
 	    cs_spfree(Grv[i]);
	    cs_spfree(pG[i]);
	    cs_spfree(CM[i]);
	    cs_sfree(GinvS[i]);
	    cs_nfree(GinvL[i]);
            cs_spfree(KGinv[i]);
        }
        
        if(nalpha>0){
	  cs_spfree(Alphainv);
	  cs_spfree(muAlpha);
          cs_spfree(Xalpha);
          cs_spfree(tXalpha);
          cs_spfree(XtmKRinv);
          cs_spfree(alphaM);
          cs_spfree(alphaMME);
          cs_spfree(Worig);
          cs_sfree(alphaS);
          cs_nfree(alphaL);
          cs_sfree(AlphainvS);                  
          cs_nfree(AlphainvL);    
          cs_spfree(alphaastar);
          cs_spfree(alphazstar);
          cs_spfree(alphapred);
          cs_spfree(alphaastar_tmp);
          cs_spfree(alphalocation);
          cs_spfree(alphalocation_tmp);
        }

        for(i=nG; i<nGR; i++){
	    cs_spfree(propC[i]);
	    cs_spfree(propCinv[i]);
	    cs_sfree(propCinvS[i]);
	    cs_nfree(propCinvL[i]);
	    cs_spfree(propC[i+nGR]);
	    cs_spfree(propCinv[i+nGR]);
	    cs_sfree(propCinvS[i+nGR]);
	    cs_nfree(propCinvL[i+nGR]);
	    cs_spfree(muC[i]);
 	    cs_spfree(muC[i+nGR]);
        }
}
}

