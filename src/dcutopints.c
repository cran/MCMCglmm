#include "dcutpoints.h"

double dcutpoints(const cs *liab, double *yP, int *observed, int start,int finish, double *oldcutpoints, double *newcutpoints, int stcutpoints, int ncutpoints, double sdcp)
{
    int i,j,w;
    double llik = 0;
 
    for (j = 2 ; j < (ncutpoints-2); j++){
        llik += log(pnorm((oldcutpoints[stcutpoints+j+1]-oldcutpoints[j])/sdcp, 0.0, 1.0, TRUE,FALSE)-pnorm((newcutpoints[stcutpoints+j-1]-oldcutpoints[j])/sdcp, 0.0, 1.0, TRUE,FALSE));
        llik -= log(pnorm((newcutpoints[stcutpoints+j+1]-newcutpoints[j])/sdcp, 0.0, 1.0, TRUE,FALSE)-pnorm((oldcutpoints[stcutpoints+j-1]-newcutpoints[j])/sdcp, 0.0, 1.0, TRUE,FALSE));
    }

   llik += log(1.0-pnorm((newcutpoints[stcutpoints+ncutpoints-3]-oldcutpoints[stcutpoints+ncutpoints-2])/sdcp, 0.0, 1.0, TRUE,FALSE));
   llik -= log(1.0-pnorm((oldcutpoints[stcutpoints+ncutpoints-3]-newcutpoints[stcutpoints+ncutpoints-2])/sdcp, 0.0, 1.0, TRUE,FALSE));


    for (i = start ; i < finish; i++){
        w = yP[i];
        if(w>1 && observed[i]==1){
          if(w==(ncutpoints-1)){
            llik += log(1.0-pnorm((newcutpoints[stcutpoints+w-1]-liab->x[i]), 0.0, 1.0, TRUE,FALSE));
            llik -= log(1.0-pnorm((oldcutpoints[stcutpoints+w-1]-liab->x[i]), 0.0, 1.0, TRUE,FALSE));
          }else{
            llik += log(pnorm((newcutpoints[stcutpoints+w]-liab->x[i]), 0.0, 1.0, TRUE,FALSE)-pnorm((newcutpoints[stcutpoints+w-1]-liab->x[i]), 0.0, 1.0, TRUE,FALSE));
            llik -= log(pnorm((oldcutpoints[stcutpoints+w]-liab->x[i]), 0.0, 1.0, TRUE,FALSE)-pnorm((oldcutpoints[stcutpoints+w-1]-liab->x[i]), 0.0, 1.0, TRUE,FALSE));
          }
        }
    }
    return llik;
}
