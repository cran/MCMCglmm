#include "rtnorm.h"

double rtnorm(double mu, double sd, double lower, double upper)
{
 
    double z,pz,u,slower,supper,tr,alpha;
    int sample=1;

    if(lower < -1e+16 || upper > 1e+16){

      if(upper > 1e+16){
        tr = (lower-mu)/sd;
      }else{
        tr = (mu-upper)/sd;
      }
      if(tr<0){                          // if sampling >0.5 of a normal density possibly quicker just to sample and reject
        while(sample==1){
          z = rnorm(0.0,1.0);
          if(z>tr){
            sample = 0;
          }
        }
      }else{
        alpha = (tr+sqrt((tr*tr)+4.0))/2.0;

        while(sample==1){
          z = rexp(alpha)+tr;
          pz = exp(-((z-alpha)*(z-alpha))/2.0);
          u = runif(0.0,1.0);
          if(u<pz){
            sample = 0;
          }
        }
      }
    }else{

      slower = (lower-mu)/sd;
      supper = (upper-mu)/sd;

      tr = pnorm(supper, 0.0,1.0, TRUE, FALSE)-pnorm(lower, 0.0,1.0, TRUE, FALSE);

      if(tr>0.5){                   // if sampling >0.5 of a normal density possibly quicker just to sample and reject
        while(sample==1){
          z = rnorm(0.0,1.0);
          if(z>slower && z<upper){
            sample = 0;
          }
        }
      }else{
        while(sample==1){
          z = runif(slower,supper);
          if(slower<=0.0 && 0.0<=supper){
            pz = exp(-z*z/2.0);
          }else{
            if(supper<0.0){
              pz = exp((supper*supper-z*z)/2.0);
            }else{
              pz = exp((slower*slower-z*z)/2.0);
            }
          }
          u = runif(0.0,1.0);
          if(u<pz){
            sample=0;
          }
        }
      }
    }
    if(lower < -1e+16){
      return(mu-z*sd);
    }else{
      return(z*sd+mu);
    }
}



