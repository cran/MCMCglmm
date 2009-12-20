#include "lpzF.h"

double lpzF(double x,  int lowertail)
{
double u;
    
                        if(x>(-7.5) && x<7.5){
	                   u= pnorm(x,  0.0, 1.0,lowertail, TRUE);
                         }else{                                                      // if l is extreme use Demidenko's Feller approximation
                           if(x<(-7.5)){
                             if(lowertail){
                               u = dnorm(-x, 0.0, 1.0, TRUE)-log(x);
                             }else{
                               u = dnorm(x, 0.0, 1.0, FALSE)/x;
                             }
                           }else{
                             if(lowertail){
                               u = -dnorm(x, 0.0, 1.0, TRUE)/x;
                             }else{
                               u = dnorm(x, 0.0, 1.0, FALSE)-log(x);
                             }
                           }
                         }   
                         return(u);

}
