#include "MCMCglmm.h"

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

cs *cs_initialize(double *x, int *p, int *i, int n, int m, int nzmax){

        cs *X;
        int j;

        X = cs_spalloc(n, m, nzmax, 1, 0); 
                                                           
        for (j = 0 ; j < nzmax; j++){
          X->i[j] = i[j];
          X->x[j] = x[j];
        }
        for (j = 0 ; j <= m ; j++){
          X->p[j] = p[j];
        }

        return (cs_done (X, NULL, NULL, 1));	/* success; free workspace, return C */

}
