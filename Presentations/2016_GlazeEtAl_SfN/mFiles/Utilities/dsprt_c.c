/*=================================================================
 *
 * DDDM ALGORITHM FOR ESTIMATION OF POSTERIORS, ETC IN HMM
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void dsprt_c(
        double *LPRpr,
        double *LLRpr,
        double H,
        int N
        )
        
{
    
    int n;
    double J,L;
    
    J = (1-H)/H;
    
    L = 0;
    
    for (n = 0; n < N; n++) {
        
        L += log(J+exp(-L)) - log(J+exp(L)) + *(LLRpr+n);
        
        *(LPRpr+n) = L;
        
    }
    
}


void mexFunction( int nlhs, mxArray*plhs[],
        int nrhs, const mxArray*prhs[])
        
        /* VARIABLES ARE N, LPr, H, H0, alpha, J, H */
        
        
{
    double *LPRpr, *LLRpr, H;
    int N;
    
    LLRpr = mxGetPr(prhs[0]);
    H = mxGetScalar(prhs[1]);
    N = mxGetM(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    
    LPRpr = mxGetPr(plhs[0]);
    
    /* Do the actual computations in a subroutine */
    dsprt_c(LPRpr,LLRpr,H,N);
    
    
}


