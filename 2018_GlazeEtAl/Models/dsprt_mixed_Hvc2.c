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

void dsprt_mixed_Hvc(
        double *q1pr,
        double *q2pr,
        double *l1pr,
        double *l2pr,
        double *Hvcpr,
        int N,
        double *musgnpr
        )
        
{
    
    int n;
    double q1,q2,l1,l2,H,p1,p2,z;
    
    q1 = .5;
    q2 = .5;
    
    for (n = 0; n < N; n++) {
        
        H = *(Hvcpr+n);
        l1 = *(l1pr+n);
        l2 = *(l2pr+n);
        
        if (*(musgnpr+n)==0) {
            p1 = (1-H)*q1 + H*q2;
            p2 = (1-H)*q2 + H*q1;
        }
        else {
            p1 = *(musgnpr+n)==1?(1-H):H;
            p2 = 1-p1;
        }
        
        q1 = p1*l1;
        q2 = p2*l2;
        z = q1+q2;
        q1 = q1/z;
        q2 = q2/z;
       
        *(q1pr+n) = q1;
        *(q2pr+n) = q2;

    }
    
}


void mexFunction( int nlhs, mxArray*plhs[],
        int nrhs, const mxArray*prhs[])
        
        /* VARIABLES ARE N, LPr, H, H0, alpha, J, H */
        
        
{
    double *q1pr, *q2pr, *l1pr, *l2pr, *Hvcpr, *musgnpr;
    int N;
    
    l1pr = mxGetPr(prhs[0]);
    l2pr = mxGetPr(prhs[1]);
    Hvcpr = mxGetPr(prhs[2]);
    musgnpr = mxGetPr(prhs[3]);
    N = mxGetM(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
    
    q1pr = mxGetPr(plhs[0]);
    q2pr = mxGetPr(plhs[1]);
    
    /* Do the actual computations in a subroutine */
    dsprt_mixed_Hvc(q1pr,q2pr,l1pr,l2pr,Hvcpr,N,musgnpr);
    
    
}


