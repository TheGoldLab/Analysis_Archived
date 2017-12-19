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

void HHMM_mixed_c(
        double *qexpr,
        double *Hexpr,
        double *l1pr,
        double *l2pr,
        double *Hpr,
        double *pHpr,
        double K,
        int N,
        int M,
        double *musgnpr
        )
        
{
    
    int n, m, m2;
    double q1,q2,H,p1,p2,l1,l2,Ev,*Qpr,*Ppr,pH,pH0,qexp,Hexp;
    
    Ppr = mxGetPr(mxCreateDoubleMatrix(M, 2, mxREAL));
    Qpr = mxGetPr(mxCreateDoubleMatrix(M, 2, mxREAL));
    
    for (m = 0; m < M; m++) {
        *(Ppr+m) = .5* *(pHpr+m);
        *(Ppr+m+M) = .5* *(pHpr+m);
    }
    
    for (n = 0; n < N; n++) {

        l1 = *(l1pr+n);
        l2 = *(l2pr+n);

        Ev = 0;
        
        for (m = 0; m < M; m++) {
            H = *(Hpr+m);

            *(Qpr+m)=l1*(*(Ppr+m)*(1-H)+*(Ppr+m+M)*H)*(1-K);
            *(Qpr+m+M)=l2*(*(Ppr+m)*H+*(Ppr+m+M)*(1-H))*(1-K);
            
            for (m2 = 0; m2 < M; m2++) {
                *(Qpr+m)+=l1*(*(Ppr+m2)*(1-H)+*(Ppr+m2+M)*H)*K* *(pHpr+m);
                *(Qpr+m+M)+=l2*(*(Ppr+m2)*H+*(Ppr+m2+M)*(1-H))*K* *(pHpr+m);
            }
            
            Ev += *(Qpr+m)+*(Qpr+m+M);

        }
        
        qexp = 0;
        for (m = 0; m < M; m++) {
            qexp+=*(Qpr+m)/Ev;
        }
        
        Hexp = 0;
        
        if (n==N | *(musgnpr+n+1)==0) {
            for (m = 0; m < M; m++) {
                *(Ppr+m) = *(Qpr+m)/Ev;
                *(Ppr+m+M) = *(Qpr+m+M)/Ev;
                Hexp += *(Hpr+m)* (*(Ppr+m)+*(Ppr+m+M));
            }
        } else if (*(musgnpr+n+1)==1) {
            Ev = 0;
            for (m = 0; m < M; m++) {
                *(Ppr+m+M)=0;
                Ev+=*(Qpr+m);
            }
        } else {
            Ev = 0;
            for (m = 0; m < M; m++) {
                *(Ppr+m)=0;
                Ev+=*(Qpr+m+M);
            }
        }
        
        if (*(musgnpr+n+1)==1) {
            for (m = 0; m < M; m++) {
                *(Ppr+m)=*(Qpr+m)/Ev;
                Hexp += *(Hpr+m)* *(Ppr+m);
            }
        } else if (*(musgnpr+n+1)==-1) {
            for (m = 0; m < M; m++) {
                *(Ppr+m+M)=*(Qpr+m+M)/Ev;
                Hexp += *(Hpr+m)* *(Ppr+m+M);
            }
        }
        
        *(Hexpr+n) = Hexp;
        *(qexpr+n) = qexp;

    }

}


void mexFunction( int nlhs, mxArray*plhs[],
        int nrhs, const mxArray*prhs[])
        
        /* VARIABLES ARE N, LPr, H, H0, alpha, J, H */
        
        
{
    double *qexpr, *Hexpr, *l1pr, *l2pr, *Hpr, *Evpr, K, *pHpr, *musgnpr;
    int N, M;
    
    l1pr = mxGetPr(prhs[0]);
    l2pr = mxGetPr(prhs[1]);
    Hpr = mxGetPr(prhs[2]);
    pHpr = mxGetPr(prhs[3]);
    K = mxGetScalar(prhs[4]);
    musgnpr = mxGetPr(prhs[5]);
    
    N = mxGetM(prhs[1]);
    M = mxGetM(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
    
    qexpr = mxGetPr(plhs[0]);
    Hexpr = mxGetPr(plhs[1]);
    
    /* Do the actual computations in a subroutine */
    HHMM_mixed_c(qexpr,Hexpr,l1pr,l2pr,Hpr,pHpr,K,N,M,musgnpr);
    
    
}


