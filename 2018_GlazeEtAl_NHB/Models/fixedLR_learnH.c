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

/*Bayesian_approx(K,pH0pr,Hpr,l1pr,l2pr,musgnpr,pHmatpr,q1matpr,q2matpr,N,M);*/

void fixedLR_learnH(
        double alpha,
        double H0,
        double *xpr,
        double *musgnpr,
        double *Lexp_pr,
        int N,
        int simN,
        double F,
        double *etapr,
        int etalen,
        double *Lvar_pr,
        double *Hexp_pr
        )
        
{
    
    int n, m, musgn, cp, musgnext, found, randm, simn, foundhigh, foundm, midpt;
    double LLR, eta, x, H, q1, q2, l1, l2, L, qmax, logqmax, neglogqmax, Omega, q1prev, q2prev;
    
    logqmax = 15;
    neglogqmax = -15;
    qmax = 1/(1+exp(neglogqmax));
    
    for (n = 0; n < N; n++) {
        *(Lexp_pr+n) = 0;
        *(Hexp_pr+n) = 0;
    }
    
    for (simn = 0; simn < simN; simn++) {
        
        q1prev = .5;
        q2prev = .5;
        H = H0;
        
        for (n = 0; n < N; n++) {
            
            x = *(xpr+n);
            randm = rand()%etalen;
            eta = *(etapr+randm);
            x += eta;
            
            LLR = x/F;
            
            l1 = 1/(1+exp(-LLR));
            l2 = 1-l1;

            q1 = l1*((1-H) * q1prev + H * q2prev);
            q2 = l2*(H * q1prev + (1-H) * q2prev);
  
            q1 = q1 / (q1 + q2);
            q2 = 1 - q1;
            
            L= q1>qmax? logqmax:(q2>qmax? neglogqmax: log(q1) - log(q2));
            
            *(Lexp_pr+n) += L;
            *(Lvar_pr+n) += pow(L,2);
            
            Omega = -1;
            if (n<N) {
                musgnext = *(musgnpr+n+1);
                
                if (musgn!=0 & musgnext!=0 & musgn==musgnext) {
                    Omega = 0;
                }
                else if (musgn!=0 & musgnext!=0 & musgn!=musgnext) {
                    Omega = 1;
                }  
                else {
                    Omega = 1/(1+((1-H)/H)*(l1*q1prev+l2*q2prev)/(l1*q2prev+l2*q1prev));
                }
            }
            
            H+=alpha*(Omega-H);
            *(Hexp_pr+n) += H;

            if (musgnext==1) {
                q1prev=1;
                q2prev=0;
            } else if (musgnext==-1) {
                q1prev=0;
                q2prev=1;
            } else {
                q1prev = q1;
                q2prev = q2;
            }
            
        }
        
    }
    
    for (n = 0; n < N; n++) {
        *(Lexp_pr+n) = *(Lexp_pr+n) / simN;
        *(Hexp_pr+n) = *(Hexp_pr+n) / simN;
        *(Lvar_pr+n) = *(Lvar_pr+n) / simN - pow(*(Lexp_pr+n),2);
    }
    
}


void mexFunction( int nlhs, mxArray*plhs[],
        int nrhs, const mxArray*prhs[])
        
        /* VARIABLES ARE N, LPr, H, H0, alpha, J, H */
        
        
{
    double F, H0, *xpr, *musgnpr, *Lexp_pr, *etapr, etalen, *Lvar_pr, alpha, *Hexp_pr;
    int N, simN;
    
    alpha = mxGetScalar(prhs[0]);
    H0 = mxGetScalar(prhs[1]);
    F = mxGetScalar(prhs[2]);
    xpr = mxGetPr(prhs[3]);
    musgnpr = mxGetPr(prhs[4]);
    simN = mxGetScalar(prhs[5]);
    etapr = mxGetPr(prhs[6]);
    
    N = mxGetM(prhs[3]);
    etalen = mxGetM(prhs[6]);
    
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(N, 1, mxREAL);
    
    Lexp_pr = mxGetPr(plhs[0]);
    Lvar_pr = mxGetPr(plhs[1]);
    Hexp_pr = mxGetPr(plhs[2]);
    
    /* Do the actual computations in a subroutine */
    fixedLR_learnH(alpha,H0,xpr,musgnpr,Lexp_pr,N,simN,F,etapr,etalen,Lvar_pr,Hexp_pr);
    
    
}


