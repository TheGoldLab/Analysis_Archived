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

void particle_filter_learnH(
        double K,
        double *H0pr,
        int H0len,
        double *xpr,
        double *musgnpr,
        double *q1samp_pr,
        int N,
        int partM,
        int simN,
        double F,
        double *etapr,
        double etalen,
        double *Hsamp_pr
        )
        
{
    
    int n, m, musgn, cp, musgnext, found, randm, simn, foundhigh, foundm;
    double LLR, eta, x, *partHpr, *partHprevpr, *cumwpr, *q1currpr, *q2currpr, H, q1samp, q2samp, *wpr, w, l1, l2, cumw, randval;
    
    partHpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    partHprevpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    q1currpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    q2currpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    wpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    cumwpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    
    for (simn = 0; simn < simN; simn++) {
        
        q1samp = .5;
        q2samp = .5;
        
        for (m=0; m<partM; m++) {
            randm = H0len * (double)rand() / (double)((unsigned)RAND_MAX + 1);
            *(partHpr+m) = *(H0pr+randm);
        }
        
        for (n = 0; n < N; n++) {
            
            x = *(xpr+n);
            randm = etalen * (double)rand() / (double)((unsigned)RAND_MAX + 1);
            eta = *(etapr+randm);
            x += eta;
            
            LLR = x/F;
            
            l1 = 1/(1+exp(-LLR));
            l2 = 1-l1;
            
            musgn = *(musgnpr+n);
            
            switch (musgn) {
                case 0:
                    for (m=0; m<partM; m++) {
                        H = *(partHpr+m);
                        *(q1currpr+m) = l1*((1-H) * q1samp + H * q2samp);
                        *(q2currpr+m) = l2*(H * q1samp + (1-H) * q2samp);
                    }
                    break;
                    
                case 1:
                    for (m=0; m<partM; m++) {
                        H = *(partHpr+m);
                        *(q1currpr+m) = l1*(1-H);
                        *(q2currpr+m) = l2*H;
                    }
                    break;
                    
                case -1:
                    for (m=0; m<partM; m++) {
                        H = *(partHpr+m);
                        *(q1currpr+m) = l1*H;
                        *(q2currpr+m) = l2*(1-H);
                    }
                    break;
                    
            }
            
            w = 0;
            
            for (m=0; m<partM; m++) {
                *(wpr+m) = *(q1currpr+m) + *(q2currpr+m);
                w+=*(wpr+m);
                *(cumwpr+m)=w;
            }
            
            
            found=0;
            randm = 0;
            randval = (double)rand() / (double)((unsigned)RAND_MAX + 1);
            
            while (found==0 & randm<partM) {
                found =  *(cumwpr+randm)/w>randval?1:0;
                randm++;
            }
            
            foundm = randm-1;

            q1samp = *(q1currpr+foundm) / *(wpr+foundm);
            q2samp = 1 - q1samp;
            
            *(q1samp_pr+simn+simN*n) = q1samp;
            *(Hsamp_pr+simn+simN*n) = *(partHpr+foundm);
            
            cp = -1;
            if (n<N) {
                musgnext = *(musgnpr+n+1);
                
                if (musgn!=0 & musgnext!=0 & musgn==musgnext) {
                    cp = 0;
                }
                else if (musgn!=0 & musgnext!=0 & musgn!=musgnext) {
                    cp = 1;
                }
                
            }
            
            /* determine new weights if there was feedback */
            
            switch (cp) {
                case 0:
                    w = 0;
                    for (m=0; m<partM; m++) {
                        H = *(partHpr+m);
                        *(wpr+m) = (1-H);
                        w += *(wpr+m);
                        *(cumwpr+m)=w;
                    }
                    break;
                    
                case 1:
                    w = 0;
                    for (m=0; m<partM; m++) {
                        H = *(partHpr+m);
                        *(wpr+m) = H;
                        w += *(wpr+m);
                        *(cumwpr+m)=w;
                    }
                    break;
                    
            }

            for (m=0; m<partM; m++) {
                *(cumwpr+m)=*(cumwpr+m)/w;
            }
            
            /* update particles */
            
            for (m=0; m<partM; m++) {
                
                randval = (double)rand() / (double)((unsigned)RAND_MAX + 1);
                
                if (randval<K){
                    randm = H0len * (double)rand() / (double)((unsigned)RAND_MAX + 1);
                    *(partHprevpr+m) = *(H0pr+randm);
                    
                } else {

                    
                    found=0;
                    randm = 0;
                    randval = (double)rand() / (double)((unsigned)RAND_MAX + 1);
                    
                    while (found==0 & randm<partM) {
                        found = *(cumwpr+randm)>randval?1:0;
                        randm++;
                    }
                    foundm = randm-1;
                    
                    *(partHprevpr+m) = *(partHpr+foundm);
                    
                }
                
            }
            
            for (m=0; m<partM; m++) {
                *(partHpr+m) = *(partHprevpr+m);
            }
            
        }
        
    }
    
}


void mexFunction( int nlhs, mxArray*plhs[],
        int nrhs, const mxArray*prhs[])
        
        /* VARIABLES ARE N, LPr, H, H0, alpha, J, H */
        
        
{
    double F, K, *H0pr, *xpr, *musgnpr, *q1samp_pr, partM, *etapr, etalen, *Hsamp_pr;
    int N, Hlen, simN;
    
    K = mxGetScalar(prhs[0]);
    H0pr = mxGetPr(prhs[1]);
    F = mxGetScalar(prhs[2]);
    xpr = mxGetPr(prhs[3]);
    musgnpr = mxGetPr(prhs[4]);
    partM = mxGetScalar(prhs[5]);
    simN = mxGetScalar(prhs[6]);
    etapr = mxGetPr(prhs[7]);
    
    N = mxGetM(prhs[3]);
    Hlen = mxGetM(prhs[1]);
    etalen = mxGetM(prhs[7]);
    
    plhs[0] = mxCreateDoubleMatrix(simN, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(simN, N, mxREAL);
    
    q1samp_pr = mxGetPr(plhs[0]);
    Hsamp_pr = mxGetPr(plhs[1]);
    
    /* Do the actual computations in a subroutine */
    particle_filter_learnH(K,H0pr,Hlen,xpr,musgnpr,q1samp_pr,N,partM,simN,F,etapr,etalen,Hsamp_pr);
    
    
}


