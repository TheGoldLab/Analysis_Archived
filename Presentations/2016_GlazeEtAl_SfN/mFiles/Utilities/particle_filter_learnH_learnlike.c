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

/*particle_filter_learnH_learnlike(K,H0pr,Hlen,xpr,q1samp_pr,q2samp_pr,Hsamp_pr,likesamp_pr,N,partM,simN,like0pr,w_pr);*/

void particle_filter_learnH_learnlike(
        double K,
        double *H0pr,
        int H0len,
        double *xpr,
        double *q1samp_pr,
        double *q2samp_pr,
        double *Hsamp_pr,
        double *likesamp_pr,
        int N,
        double partM,
        int simN,
        double *like0pr,
        double *w_pr
        )
        
{
    
    int n, m, found, randm, simn;
    double eta, x, *partHpr, *partHprevpr, *partlikepr, *partlikeprevpr, *q1currpr, *q2currpr, *q1prevpr, *q2prevpr, H, q1samp, q2samp, *wpr, w, l1, l2, cumw, randval, like;
    
    for (simn = 0; simn < simN; simn++) {
        
        partHpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
        partlikepr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
        partHprevpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
        partlikeprevpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
        q1currpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
        q2currpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
        q1prevpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
        q2prevpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
        
        wpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
        
        q1samp = .5;
        q2samp = .5;
        
        for (m=0; m<partM; m++) {
            randm = H0len * (double)rand() / (double)((unsigned)RAND_MAX + 1);
            *(partHpr+m) = *(H0pr+randm);
            *(partlikepr+m) = *(like0pr+randm);
            *(q1prevpr+m) = .5;
            *(q2prevpr+m) = .5;
        }
        
        for (n = 0; n < N; n++) {
            
            x = *(xpr+n);
            
            for (m=0; m<partM; m++) {
                H = *(partHpr+m);
                like = *(partlikepr+m);
                l1 = x==1?like:1-like;
                l2 = 1-l1;
                
                *(q1currpr+m) = l1*((1-H) * *(q1prevpr+m) + H * *(q2prevpr+m));
                *(q2currpr+m) = l2*(H * *(q1prevpr+m) + (1-H) * *(q2prevpr+m));
            }
            
            w = 0;
            
            for (m=0; m<partM; m++) {
                *(wpr+m) = *(q1currpr+m) + *(q2currpr+m);
                w+=*(wpr+m);
            }
            
            cumw=0;
            found=0;
            randm = 0;
            randval = (double)rand() / (double)((unsigned)RAND_MAX + 1);
            
            while (found==0 & randm<partM) {
                cumw+=*(wpr+randm)/w;
                found = cumw>randval?1:0;
                randm++;
            }
            
            q1samp = *(q1currpr+randm-1) / *(wpr+randm-1);
            q2samp = *(q2currpr+randm-1) / *(wpr+randm-1);
            
            *(q1samp_pr+simn+simN*n) = q1samp;
            *(q2samp_pr+simn+simN*n) = q2samp;
            *(Hsamp_pr+simn+simN*n) = *(partHpr+randm-1);
            *(likesamp_pr+simn+simN*n) = *(partlikepr+randm-1);
            *(w_pr+simn+simN*n) = *(wpr+randm-1);
            
            /* update particles */
            
            
            for (m=0; m<partM; m++) {
                randval = (double)rand() / (double)((unsigned)RAND_MAX + 1);
                
                if (randval<K){
                    randm = H0len * (double)rand() / (double)((unsigned)RAND_MAX + 1);
                    *(partHprevpr+m) = *(H0pr+randm);
                    *(q1prevpr+m) = q1samp;
                    *(q2prevpr+m) = q2samp;
                    *(partlikeprevpr+m) = *(like0pr+randm);
                }
                
                else {
                    
                    cumw=0;
                    found=0;
                    randm = 0;
                    
                    while (found==0 & randm<partM) {
                        cumw+=*(wpr+randm)/w;
                        found = cumw>randval?1:0;
                        randm++;
                    }
                    
                    *(partHprevpr+m) = *(partHpr+randm-1);
                    *(q1prevpr+m) = *(q1currpr+randm-1)  / *(wpr+randm-1);
                    *(q2prevpr+m) = *(q2currpr+randm-1)  / *(wpr+randm-1);
                    *(partlikeprevpr+m) = *(partlikepr+randm-1);
                    
                }
            }
            
            for (m=0; m<partM; m++) {
                *(partHpr+m) = *(partHprevpr+m);
                *(partlikepr+m) = *(partlikeprevpr+m);
                
            }
            
        }
        
    }
    
}


void mexFunction( int nlhs, mxArray*plhs[],
        int nrhs, const mxArray*prhs[])
        
        /* VARIABLES ARE N, LPr, H, H0, alpha, J, H */
        
        
{
    double K, *H0pr, *xpr, *musgnpr, *q1samp_pr, *q2samp_pr, *Hsamp_pr, partM, *likesamp_pr, *like0pr, *w_pr;
    int N, Hlen, simN;
    
    K = mxGetScalar(prhs[0]);
    H0pr = mxGetPr(prhs[1]);
    like0pr = mxGetPr(prhs[2]);
    xpr = mxGetPr(prhs[3]);
    partM = mxGetScalar(prhs[4]);
    simN = mxGetScalar(prhs[5]);
    
    N = mxGetM(prhs[3]);
    Hlen = mxGetM(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(simN, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(simN, N, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(simN, N, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(simN, N, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(simN, N, mxREAL);
    
    q1samp_pr = mxGetPr(plhs[0]);
    q2samp_pr = mxGetPr(plhs[1]);
    Hsamp_pr = mxGetPr(plhs[2]);
    likesamp_pr = mxGetPr(plhs[3]);
    w_pr = mxGetPr(plhs[4]);
    
    /*
            double K,
        double *H0pr,
        int H0len,
        double *xpr,
        double *q1samp_pr,
        double *q2samp_pr,
        double *Hsamp_pr,
        double *likesamp_pr,
        int N,
        double partM,
        int simN,
        double *like0pr,
        double *w_pr
     **/
    
    /* Do the actual computations in a subroutine */
    particle_filter_learnH_learnlike(K,H0pr,Hlen,xpr,q1samp_pr,q2samp_pr,Hsamp_pr,likesamp_pr,N,partM,simN,like0pr,w_pr);
    
    
}


