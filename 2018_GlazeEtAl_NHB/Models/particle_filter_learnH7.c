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

void particle_filter_learnH6(
        double K,
        double *H0pr,
        int H0len,
        double *xpr,
        double *musgnpr,
        double *Lexp_pr,
        int N,
        int partM,
        int simN,
        double F,
        double *etapr,
        int etalen,
        double *Lvar_pr,
        double *Hexpr
        )
        
{
    
    int n, m, musgn, cp, musgnext, found, randm, simn, foundhigh, foundm, midpt;
    double randmax, cumpmid, LLR, eta, x, *partHpr, *partHprevpr, *cumwpr, *q1currpr, *q2currpr, H, q1samp, q2samp, *wpr, w, l1, l2, cumw, randval, randvalprev, Lsamp, qmax, logqmax, neglogqmax;
    
    logqmax = 30;
    neglogqmax = -30;
    qmax = 1/(1+exp(neglogqmax));
    
    partHpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    partHprevpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    q1currpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    q2currpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    wpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    cumwpr = mxGetPr(mxCreateDoubleMatrix(partM, 1, mxREAL));
    
    midpt = partM/2;
    randmax = (double)((unsigned)RAND_MAX + 1);
    
    for (n = 0; n < N; n++) {
        *(Lexp_pr+n) = 0;
        *(Lvar_pr+n) = 0;
    }
    
    for (simn = 0; simn < simN; simn++) {
        
        q1samp = .5;
        q2samp = .5;
        
        for (m=0; m<partM; m++) {
            randm = rand()%H0len;
            *(partHpr+m) = *(H0pr+randm);
        }
        
        for (n = 0; n < N; n++) {
            
            x = *(xpr+n);
            randm = rand()%etalen;
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
            randval = (double)rand() / randmax;
            
            while (found==0 & randm<partM) {
                found =  *(cumwpr+randm)/w>randval?1:0;
                randm++;
            }
            
            foundm = randm-1;

            q1samp = *(q1currpr+foundm) / *(wpr+foundm);
            q2samp = 1 - q1samp;
            *(Hexpr+n) += *(partHpr+foundm) / simN;
            
            Lsamp = q1samp>qmax? logqmax:(q2samp>qmax? neglogqmax: log(q1samp) - log(q2samp));
            
            *(Lexp_pr+n) += Lsamp / simN;
            *(Lvar_pr+n) += pow(Lsamp,2) / simN;

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
            
            cumpmid = *(cumwpr+midpt);
            
            /* update particles */
            
            for (m=0; m<partM; m++) {
                
                randval = (double)rand() / randmax;
                
                if (randval<K){
                    randm = rand()%H0len;
                    *(partHprevpr+m) = *(H0pr+randm);
                    
                } else {

                    
                    found=0;
                    randval = (double)rand() / randmax;
                    randm = randval>cumpmid?midpt:0;
                    
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
    
    for (n = 0; n < N; n++) {
        *(Lvar_pr+n) = *(Lvar_pr+n) - pow(*(Lexp_pr+n),2);
    }

}


void mexFunction( int nlhs, mxArray*plhs[],
        int nrhs, const mxArray*prhs[])
        
        /* VARIABLES ARE N, LPr, H, H0, alpha, J, H */
        
        
{
    double F, K, *H0pr, *xpr, *musgnpr, *Lexp_pr, partM, *etapr, etalen, *Lvar_pr, maxL, *Hexpr;
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
    
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(N, 1, mxREAL);
    
    Lexp_pr = mxGetPr(plhs[0]);
    Lvar_pr = mxGetPr(plhs[1]);
    Hexpr = mxGetPr(plhs[2]);
    
    /* Do the actual computations in a subroutine */
    particle_filter_learnH6(K,H0pr,Hlen,xpr,musgnpr,Lexp_pr,N,partM,simN,F,etapr,etalen,Lvar_pr,Hexpr);
    
    
}


