function [B, totSig, R, pCha, sigmaU] = frugFun5_wRange(data, Hazard, noise, ...
    drift, likeWeight, trueRun, initGuess, inRun, dataRange)  

    %frugal bayesian belief updating MRN, UPENN, 2009
    %returns beliefs (B), uncertainty (totSig), run-length estimates (R), 
    %and change-point probability estimates (pCha) made by a reduced Bayesian model 
    %about the mean of the distribution that produced data (data).
 
    % notes for katie:  set trueRun = 0, noise = std 
    
    
    
    % the model expects the width of the generative distribution to be equal to noise,
    % , expects a (mean zero) drift with a rate = drift, and expects the
    % mean of the distribution to be repicked from a uniform distribution
    % with probability = Hazard.  The model can range from near-optimal
    % (likelihood weight = 1) to a fixed alpha delta rule (likelihood
    % weight = 0).  The model can track uncertainty by taking the mean of
    % the run length distribution (trueRun=1) or by matching the second
    % moment of the predictive distribution (trueRun=0).  The model also
    % takes an initial guess and updates beliefs from there.

    % MRN modified on 4/23/15 
    % Goals: 1) allow the model to take different ranges... ie 0-180
    % degrees.
    %        2) also, normalize the likelihood according to "out of range"
    %        part of CDF.
    
    
    
    
    
    sigmaE = noise;
    sigmaU=nan(1, length(data));
    totSig=nan(1, length(data));
    pCha=nan(1, length(data));
    B=nan(1, length(data));
    if nargin<7 || isempty(initGuess)
        if nargin>=9
            
            B(1)=nanmean(dataRange);
        else
            B(1)=150;
        end
    else
        B(1)=initGuess; 
    end
    
    
    R=nan(1, length(data));
    if exist('inRun', 'var')
        R(1)=inRun;
    else
        R(1)=1;
    end
    
    
    
    if nargin < 6 || isempty(trueRun)
        trueRun=0;
    end
    
    if nargin < 5 || isempty(likeWeight)
        likeWeight=1;
    end
    
    if nargin < 4 || isempty(drift)
        drift=0;
    end
    
    if Hazard >1
        Hazard=1;
    end
    
    if Hazard<0
        Hazard=0;
    end

%% loop through data making sequential predictions
        
        for i = 1:length(data)
            % part 1 get the expected distribution
            sigmaU(i)=sqrt((sigmaE./sqrt(R(i))).^2+drift.^2);  % changed things to squared and added in drift uncertainty
            R(i)=noise.^2 ./sigmaU(i).^2;                       % recompute R including drift uncertainty
            
            totSig(i)=sqrt(sigmaE.^2+sigmaU(i).^2);             % same deal
            % part 2 calculate probability of change
            pI=normpdf(data(i),B(i),totSig(i));
            
            
            
            % normalize to correct for probability outside of range.  fixed
            % an error in this code on 9-25-09...
            
            
            
            % compute total probability density within specified range:
            normalize=normcdf(dataRange(2), B(i), totSig(i)) - normcdf(dataRange(1), B(i), totSig(i)) ;
            % divide by total probability in range to deal with edge effects:         
            pI=pI./normalize;
            changLike=1./(dataRange(2)-dataRange(1));
            changeRatio=exp(likeWeight.*log(changLike./pI)+log(Hazard ./(1-Hazard)));
            
            if ~isfinite(changeRatio)
                pCha(i)=1;
                pNoCha=0;
            else
                pCha(i)=changeRatio./(changeRatio+1);
                pNoCha=1-pCha(i);
            end

            % part 3 update belief about mean
          
            yInt  = 1./(R(i)+1);    % now R can be really small if there is a big drift...
            slope = (1-yInt);
            Alph  = yInt+pCha(i).*slope;
            Delta = data(i)-B(i);
            B(i+1)=B(i)+Alph.*Delta;
            % part 4 update run length    
            if trueRun==1
            R(i+1)=(R(i)+1).*pNoCha + pCha(i);
            else 
            ss=pCha(i).*((sigmaE.^2)./1)+pNoCha.*((sigmaE.^2)./(R(i)+1))+pCha(i).*pNoCha.*((B(i)+yInt.*Delta)-data(i)).^2;
            R(i+1)=(sigmaE.^2)./ss;
            end
        
        end
end