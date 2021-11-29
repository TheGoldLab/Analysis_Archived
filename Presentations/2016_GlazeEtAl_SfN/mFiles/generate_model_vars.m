clear
clc
close all

% ROUTINE USED TO GENERATE MODEL EVIDENCE AND OTHER LATENT VARIABLES FOR
% PILOT ANALYSIS COMPARING WITH BEHAVIOR, PUPIL DIAMETER.

% NOTE: Throughout Chris Glaze has defined trial-by-trial model evidence as the
% normalization constant from the Baum-Welch algorithm, i.e. the sum of the
% unnormalized posteriors in state-space. The final value of this term at
% the end of an entire data sequence is what is effectively maximized in the
% Expectation-Maximization algorithm for learning the parameters of an HMM.
% Here we see trial-by-trial model evidence nicely fluctuate with changes
% in generative hazard rate, but CG never rigorously proved that this is
% in fact our most principled estimate of model evidence for online
% inference.

% THE THREE DIFFERENT STIMULUS SEQUENCES WE'VE COLLECTED DATA ON
stimfls = {'stimPilot_auditory_2AFC_071116.mat',...
    'stimPilot_auditory_2AFC_071816.mat',...
    'stimPilot_auditory_2AFC_080116.mat'};

for i = 1:3
    
    % stimvc contains binary vector indicating left/right, Hvc is
    % generative hazard rate and lvc is likelihood (Pr(stimulus|state))
    load(stimfls{i},'lvc','stimvc','Hvc');

    simN = 50;
    partM = 1000;
    K = .005;

    alpha = .5;
    
    priorM = 100000;

    likeweightvc = .8+.2*rand(priorM,1);
    
    H0vc = betarnd(alpha,alpha,priorM,1);
    
    x = sign(stimvc-mean(stimvc));

    % c-mex routine with sampling algorithm that learns and adapts to
    % changes in both hazard rate and the likelihood function. partM is
    % very high so we can get close to an ideal observer solution. Routine
    % generates same process simN times and we can just average over that.
    
    [q1,q2,Hsamp,likesamp,w] = particle_filter_learnH_learnlike(K,H0vc,likeweightvc,x,partM,simN);
    
    % To compare with task parameters you could plot something like
    
    %     figure;
    %     subplot(2,1,1)
    %     plot(Hvc)
    %     hold on
    %     plot(mean(Hsamp))
    %
    %     subplot(2,1,2)
    %     plot(lvc); hold on;
    %     plot(mean(likesamp))
    
    
    likesamp = mean(likesamp);
    Hsamp = mean(Hsamp);
    q1 = mean(q1);
    q2 = mean(q2);
    
    % Log posterior odds
    L = log(q1)-log(q2);
    
    % Log-likelihood ratios based on expected likelihoods
    LLR = x.*log(likesamp(:)./(1-likesamp(:)));
    
    % Model evidence
    logw = mean(log(w));
    Hsamp = Hsamp(:);
    L = L(:);
    logw = logw(:);
    
    % Changepoint probability based on expected hazard rate and 2nd moment
    % of posterior
    CPP =  1./(1+(1-Hsamp(1:end-1)).*(exp(LLR(2:end)+L(1:end-1))+1)./(Hsamp(1:end-1).*(exp(LLR(2:end))+exp(L(1:end-1)))));
    CPP = [1;CPP];
    
    eval(['Hsamp' num2str(i) '=Hsamp;'])
    eval(['LLR' num2str(i) '=LLR;'])
    eval(['L' num2str(i) '=L;'])
    eval(['logw' num2str(i) '=logw;'])
    eval(['CPP' num2str(i) '=CPP;'])
end

save modelevs.mat LLR1 LLR2 LLR3 L1 L2 L3 Hsamp1 Hsamp2 Hsamp3 logw1 logw2 logw3 CPP1 CPP2 CPP3
