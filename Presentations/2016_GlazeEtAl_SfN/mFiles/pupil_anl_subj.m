function [rmn,rchange,rbase] = pupil_anl_subj(subjind,stimind)
% THIS TAKES PUPIL DATA FROM A SINGLE SUBJECT
% AND ASKS WHETHER VARIOUS FEATURES (CHANGE FROM
% PRE-TRIAL BASELINE, PRE-TRIAL BASELINE OR SIMPLY MEAN) ARE CORRELATED
% WITH DIFFERENT LATENT VARIABLES FROM THE APPROXIMATED IDEAL OBSERVER MODEL
% IN modelevs.mat 

% See generate_mod_vars.m for description of model variables

load pupildata.mat

load modelevs.mat

switch stimind
    case 1
        CPPall = CPP1(100:end);
        logwall = logw1(100:end-1);
        logwall = abs(Hsamp1(101:end)-Hsamp1(100:end-1));
        pupilmnall = pupilmnnorm1(subjind,101:end)';
        pupilchangeall = pupilchangenorm1(subjind,101:end)';
        pupilbaseall = pupilbasenorm1(subjind,101:end)';
        Hsampall = Hsamp1(100:end);
        L = L1;
        LLR = LLR1;
        
    case 2
        CPPall = CPP2(100:end);
        logwall = logw2(100:end-1);
        logwall = abs(Hsamp2(101:end)-Hsamp2(100:end-1));
        pupilmnall = pupilmnnorm2(subjind,101:end)';
        pupilchangeall = pupilchangenorm2(subjind,101:end)';
        pupilbaseall = pupilbasenorm2(subjind,101:end)';
        Hsampall = Hsamp2(100:end);
        L = L2;
        LLR = LLR2;
        
    case 3
        CPPall = CPP3(100:end);
        logwall = logw3(100:end-1);
        logwall = abs(Hsamp3(101:end)-Hsamp3(100:end-1));
        pupilmnall = pupilmnnorm3(subjind,101:end)';
        pupilchangeall = pupilchangenorm3(subjind,101:end)';
        pupilbaseall = pupilbasenorm3(subjind,101:end)';
        Hsampall = Hsamp3(100:end);
        L = L3;
        LLR = LLR3;
end

Psi = L-LLR;
PriorCertAll = abs(Psi(100:end));
PostCertAll = abs(L(100:end));

Zmat = [CPPall,PriorCertAll];

[rmn,pmn] = corr([pupilmnall,logwall],'type','Spearman');
[rchange,pchange] = corr([pupilchangeall,logwall],'type','Spearman');
[rbase,pbase] = corr([pupilbaseall,logwall],'type','Spearman');
