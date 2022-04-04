clear
clc
close all

% INDIVIDUAL DIFFERENCES ANALYSIS. POORLY WRITTEN CODE THAT WAS
% EXPLORATORY, SHOULD BE CLEANED UP IF USED MORE.

% The SFN analysis did not rely on this except for some kind of compilation
% of the regcoeff variable. The following statement:
% [rmn,rchange,rbase] = pupil_anl_subj(subjind,1);
% measures trial-wise correlation between pupil average, baseline, etc. and
% model evidence. In theory we could see a subject-by-subject relationship
% between this correlation and learning.

% GENERATED BY behavioral_analysis.m
load behavioraldata.mat
rmat = zeros(size(Hmat1,1),2);
Hdev = zeros(size(Hmat1,1),1);

regcoeff = zeros(size(Hmat1,1)+size(Hmat2,1)+size(Hmat3,1),1);
psi = regcoeff;

regcoeff1 = zeros(size(Hmat1,1),1);
regcoeff2 = zeros(size(Hmat2,1),1);
regcoeff3 = zeros(size(Hmat3,1),1);

% ESTIMATE LEARNING AS REGRESSION OF FIT BLOCKWISE HAZARD RATE ON TASK

% DO IT FOR 1ST STIMULUS SEQUENCE
for subjind = 1:size(Hmat1,1)
    
    [rmn,rchange,rbase] = pupil_anl_subj(subjind,1);
    
    rmat(subjind,:) = rmn(1,:);
    Hdev(subjind) = median(abs(Hmat1(subjind,:)-Hobjmat1(subjind,:)));
    
    btmp = robustfit(Hobjmat1(subjind,:)',Hmat1(subjind,:)');
    regcoeff(subjind) = btmp(2);
    regcoeff1(subjind) = btmp(2);
%    psi(subjind) = log((1-pfollowedata1(subjind))/pfollowedata1(subjind));
    psi(subjind) = 1-pfollowedata1(subjind);
    
end

% DO IT FOR 2ND STIMULUS SEQUENCE
for subjind = 1:size(Hmat2,1)
    [rmn,rchange,rbase] = pupil_anl_subj(subjind,2);
    rmat(subjind+size(Hmat1,1),:) = rmn(1,:);
    Hdev(subjind+size(Hmat1,1)) = median(abs(Hmat2(subjind,:)-Hobjmat2(subjind,:)));
    
    btmp = robustfit(Hobjmat2(subjind,:)',Hmat2(subjind,:)');
    regcoeff(subjind+size(Hmat2,1)) = btmp(2);
    regcoeff2(subjind) = btmp(2);
 %   psi(subjind+size(Hmat1,1)) = log((1-pfollowedata2(subjind))/pfollowedata2(subjind));
    psi(subjind+size(Hmat1,1)) = 1-pfollowedata2(subjind);
    
end


% DO IT FOR 3RD STIMULUS SEQUENCE
for subjind = 1:size(Hmat3,1)
    [rmn,rchange,rbase] = pupil_anl_subj(subjind,3);
    rmat(subjind+size(Hmat1,1)+size(Hmat2,1),:) = rmn(1,:);
    Hdev(subjind+size(Hmat1,1)+size(Hmat2,1)) = median(abs(Hmat3(subjind,:)-Hobjmat3(subjind,:)));
    
    btmp = robustfit(Hobjmat3(subjind,:)',Hmat3(subjind,:)');
    regcoeff(subjind+size(Hmat1,1)+size(Hmat2,1)) = btmp(2);
    regcoeff3(subjind) = btmp(2);
   % psi(subjind+size(Hmat1,1)+size(Hmat2,1)) = log((1-pfollowedata3(subjind))/pfollowedata3(subjind));
    psi(subjind+size(Hmat1,1)+size(Hmat2,1)) = 1-pfollowedata3(subjind);
    
end


Hsubj = [Hmat1(:);Hmat2(:);Hmat3(:)];
Hobj = [Hobjmat1(:);Hobjmat2(:);Hobjmat3(:)];

figure
plot(Hobj,Hsubj,'ko')
xlim([0 1])
ylim([0 1])
hold on
plot([0 1],[0 1],'k--')

save indiffdata.mat
