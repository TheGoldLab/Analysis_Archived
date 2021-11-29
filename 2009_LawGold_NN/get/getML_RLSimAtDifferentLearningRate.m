%% getML_RLSimAtDifferentLearningRate.m
%  Simulatie learning at different learning rates

%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [2.75]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:85];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA, W, MTind] = getML_RLSimPerceptualLearning11(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end

%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [91];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA, W, MTind] = getML_RLSimPerceptualLearning12(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end
%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81:85];
DFLAG     = 1;
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA, W, MTind] = getML_RLSimPerceptualLearningEasyDiff(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), DFLAG,RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end


%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-10 10];
SIMNUM    = [71:2];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA, W, MTind] = getML_RLSimPerceptualLearningFine(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end


%%

Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [0.2 0.6 0.8]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
nl        = [1.19]
RECOMPUTE = 1;



% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        for j = 1:length(nl)
            SIMNUM(kk)
            eln_(i)
            nl(j)
            time1            = clock;
            [FIRA WP WN]     = getML_RLSimPerceptualLearningNonLinPool(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), nl(j), RECOMPUTE);
            time2            = clock;
            time2-time1
            %clear FIRA W MTind
            pack
        end
    end
end





%%

Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [82];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA WP WN]           = getML_RLSimPerceptualLearningTwoLIPPools(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end





%%

Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [1 3 5 7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA WP WN]           = getML_RLSimPerceptualLearningNoNorm(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end

%%

Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [1.8]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA WP WN]           = getML_RLSimPerceptualLearningFitSig(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end

%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [0.4 0.6]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA WP WN]           = getML_RLSimPerceptualLearningGenRL01(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end


%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [0.55 0.6]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA WP WN]           = getML_RLSimPerceptualLearningGenRL00(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end



%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [2 4]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA WP WN]           = getML_RLSimPerceptualLearningAddNoise(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end


%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [0.2 0.7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA WP WN]           = getML_RLSimPerceptualLearningMulNoise(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end


%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [0.3 0.5 0.7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA WP WN]           = getML_RLSimPerceptualLearningOja(Monk, NSEN, NTUNE, eln_(i), DIRS, DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end


%%
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [20 40 60 80]; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
RECOMPUTE = 1;
% 25 7
for kk = 1:length(SIMNUM)
    for i = 1:length(eln_)
        SIMNUM(kk)
        eln_(i)
        time1            = clock;
        [FIRA WP WN]     = getML_RLSimPerceptualLearningTestBound(Monk, NSEN, NTUNE, eln_(i), DIRS, SIMNUM(kk), RECOMPUTE);
        time2            = clock;
        time2-time1
        %clear FIRA W MTind
        pack
    end
end

