clear
clc
close all

% COMPILES MODEL-FREE ANALYSIS AND FITS CHOICE DATA TO FIXED HAZARD RATE
% MODEL FOR EACH HAZARD-SPECIFIC
% BLOCK OF TRIALS FOR A GIVEN TASK SEQUENCE. A BUNCH OF OTHER JUNK IN HERE,
% THIS IS BADLY WRITTEN, EXPLORATORY CODE.

% Set sequenceID to the stimulus sequence you want to process data for,
% and it will append processed data to behavioraldata.maT

sequenceID = 3;

% KEEP THIS = 0
within = 0;

load stimPilot_auditory_2AFC_071116.mat stimvc
stimvc1 = stimvc;

load stimPilot_auditory_2AFC_071816.mat stimvc
stimvc2 = stimvc;

load stimPilot_auditory_2AFC_080116.mat stimvc
stimvc3 = stimvc;

switch sequenceID
    
    case 2
        load stimPilot_auditory_2AFC_071816.mat z Hvc lvc stimvc
        
    case 1
        load stimPilot_auditory_2AFC_071116.mat z Hvc lvc stimvc
        
    case 3
        load stimPilot_auditory_2AFC_080116.mat z Hvc lvc stimvc
        
end

Hlen = 200*ones(1,5);
blockbegin = [1   200   400   600   800];
blockend = [blockbegin(2:end)+1 999];

dirmnarr = cell(1,1);

dirmnarr{1} = ['Pupil Data'];

Rbins = [1 2 3 6 inf];
Rbins = [1 2 3 4 5];
Rmin = 0;

stimspace = 0;

Rlabs = {'0','1','2-5','5+'};

subjindall = 0;
subjN = 0;

isvalid = cell(numel(dirmnarr),1);

% LOOP THROUGH ALL DIRECTORIES AND FILES AND PRE-DETERMINE WHICH ARE VALID
for i = 1:numel(dirmnarr)
    dirnm = dirmnarr{i};
    
    fldir = dir(dirnm);
    flist = {fldir.name};
    
    isvalid{i} = zeros(numel(flist),1);
    
    for flind = 1:numel(flist)
        if length(flist{flind})>3 & strcmp(flist{flind}(end-3:end),'.mat')
            isvalid{i}(flind) = 1;
        end
    end
    
    subjN = subjN + sum(isvalid{i});
    
end
correctmn = zeros(subjN,1);
payment = correctmn;

correctmat = zeros(subjN,numel(Rbins)-1);
correctarr = cell(numel(Hlen)-1,1);
Nmat = zeros(numel(Hlen),numel(Rbins)-1);

for Hind = 1:numel(correctarr)
    correctarr{Hind} = correctmat;
end

RTmed = zeros(subjN,1);
RThigh = RTmed;

pfollowedata = zeros(subjN,numel(Hlen));

correct_cp = zeros(subjN,1);
correct_ss = correct_cp;

Hmat = zeros(subjN,numel(Hlen));
psimat = Hmat;
psistimvc = zeros(subjN,1);
lmat = Hmat;
Hobjmat = Hmat;
H0subj = zeros(subjN,1);
r0subj = H0subj;

totaldur = zeros(subjN,1);

% LOOP THROUGH ALL RELEVANT DIRECTORIES (ONLY ONE IN CURRENT VERSION)
for i = 1:numel(dirmnarr)
    
    dirnm = dirmnarr{i};
    
    fldir = dir(dirnm);
    flist = {fldir.name};
    
    % LOOP THROUGH ALL FILES)
    for flind = 1:numel(flist)
        flind
        
        if isvalid{i}(flind)
            fl = flist{flind};
            flnm = [dirnm '/' fl];
            vars = load(flnm);
            
            subjbegin = findstr(flnm,'Data_')+5;
            subjend = findstr(flnm,'.mat')-1;
            
            Datarrnm = fieldnames(vars);
            
            % GET ALL STRUCTURES FROM THE FILE.
            eval(['structarr=vars.' Datarrnm{1} ';'])
            
            if ~iscell(structarr)
                structarr2 = cell(1,1);
                structarr2{1} = structarr;
                structarr = structarr2;
                clear structarr2
            end
            
            clear Data
            
            typelist = zeros(numel(structarr),1);
            
            % FOR EACH STRUCTURE FIGURE OUT WHICH STIMULUS SEQUENCE WAS USED FOR
            
            for structind = 1:numel(structarr)
                
                Data = structarr{structind};
                
                if isfield(Data,'StimList')
                    
                    if max(abs(stimvc1-Data.StimList))==0
                        typelist(structind) = 1;
                    elseif max(abs(stimvc2-Data.StimList))==0
                        typelist(structind) = 2;
                    elseif max(abs(stimvc3-Data.StimList))==0
                        typelist(structind) = 3;
                    end
                end
            end
            
            if numel(structarr)>1
                if typelist(2)==typelist(1)
                    structarr = structarr(1);
                    typelist = typelist(1);
                end
            end
            
            emptyinds = typelist==0;
            structarr = structarr(~emptyinds);
            
            if (within & numel(structarr)>=2) | (~within & numel(structarr)>0)
                
                for structind = 1:numel(structarr)
                    
                    Data = structarr{structind};
                    
                    % IF THE STRUCTURE STIMULUS SEQUENCE MATCHES THE TARGET
                    % AND THERE IS A VALID PUPIL DATA FIELD
                    
                    if max(abs(stimvc-Data.StimList))==0 & isfield(Data,'EyeL')
                        
                        Data = structarr{structind};
                        
                        subjindall = subjindall + 1;
                        Choices = Data.Choices(2:end-1);
                        
                        Stims = Data.StimList';
                        Stims = Stims(2:end);
                        
                        % stimspace is 0 by default. This determines
                        % whether to consider changepoints in stimulus
                        % space or state space.
                        if stimspace
                            Targets = Stims;
                            CPs = Data.StimList(2:end)~=Data.StimList(1:end-1);
                        else
                            Targets = z(2:end);
                            CPs = z(2:end)~=z(1:end-1);
                        end
                        
                        
                        % DETERMINE RUN LENGTH, WHICH TRIAL HAD A
                        % CORRECT PREDICTION
                        R = cp2run(CPs);
                        
                        correct = Choices==Targets;
                        followedata = Choices(2:end)==Stims(1:end-1);
                        pfollowedata(subjindall) = mean(followedata);
                        
                        payment(subjindall) = 5*(sum(correct)-sum(~correct));
                        
                        cpanl = find(R(1:end-1)>=Rmin & R(2:end)==1)+1;
                        
                        Hbounds = [1 cumsum(Hlen)];
                        
                        correct_cp(flind) = mean(correct(R==2));
                        correct_ss(flind) = mean(correct(R>5));
                        
                        Stims2 = zeros(numel(Stims),1);
                        
                        Stims2(Stims==1) = 1;
                        Stims2(Stims==2) = -1;
                        Choices2 = double(Choices==1);
                        
                        Choices2 = Choices2(:);
                        Stims2 = Stims2(:);
                        
                        
                        
                        % BLOCK-WISE FITS TO SINGLE HAZARD RATE
                        % MODEL TO ESTIMATE LEARNING. EMPIRICALLY,
                        % THIS HAS BEEN GIVING WEIRD RESULTS IN
                        % WHICH SOME SUBJECTS ARE INTEGRATING
                        % (USING LOW H) EVEN DURING H=0.5 AND
                        % H>>0.5. CG THINKS THIS IS BECAUSE IN THIS
                        % TASK SOME ARE SLOWLY ADAPTING A SINGLE
                        % HAZARD RATE ESTIMATE BASED ON ENTIRE
                        % HISTORY OF SESSION DATA, SO FINAL
                        % DECISION POLICY REFLECTS AVERAGE HAZARD
                        % RATE OVER ENTIRE SESSION (OR EXPERIMENT).
                        
                        
                        
                        for blockind = 1:numel(Hlen)
                            bininds = blockbegin(blockind):blockend(blockind);
                            [correctmatmp,Nmatmp,correctmp] = ...
                                calculate_correct_mat(R(bininds),correct(bininds),Rbins,Rmin);
                            
                            correctmp(sum(Nmatmp)<2)=nan;
                            
                            correctarr{blockind}(subjindall,:) = correctmp;
                            
                            Nmat(blockind,:) = sum(Nmatmp);
                            pfollowedata(subjindall,blockind) = mean(followedata(bininds(1:end-1)));
                            
                            
                            LLR = Stims2(bininds).*log(lvc(bininds)./(1-lvc(bininds)));
                            
                            params = fit_dsprt_pred(LLR,Choices2(bininds));
                            
                            Hmat(subjindall,blockind) = params(1);
                            psimat(subjindall,blockind) = params(2);
                            Hobjmat(subjindall,blockind) = mean(CPs(bininds));
                            
                            
                        end
                        
                        % THIS IS THE MAIN MODEL-FREE DATA OF INTEREST:
                        % PROBABILITY CORRECT AS A FUNCTION OF RUN
                        % LENGTH AND HAZARD RATE BLOCK
                        [correctmatmp,Nmatmp,correctmp] = calculate_correct_mat(R,correct,Rbins,Rmin);
                        
                        correctmat(subjindall,:) = correctmp;
                        
                        RTs = diff(Data.ChoiceTimes);
                        
                        RTmed(subjindall) = median(RTs);
                        RThigh(subjindall) = prctile(RTs,95);
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

payment = 8 + payment / 100;

% NOW APPEND UPDATED BEHAVIORAL ANALYSIS

switch sequenceID
    
    case 1
        correctmat1 = correctmat(1:subjindall,:);
        Hmat1 = Hmat(1:subjindall,:);
        Hobjmat1 = Hobjmat(1:subjindall,:);
        psimat1 = psimat(1:subjindall,:);
        pfollowedata1 = pfollowedata(1:subjindall);
        Hnum1 = numel(Hlen);
        
    case 2
        correctmat2 = correctmat(1:subjindall,:);
        Hmat2 = Hmat(1:subjindall,:);
        Hobjmat2 = Hobjmat(1:subjindall,:);
        psimat2 = psimat(1:subjindall,:);
        pfollowedata2 = pfollowedata(1:subjindall);
        Hnum2 = numel(Hlen);
        
        
    case 3
        correctmat3 = correctmat(1:subjindall,:);
        Hmat3 = Hmat(1:subjindall,:);
        Hobjmat3 = Hobjmat(1:subjindall,:);
        psimat3 = psimat(1:subjindall,:);
        pfollowedata3 = pfollowedata(1:subjindall);
        Hnum3 = numel(Hlen);
        
end

save behavioraldata.mat -append


