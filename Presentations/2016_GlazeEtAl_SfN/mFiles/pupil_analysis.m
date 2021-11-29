clear
clc
close all

% MAIN PUPIL ANALYSIS ROUTINE. 

% Set sequenceID to the stimulus sequence you want to process data for,
% and it will append processed data to pupildata.mat or
% pupildata_within.mat
%
% Note: this is poorly written exploratory code that should be cleaned up
% for future use if used.
%
%

load stimPilot_auditory_2AFC_071816.mat stimvc
stimvc2 = stimvc;

load stimPilot_auditory_2AFC_071116.mat stimvc
stimvc1 = stimvc;

load stimPilot_auditory_2AFC_080116.mat stimvc
stimvc3 = stimvc;


within = 0;
sequenceID = 3;


switch sequenceID
    
    case 2
        
        load stimPilot_auditory_2AFC_071816.mat z Hvc lvc stimvc
        
        titlelabs = {'H=0.3','H=0.99','H=0.3'};
        
        Hlen = [200 600 200];
        blockbegin = [100 300 850];
        blockend = [200 800 999];
        
    case 1
        
        
        titlelabs = {'H=0.01','H=0.3','H=0.01'};
        load stimPilot_auditory_2AFC_071116.mat z Hvc lvc stimvc
        
        Hlen = [400 200 400];
        blockbegin = [100 500 700];
        blockend = [400 600 999];
    case 3
        
        titlelabs = {'H=0.99','H=0.01','H=0.99','H=0.01'};
        load stimPilot_auditory_2AFC_080116.mat z Hvc lvc stimvc
        
        Hlen = [200 500 100 200];
        blockbegin = [100 600 750 850];
        blockend = [200 700 800 999];
        
        
        
end


dirnmarr{1} = ['.' filesep 'Pupil Data'];
Rbins = [1 2 3 6 inf];
%Rbins = [1 2 3 inf];
Rmin = 0;

stimspace = 1;

Rlabs = {'0','1','2-5','5+'};

subjindall = 0;
subjN = 0;

% LOOP THROUGH ALL DIRECTORIES AND FILES AND PRE-DETERMINE WHICH ARE VALID
isvalid = cell(numel(dirnmarr),1);
subjarr = cell(subjN,1);

for i = 1:numel(dirnmarr)
    dirnm = dirnmarr{i};
    
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

correctmat = zeros(subjN,1000);
correctarr = cell(numel(Hlen)-1,1);
Nmat = zeros(numel(Hlen),numel(Rbins)-1);

for Hind = 1:numel(correctarr)
    correctarr{Hind} = correctmat;
end

pupilarr = correctarr;

RTmed = zeros(subjN,1);
RThigh = RTmed;

pfollowedata = zeros(subjN,numel(Hlen));

correct_cp = zeros(subjN,1);
correct_ss = correct_cp;

Hmat = zeros(subjN,numel(Hlen));
psimat = Hmat;
lmat = Hmat;

totaldur = zeros(subjN,1);

N = numel(stimvc);

pupilmnmat = zeros(subjN,N);
pupilbasemat = pupilmnmat;
pupilchangemat = pupilmnmat;

dropoutvc = zeros(subjN,1);

pupilwin = 1.5;

minRTsubj = zeros(subjN,1);

RTall = [];
RTall2 = [];

% LOOP THROUGH ALL RELEVANT DIRECTORIES (ONLY ONE IN CURRENT VERSION)

for i = 1:numel(dirnmarr)
    
    dirnm = dirnmarr{i};
    
    fldir = dir(dirnm);
    flist = {fldir.name};
    
    % LOOP THROUGH ALL FILES
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
            
            if (within & numel(structarr)==2) | (~within & numel(structarr)>0)
                
                for structind = 1:numel(structarr)
                    
                    Data = structarr{structind};
                    
                    % IF THE STRUCTURE STIMULUS SEQUENCE MATCHES THE TARGET
                    % AND THERE IS A VALID PUPIL DATA FIELD
                    
                    if max(abs(stimvc-Data.StimList))==0 & isfield(Data,'EyeL')
                        
                        subjindall = subjindall + 1;
                        Choices = Data.Choices(1:end-1);
                        
                        subjarr{subjindall} = flnm;
                        
                        Stims = Data.StimList';
                        Stims = Stims(1:end);
                        
                        % stimspace is 0 by default. This determines
                        % whether to consider changepoints in stimulus
                        % space or state space.
                        
                        if stimspace
                            Targets = Stims;
                            CPs = Data.StimList(2:end)~=Data.StimList(1:end-1);
                        else
                            Targets = z;
                            CPs = z(2:end)~=z(1:end-1);
                        end
                        
                        %     CPs = z(2:end)~=z(1:end-1);
                        R = cp2run(CPs);
                        
                        %     Choices(2:end)=Stims(1:end-1);
                        
                        correct = Choices==Targets;
                        correctmat(subjindall,:) = correct;
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
                        
                        
                        % TAKE THE AVERAGE OF THE LEFT AND RIGHT EYE
                        pupil = .5*(Data.EyeL(:,3)+Data.EyeR(:,3));
                        
                        % CLEAN DATA AND CALCULATE % DROPOUT
                        [pupil,dropout] = pupilclean(pupil,2.2);
                        dropoutvc(subjindall) = dropout;
                        
                        pupilTimes = double(Data.EyeTime)/1000000;
                        dataOnset = Data.StimTimes(1);
                        
                        % FOR EACH TRIAL FIND PUPIL TRACE FOR PRE-TRIAL
                        % ONSET BASELINE AND WITHIN-TRIAL. 
                        for trialind = 1:N
                            stimOnset = Data.StimTimes(trialind);
                            pupilinds = find(pupilTimes>stimOnset & pupilTimes<(stimOnset+pupilwin));
                            
                            % PRE-TRIAL WINDOW IS 1 SECOND
                            baseinds = find(pupilTimes>stimOnset-1 & pupilTimes<stimOnset);
                            %                        pupilmat(subjindall,trialind) = max(pupil(pupilinds)) - pupil(pupilinds(1));
                            pupilmnmat(subjindall,trialind) = mean(pupil(pupilinds));
                            pupilbasemat(subjindall,trialind) = mean(pupil(baseinds));
                            
                            % DEFINING PUPIL CHANGE HERE AS AVERAGE
                            % DERIVATIVE. 
                            pupilchangemat(subjindall,trialind) = mean(diff(pupil(pupilinds)));
                        end
                        
                        trialpupil = pupilmnmat(subjindall,1:end-1);
                        trialpupil = trialpupil - mean(trialpupil);
          
                        RTs = diff(Data.ChoiceTimes);
                        RTall = [RTall;RTs];
                        
                        RTs = Data.ChoiceTimes - Data.StimTimes;
                        RTall2 = [RTall2;RTs];
                        
                        minRTsubj(subjindall) = min(RTs);
                        
                        RTmed(subjindall) = median(RTs);
                        RThigh(subjindall) = prctile(RTs,95);
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

minRTsubj = minRTsubj(1:subjindall);

payment = 8 + payment / 100;

pupilmnmat = pupilmnmat(1:subjindall,:);
pupilchangemat = pupilchangemat(1:subjindall,:);
pupilbasemat = pupilbasemat(1:subjindall,:);
subjarr = subjarr(1:subjindall);
dropoutvc = dropoutvc(1:subjindall);

% NORMALIZE MEASUREMENTS BY SUBTRACTING AVERAGES FROM FIRST 100 TRIALS
% ONWARD (OMITTING FIRST 1O0 TO ACCOUNT FOR TRANSIENT DIALATION DUE TO
% NOVELTY OF TASK ETC)
pupilmns = mean(pupilmnmat(:,100:end)')';
pupilmnnorm = pupilmnmat - repmat(pupilmns,1,N);

pupilmns = mean(pupilchangemat(:,100:end)')';
pupilchangenorm = pupilchangemat - repmat(pupilmns,1,N);

pupilmns = mean(pupilbasemat(:,100:end)')';
pupilbasenorm = pupilbasemat - repmat(pupilmns,1,N);

% AVERAGE ACROSS SUBJECTS
pupilmnmn = mean(pupilmnnorm);
pupilbasemn = mean(pupilbasenorm);
pupilchangemn = mean(pupilchangenorm);

correctmn = mean(correctmat);

% NOW APPEND DATA

switch sequenceID
    case 1
        pupilmnmn1 = mean(pupilmnnorm);
        pupilbasemn1 = mean(pupilbasenorm);
        pupilchangemn1 = mean(pupilchangenorm);
        pupilmnnorm1 = pupilmnnorm;
        pupilbasenorm1 = pupilbasenorm;
        pupilchangenorm1 = pupilchangenorm;
        dropoutvc1 = dropoutvc;
        
        correctmn1 = correctmn;
        correctmat1 = correctmat;
        
        if within
            save pupildata_within.mat -append correctmn1 correcmat1 pupilmnmn1 pupilchangemn1 pupilbasemn1 pupilmnnorm1 pupilchangenorm1 pupilbasenorm1
        else
            save pupildata.mat -append dropoutvc1 correctmn1 correctmat1 pupilmnmn1 pupilchangemn1 pupilbasemn1 pupilmnnorm1 pupilchangenorm1 pupilbasenorm1
        end
    case 2
        pupilmnmn2 = mean(pupilmnnorm);
        pupilbasemn2 = mean(pupilbasenorm);
        pupilchangemn2 = mean(pupilchangenorm);
        pupilmnnorm2 = pupilmnnorm;
        pupilbasenorm2 = pupilbasenorm;
        pupilchangenorm2 = pupilchangenorm;
        correctmn2 = correctmn;
        correctmat2 = correctmat;
        dropoutvc2 = dropoutvc;
        
        if within
            save pupildata_within.mat -append correctmn2 correcmat2 pupilmnmn2 pupilchangemn2 pupilbasemn2 pupilmnnorm2 pupilchangenorm2 pupilbasenorm2
        else
            save pupildata.mat -append dropoutvc2 correctmn2 correctmat2 pupilmnmn2 pupilchangemn2 pupilbasemn2 pupilmnnorm2 pupilchangenorm2 pupilbasenorm2
        end
        
    case 3
        pupilmnmn3 = mean(pupilmnnorm);
        pupilbasemn3 = mean(pupilbasenorm);
        pupilchangemn3 = mean(pupilchangenorm);
        pupilmnnorm3 = pupilmnnorm;
        pupilbasenorm3 = pupilbasenorm;
        pupilchangenorm3 = pupilchangenorm;
        correctmn3 = correctmn;
        correctmat3 = correctmat;
        dropoutvc3 = dropoutvc;
        
        if within
            save pupildata_within.mat -append correctmn3 correcmat3 pupilmnmn3 pupilchangemn3 pupilbasemn3 pupilmnnorm3 pupilchangenorm3 pupilbasenorm3
        else
            save pupildata.mat -append dropoutvc3 correctmn3 correctmat3 pupilmnmn3 pupilchangemn3 pupilbasemn3 pupilmnnorm3 pupilchangenorm3 pupilbasenorm3
        end
        
        
end
