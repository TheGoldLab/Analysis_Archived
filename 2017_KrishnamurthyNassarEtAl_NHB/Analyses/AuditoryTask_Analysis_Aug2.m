%% MAIN ANALYSIS PIPELINE FOR:
%       Krishnamurthy, Nassar et al, 2017
%
% ----- INCLUDES THE 9 ADDITIONAL SUBJECTS (17th MAY 2016)
% The file includes analyses from both subject and model data.
% Trial analyses and model fitting are done elsewhere; let's maintain this
% file for generating the most current version of figures

%% FILE PATHS

clear
disp('Clearing Workspace!')
clc

%Base directory
base_dir = '~/Dropbox/auditoryTaskManuscript/forArchiving/';
% make sure these are added to your path
addpath(genpath(base_dir));

% input
pupilFileDir    = [base_dir 'Data/ProcessedPupil/'];  %directory with all the main task data files
trainingFileDir = [base_dir 'Data/trainingData/'];  %directory which has training task data file
modelData_dir = [base_dir 'Data/modelData/'];

traininFN = 'train_likeWidth(17-May-2016).mat'; % single auditory likelhihood?
%traininFN='train_likeWidth(20-Sep-2015).mat'; % in multiple bins??
trainingFile = [trainingFileDir traininFN];

modelData_file = [modelData_dir 'modelSimulation_wPredNoise(28-Jul-2016).mat'];


% output
% [We currently don't write outputs]


%%  ------------------- INITIALISATION CELL (RUN THIS BEFORE ANY ANALYSES)  ----------------------

% This cell loads all the subject data and populates the containers
% required for subsequent analysis. All common quantities are calculated in
% this cell.


%*********************** FLAG TO PROCESS DAY 1 DATA
include_day1_data = true;  % Do not change this. You can select/discard Day1 data later using a variable
%***********************

%*********************** FLAG TO RUN MODEL SIMULATIONS
forceRunModelSims=false;
%***********************

% load training (likelihood) data
load(trainingFile);

% This is the list of all the subjects  -- do not edit this list
% You can select a subset of subjects by means of includeID and excludeID
% (see below)
subjectList = {...
    'processedAbhinavD2S3fullTask(18-Jul-2014).mat',...
    'processedAbhinavD3S3fullTask(23-Jul-2014).mat',...
    'processedAbhinavD4S3fullTask(24-Jul-2014).mat',...
    'processedAbhinavD5S3fullTask(25-Jul-2014).mat',...
    'processedAbhinavD6S3fullTask(29-Jul-2014).mat',...
    'processedAmiraD4S3fullTask(09-Jul-2014).mat',...
    'processedAmiraD5S3fullTask(21-Jul-2014).mat',...
    'processedAmiraD6S3fullTask(22-Jul-2014).mat',...
    'processedAmiraD7S3fullTask(25-Jul-2014).mat',...
    'processedAmiraD8S3fullTask(28-Jul-2014).mat',...
    'processedArchanaD2S3fullTask(09-Jul-2014).mat',...
    'processedArchanaD3S3fullTask(10-Jul-2014).mat',...
    'processedArchanaD4S3fullTask(21-Jul-2014).mat',...
    'processedArchanaD5S3fullTask(22-Jul-2014).mat',...
    'processedArchanaD6S3fullTask(28-Jul-2014).mat',...
    'processedAyobamiD4S3fullTask(09-Jul-2014).mat',...
    'processedAyobamiD5S3fullTask(10-Jul-2014).mat',...
    'processedAyobamiD6S3fullTask(11-Jul-2014).mat',...
    'processedAyobamiD7S3fullTask(23-Jul-2014).mat',...
    'processedAyobamiD8S3fullTask(24-Jul-2014).mat',...
    'processedChristinaD2S3fullTask(30-Jul-2014).mat',...
    'processedChristinaD3S3fullTask(31-Jul-2014).mat',...
    'processedDawnD2S3fullTask(04-Aug-2014).mat',...
    'processedDawnD3S3fullTask(05-Aug-2014).mat',...
    'processedDawnD4S3fullTask(06-Aug-2014).mat',...
    'processedDawnD5S3fullTask(08-Aug-2014).mat',...
    'processedDhruvaD2S2fullTask(17-Jul-2014).mat',...
    'processedDhruvaD3S3fullTask(21-Jul-2014).mat',...
    'processedHaftomD2S3fullTask(16-Jul-2014).mat',...
    'processedHaftomD3S3fullTask(17-Jul-2014).mat',...
    'processedHaftomD4S3fullTask(23-Jul-2014).mat',...
    'processedHaftomD5S3fullTask(25-Jul-2014).mat',...
    'processedHaftomD6S3fullTask(28-Jul-2014).mat',...
    'processedJenniferD2S3fullTask(01-Aug-2014).mat',...
    'processedJenniferD3S3fullTask(02-Aug-2014).mat',...
    'processedJenniferD4S3fullTask(05-Aug-2014).mat',...
    'processedJessicaD2S3fullTask(29-Jul-2014).mat',...
    'processedJessicaD3S3fullTask(03-Aug-2014).mat',...
    'processedJessicaD4S3fullTask(04-Aug-2014).mat',...
    'processedJulieD2S3fullTask(17-Jul-2014).mat',...
    'processedJulieD3S3fullTask(18-Jul-2014).mat',...
    'processedJustynD2S3fullTask(25-Jul-2014).mat',...
    'processedJustynD3S3fullTask(28-Jul-2014).mat',...
    'processedJustynD4S3fullTask(29-Jul-2014).mat',...
    'processedKatelynD2S3fullTask(10-Jul-2014).mat',...
    'processedKatelynD3S3fullTask(11-Jul-2014).mat',...
    'processedKatelynD4S3fullTask(16-Jul-2014).mat',...
    'processedKatelynD5S3fullTask(17-Jul-2014).mat',...
    'processedKatelynD6S3fullTask(19-Jul-2014).mat',...
    'processedKelseyD6S3fullTask(09-Jul-2014).mat',...
    'processedKelseyD7S3fullTask(10-Jul-2014).mat',...
    'processedKelseyD8S3fullTask(17-Jul-2014).mat',...
    'processedKelseyD9S3fullTask(04-Aug-2014).mat',...
    'processedKenyaD2S3fullTask(01-Aug-2014).mat',...
    'processedKenyaD3S3fullTask(03-Aug-2014).mat',...
    'processedKenyaD4S3fullTask(05-Aug-2014).mat',...
    'processedKenyaD5S3fullTask(07-Aug-2014).mat',...
    'processedKenyaD6S3fullTask(08-Aug-2014).mat',...
    'processedMarinaD2S3fullTask(01-Aug-2014).mat',...
    'processedMarinaD3S3fullTask(02-Aug-2014).mat',...
    'processedMarinaD4S3fullTask(03-Aug-2014).mat',...
    'processedMarinaD5S3fullTask(05-Aug-2014).mat',...
    'processedMarinaD6S3fullTask(08-Aug-2014).mat',...
    'processedRafaelD2S3fullTask(31-Jul-2014).mat',...
    'processedRafaelD3S3fullTask(01-Aug-2014).mat',...
    'processedRafaelD4S3fullTask(04-Aug-2014).mat',...
    'processedRafaelD5S3fullTask(05-Aug-2014).mat',...
    'processedRobD2S3fullTask(23-Jul-2014).mat',...
    'processedRobD3S3fullTask(24-Jul-2014).mat',...
    'processedSharonD2S3fullTask(31-Jul-2014).mat',...
    'processedSharonD3S3fullTask(01-Aug-2014).mat',...
    'processedSharonD4S3fullTask(02-Aug-2014).mat',...
    'processedSharonD5S3fullTask(03-Aug-2014).mat',...
    'processedSharonD6S3fullTask(04-Aug-2014).mat',...
    'processedValerieD2S3fullTask(09-Jul-2014).mat',...
    'processedValerieD3S3fullTask(10-Jul-2014).mat',...
    'processedValerieD4S3fullTask(15-Jul-2014).mat',...
    'processedValerieD5S3fullTask(16-Jul-2014).mat',...
    'processedValerieD6S3fullTask(18-Jul-2014).mat',...
    };

%ID = 1:length(subjectList);
ID = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 ...
    6 6 6 6 7 7 8 8 8 8 8 9 9 9 10 10 10 11 11 12 12 12,...
    13 13 13 13 13 14 14 14 14 15 15 15 15 15 16 16 16 16 16,...
    17 17 17 17 18 18 19 19 19 19 19 20 20 20 20 20,...
    ];

% Adding the subjects from the new data collection (May 2016)
% Keeping this separate for easier bookeeping
subjectList = {subjectList{:},...
    'processedCarolineD2S3fullTask(03-May-2016).mat',...
    'processedCarolineD3S3fullTask(04-May-2016).mat',...
    'processedCarolineD4S3fullTask(06-May-2016).mat',...
    'processedDiarraD2S3fullTask(05-May-2016).mat',...
    'processedDiarraD3S3fullTask(06-May-2016).mat',...
    'processedDiarraD4S3fullTask(09-May-2016).mat',...
    'processedHarrisonD2S3fullTask(21-Apr-2016).mat',...
    'processedHarrisonD3S3fullTask(26-Apr-2016).mat',...
    'processedHarrisonD4S3fullTask(03-May-2016).mat',...
    'processedMACD2S3fullTask(12-May-2016).mat',...
    'processedMACD3S3fullTask(13-May-2016).mat',...
    'processedMACD4S3fullTask(14-May-2016).mat',...
    'processedMauriceD2S3fullTask(19-Apr-2016).mat',...
    'processedMauriceD3S3fullTask(20-Apr-2016).mat',...
    'processedMauriceD4S3fullTask(22-Apr-2016).mat',...
    'processedMBD2S3fullTask(12-May-2016).mat',...
    'processedMBD3S3fullTask(13-May-2016).mat',...
    'processedMBD4S3fullTask(14-May-2016).mat',...
    'processedMKD2S3fullTask(11-May-2016).mat',...
    'processedMKD3S3fullTask(16-May-2016).mat',...
    'processedMKD4S3fullTask(17-May-2016).mat',...
    'processedSylviaD2S3fullTask(03-May-2016).mat',...
    'processedSylviaD3S3fullTask(05-May-2016).mat',...
    'processedSylviaD4S3fullTask(06-May-2016).mat',...
    'processedZBD2S3fullTask(13-May-2016).mat',...
    'processedZBD3S3fullTask(16-May-2016).mat',...
    'processedZBD4S3fullTask(17-May-2016).mat',...
    };
ID = [ID 21 21 21 22 22 22 23 23 23 24 24 24 25 25 25 ...
    26 26 26 27 27 27 28 28 28 29 29 29];


% ----------------- Do you want to include DAY1 data?
if(include_day1_data)
    
    subjectList = {subjectList{:},...
        'processedAbhinavD1S3fullTask(15-Jul-2014).mat',...
        'processedChristinaD1S3fullTask(29-Jul-2014).mat',...
        'processedDawnD1S3fullTask(31-Jul-2014).mat',...
        'processedDhruvaD1S3fullTask(16-Jul-2014).mat',...
        'processedJenniferD1S3fullTask(30-Jul-2014).mat',...
        'processedJessicaD1S3fullTask(28-Jul-2014).mat',...
        'processedJulieD1S3fullTask(15-Jul-2014).mat',...
        'processedJustynD1S3fullTask(24-Jul-2014).mat',...
        'processedKatelynD1S3fullTask(09-Jul-2014).mat',...
        'processedKenyaD1S3fullTask(31-Jul-2014).mat',...
        'processedMarinaD1S3fullTask(24-Jul-2014).mat',...
        'processedRafaelD1S3fullTask(30-Jul-2014).mat',...
        'processedRobD1S3fullTask(21-Jul-2014).mat',...
        'processedSharonD1S3fullTask(30-Jul-2014).mat',...
        'processedCarolineD1S3fullTask(28-Apr-2016).mat',...
        'processedDiarraD1S3fullTask(04-May-2016).mat',...
        'processedHarrisonD1S3fullTask(20-Apr-2016).mat',...
        'processedMACD1S3fullTask(10-May-2016).mat',...
        'processedMauriceD1S3fullTask(18-Apr-2016).mat',...
        'processedMBD1S3fullTask(11-May-2016).mat',...
        'processedMKD1S3fullTask(10-May-2016).mat',...
        'processedSylviaD1S3fullTask(28-Apr-2016).mat',...
        'processedZBD1S3fullTask(12-May-2016).mat',...
        };
    
    ID = [ID, 1 5 6 7 9 10 11 12 13 15 16 17 18 19 21 22 23 24 25 26 27 28 29];
end

% ----- Include 2 subjects with bad pupil for purely behavioural analysis?
uniqueID = unique(ID);
uniqueStds = [10 20];

% OK, so here's how to do the analysis on a subset of subjects:
% If you want to do the analsysi on, say, subjects  1, 4  and 5, then set ifReject to
% false and set includeID = [ 1 4 5].
% If you want to do the analysis on all subjects other than 1 4 & 5, the
% set ifReject to true, and set excludeID = [1 4 5]

ifReject = false;
includeID = uniqueID;
%includeID =  [10:20];%[14 15 18 19 20]; %10 11 12 13
%includeID = [1];

%excludeID = [17 16 14 3 7];

if(ifReject)
    for i1 = 1:length(excludeID)
        subjectList(ID == excludeID(i1)) = [];
        ID(ID == excludeID(i1)) = [];
    end
    
else
    for i1 = 1:length(uniqueID)
        if(~any(includeID == uniqueID(i1)))
            subjectList(ID == uniqueID(i1)) = [];
            ID(ID == uniqueID(i1)) = [];
        end
    end
    
end

if(length(ID) ~= length(subjectList))
    error('Enter the correct ID list');
end
subjectCount = zeros(length(uniqueID),1);


% Number of trials-after-CP
numTAC = 6;


% OK, now we initialise containers for the subject data.
% Here's the naming convention we use:
% variables with the prefix 'all' are "sound-level" containers
% variables with the prefix 'fix' are containers for data from properly
% fixated probe trials
allOutcome = [];
allWagerOutcomes = [];
allEstimate = [];
allMean = [];
allStd = [];
allCPtrialNo = [];
allIsCPtrials = [];
allWagerTrialIdx = [];
allPredictions = [];
allUpdate = [];
allTAC = [];
allWagerTAC = [];
allWager = [];
allWagerMeans = [];
allTA_bCP=[];        %trials after "big CP"
allbCP=[];           % logical for "big CP"

%adding these to merge model fitting code (Kamesh)
allSubjNum=[];
allWager=[];
allSesNum=[];
allTrialNum = [];
allIsProbe = [];
allSoundBlkNum = [];
allSessID = [];
allDay1Sess = []; %is this a Day1 session

%containers for probe level info, but includes *all* probes, not just fixated ones!!
probeBlkNum = [];
probeSessionID = [];
probeOutcomes = [];
probeEstimate = [];
probePrediction = [];
probeMean = [];
probeStd = [];
probeTAC = [];    % TAC of the probe trial
probeWager = [];
probeSubjID = [];
probeIsFixed = [];
probebCPTAC = [];  % big change-points
probeDay1Sess = [];
probeCentralOutcomes = [];


%subjID = [];

%Pupil related variables
allZSFilteredPupil = [];  %z-scored pupil
allFilteredPupil  = [];   %filtered pupil
allIsFixed = [];      % was fixation broken in the trial(Boolean)
soundOffset = 0.60; %-- time delay(in s) between calling play() and hearding the sound

spatialNormPercpErr = [];
spatialNormPredErr = [];
allRawPredErr = [];
allRawPercpErr = [];

fixDay1Sess = [];
fixSessionID = [];  % what's the serial number of the session (each subject did between 3-5 sessions)
fixOutcomes = [];
fixEstimate = [];
fixPrediction = [];
fixMean = [];
fixStd = [];
fixTAC = [];    % TAC of the probe trial
fixbCPTAC = [];  %TAC of probe trial considering only "big" change-points
fixIsFixed = []; % logical array indicating which trials are properly fixated trials
fixRunLengthSTD = [];  %STD for the current runlength
fixCentralOutcomes = [];

% Various pupil containers  (don't get rid of these w/o checking)
fixRawPupil = [];   %pupil before any further processing
fixRawTrendSubPupil = []; %pupil with linear baseline trend removed in a day
fixMeanSubtrPupil = [];  %mean subtracted per block
fixFilteredPupil  = [];  %mean subtracted per block w/ baseline trend removed
fixPrevTrialRawPupil = [];  %traces from trial before and after
fixNextTrialRawPupil = [];
fixPrevTrialRawTrendSubPupil = [];  %traces from trial before and after
fixNextTrialRawTrendSubPupil = [];
fixPrevTrialMeanSubPupil = [];  %traces from trial before and after
fixNextTrialMeanSubPupil = [];
fixPrevTrialMSZSPupil = [];
fixNextTrialMSZSPupil = [];


fixRunLengths = [];      % number of sounds before a probe trial
fixLastOutcomeTAC = [];  %TAC of the outcome just before the current one
fixSubjID = [];
fixNormPercpErr = [];    %perceptual error after centering and normalisation
fixNormPredErr = [];     %prediction error after centering and normalisation
fixWager = [];
fixExtremeWager = [];    %logical indicating sessions where wager was extreme
fixTrialNo = [];
fixTsinceLastProbe = [];
fixWagerOutcomeTimes = [];  %time of probe trial since the block beginning
fixTimeToPredict = [];  %time spent in the prediction stage
fixTimeToEstimate = []; %time spent in estimation stage
fixTimeToBet = []; %time to bet
fixInterProbeInterval = [];
fixBaseline = [];
fixRawBaseline = [];



numPupSamples = 150;  %How many pupil samples do you want?

d=fdesign.lowpass('Fp,Fst,Ap,Ast',2,4,1,50,60);
myFilt1 = design(d,'butter','MatchExactly','passband');

d2=fdesign.lowpass('Fp,Fst,Ap,Ast',1,3,1,50,60);
myFilt2 = design(d2,'butter','MatchExactly','passband');
%***** Correcting perceptual and prediction error for spatial effects***********
% Method: divide the arc into bins of 30 degrees and calculate the bias
% and variance of the perceptual and prediction errors in these bins.
% Normalized the prediction and perceptual errors in all bins to have
% zero bias (what about the variance?)
angleBins = 0:30:180;
spatialPredErrBias = zeros(length(subjectList),length(angleBins)-1);
spatialPredErrStd = zeros(length(subjectList),length(angleBins)-1);
spatialPercpErrBias = zeros(length(subjectList),length(angleBins)-1);
spatialPercpErrStd = zeros(length(subjectList),length(angleBins)-1);


for i1 = 1:length(subjectList)
    
    subjectCount(uniqueID == ID(i1)) = subjectCount(uniqueID == ID(i1)) + 1;
    
    
    aux = load([pupilFileDir subjectList{i1}]);
    outcome = aux.outcome;
    allOutcome = [allOutcome outcome];
    
    if(i1 > 104)
        allDay1Sess = [allDay1Sess ; 1*ones(length(outcome),1)];
    else
        allDay1Sess = [allDay1Sess ; 0*ones(length(outcome),1)];
    end
    
    
    %get the indicies of wager outcomes (which of the indices in
    %allOutcomes are wager trials?)
    wagerTrialIdx = aux.wagerTrialIdx;
    allWagerTrialIdx = [allWagerTrialIdx, ...
        wagerTrialIdx + length(allStd)];        %THIS IS WRONG. DON'T USE!!
    
    %now record the outcomes corresponding to wager trials
    wagerOutcomes = outcome(wagerTrialIdx);
    allWagerOutcomes = [allWagerOutcomes wagerOutcomes];
    
    estimate = aux.estimate;
    allEstimate = [allEstimate estimate];
    
    prediction = aux.prediction;
    allPredictions = [allPredictions prediction];
    
    % THIS DOESN'T MAKE SENSE -- don't use!
    %update = [diff(prediction) 0];
    %allUpdate = [allUpdate update];
    
    Mean = aux.Mean;
    allMean = [allMean Mean];
    
    wagerMeans = aux.wagerMeans;
    allWagerMeans = [allWagerMeans wagerMeans];
    
    Std = aux.Std;
    allStd = [allStd Std];
    
    
    isCpTrial = aux.isCpTrial;
    allIsCPtrials = [allIsCPtrials isCpTrial];
    
    
    cpTrialNo = find(isCpTrial);
    allCPtrialNo = [allCPtrialNo cpTrialNo];
    
    wager = aux.wager;
    %DEBUG -- this is strange but for some reason Kenya DAY 1 has an extra
    %wager entry. No idea how it got there, but should be harmless
    
    if(strcmp(subjectList{i1},'processedKenyaD1S3fullTask(31-Jul-2014).mat'))
        wager = wager(1:end-1);
    end
    allWager = [allWager wager];
    
    
    wagerTAC = aux.wagerTAC;
    wagerTrialIdx = aux.wagerTrialIdx;
    
    allWagerTAC = [allWagerTAC wagerTAC];
    
    
    runLengthStd = nan(length(wagerTAC),1);
    for i00 = 1:length(wagerTAC)-1
        runLengthStd(i00+1) = std(outcome(wagerTrialIdx(i00):wagerTrialIdx(i00+1)));
    end
    
    
    % ADDITIONAL VARIABLES ADDED BY MRN 10-12-14 **********************
    % get the big change points -- where the mean changes by more than
    % 2*std
    
    bThresh= 2;
    bCP= isCpTrial&[0 abs(diff(Mean))./Std(1:end-1)] > bThresh;
    allbCP=[allbCP bCP];
    
    TA_bCP = nan(1,length(isCpTrial));
    TA_bCP(1) = 0;
    for i = 2:length(isCpTrial)
        if bCP(i)==1
            TA_bCP(i) = 0;
        else
            TA_bCP(i) = TA_bCP(i-1)+1;
        end
    end
    allTA_bCP = [allTA_bCP TA_bCP+1];
    
    %******************************************************
    
    
    % running total of number of trials after change point, where cp = 0
    TAC = nan(1,length(isCpTrial));
    TAC(1) = 0;
    for i = 2:length(isCpTrial)
        if isCpTrial(i)==1
            TAC(i) = 0;
        else
            TAC(i) = TAC(i-1)+1;
        end
    end
    allTAC = [allTAC TAC];
    
    
    %Calculate the wagerTAC
    wagerTAC1 = TAC(wagerTrialIdx)+1;
    
    
    
    
    
    %***** Correcting perceptual and prediction error for spatial effects***********
    
    binLabel = -1*ones(length(prediction),1);
    predErr = wagerOutcomes - prediction;
    percpErr = wagerOutcomes - estimate;
    %predErr = wagerMeans - prediction;
    %percpErr = wagerMeans - estimate;
    auxNormPercpErr = percpErr;
    auxNormPredErr = predErr;
    
    %DEBUG -- normalize wrt to total std of pred and percep err
    predIdxNonCP = wagerTAC ~= 1;
    stdRawPredErr = std(predErr(predIdxNonCP));
    
    %DEBUG --- MAKE SURE TO CHECK!! (include CP in std of predErr)
    %stdRawPredErr = std(predErr);
    
    stdRawPercpErr = std(percpErr);
    
    %Discard outlier points with very large error (only for spatial bias
    %analysis)  -- these are usually points where subject wasn't paying
    %attention
    predErrLimit = 60;
    percpErrLimit = 50;
    predErr(abs(predErr)> predErrLimit & predIdxNonCP) = nan;
    percpErr(abs(percpErr)>percpErrLimit) = nan;
    
    % DEBUG -- DO YOU WANT TO INCLUDE OUTLIER PRED ERRORS IN THE FITS?
    %auxNormPercpErr = percpErr;
    %auxNormPredErr = predErr;
    
    
    fprintf('Subj %d Discarding %d pred points for angular correction\n',ID(i1),...
        sum(abs(predErr)> predErrLimit & predIdxNonCP));
    fprintf('Subj %d Discarding %d percept points for angular correction\n',ID(i1),...
        sum(abs(percpErr)>percpErrLimit));
    
    for i2 = 1:length(angleBins)-1
        binIdx = wagerOutcomes >= angleBins(i2) & wagerOutcomes <= angleBins(i2+1);
        spatialPercpErrBias(i1,i2) = nanmedian(percpErr(binIdx));
        spatialPercpErrStd(i1,i2) = nanstd( percpErr(binIdx));
        binIdxPred = wagerOutcomes >= angleBins(i2) & ...
            wagerOutcomes <= angleBins(i2+1) & wagerTAC ~= 1;
        spatialPredErrBias(i1,i2) = nanmedian(predErr(binIdxPred));
        spatialPredErrStd(i1,i2) = nanstd( predErr(binIdxPred));
        binLabel(binIdx) = i2;
        % Should we equalize the spatial variance?
        auxNormPercpErr(binIdx) = ((auxNormPercpErr(binIdx) - spatialPercpErrBias(i1,i2))/...
            spatialPercpErrStd(i1,i2));
        auxNormPredErr(binIdx) = ((auxNormPredErr(binIdx) - spatialPredErrBias(i1,i2))/...
            spatialPredErrStd(i1,i2));
        %auxNormPercpErr(binIdx) = (auxNormPercpErr(binIdx) - spatialPercpErrBias(i1,i2));
        %auxNormPredErr(binIdx) = (auxNormPredErr(binIdx) - spatialPredErrBias(i1,i2));
    end
    
    %Normalize errors
    auxNormPercpErr = auxNormPercpErr*stdRawPercpErr;
    auxNormPredErr = auxNormPredErr*stdRawPredErr;
    
    %DEBUG
    %auxNormPredErr = auxNormPredErr*stdRawPercpErr;
    
    
    fprintf('std percp err ---> %f\n',stdRawPercpErr);
    fprintf('std pred err ---> %f\n',stdRawPredErr);
    
    
    %DEBUG ---- MAKE SURE TO CHECK (BYPASS BINNING)
    %auxNormPercpErr = percpErr;
    %auxNormPredErr = predErr;
    
    spatialNormPercpErr = [spatialNormPercpErr auxNormPercpErr];
    spatialNormPredErr = [spatialNormPredErr auxNormPredErr];
    
    
    
    % ------- some additional containers required for model fitting
    % (Kamesh)
    
    %which of the wager trials (probes) were proper fixations?
    isFixed = aux.isFixed;
    
    
    
    
    % get block and trial numbers:  %hard coding
    % NOTE : don't hardcode numSoundsPerBlock if you want to use Day1 data
    
    if(length(outcome)>2200)
        numSoundsPerBlock = 600;
    elseif(length(outcome)<2050)
        numSoundsPerBlock = 500;
    else
        error('something fishy in numSoundsPerBlock!!!')
    end
    
    newBlock=1:numSoundsPerBlock:length(outcome);
    
    trialNum=nan(length(outcome), 1);
    trialNum(1)=1;
    blkNum=nan(length(outcome), 1);
    blkNum(1)=1;
    for i = 2:length(Std)
        if ~any(i==newBlock);
            blkNum(i)=blkNum(i-1);
        else
            blkNum(i)=blkNum(i-1)+1;
        end
        
        if  ~any(i-1==wagerTrialIdx);
            trialNum(i)=trialNum(i-1);
        else
            trialNum(i)=trialNum(i-1)+1;
        end
    end
    
    soundBlkNum = blkNum;   %for each sound
    auxProbeBlkNum = blkNum(wagerTrialIdx); %for each probe
    
    isProbe = false(size(trialNum));
    isProbe(wagerTrialIdx) = true;
    
    allTrialNum = [allTrialNum; trialNum];
    allIsProbe = logical([allIsProbe; isProbe]);
    allSoundBlkNum = [allSoundBlkNum; soundBlkNum];
    allSessID = [allSessID;  i1*ones(length(soundBlkNum),1)];
    
    probeBlkNum = [probeBlkNum; auxProbeBlkNum];
    
    
    
    
    
    % -------------------------- TOBII PUPIL FILTERING ---------------
    
    
    %2 - remove the first wager trial of each block
    wagerStd = Std(wagerTrialIdx);
    %     blkStart = find(diff([1 wagerStd]));
    %     blkEnd = [find(diff(wagerStd)) length(wagerStd)];
    %     isFixed([blkStart]) = 0;
    
    fP = aux.fP;         %get the false positive trial indices
    isFixed(fP) = 0;
    numFixed = sum(isFixed==1);
    idxFixed = find(isFixed);
    
    probeIsFixed=[probeIsFixed, isFixed];
    
    
    
    fixTimes = aux.pupTimes;
    outcomeTimes = aux.outcomeTimes;
    FilteredAlignedArray = nan(numFixed,numPupSamples);
    
    %     % average pupil dimensions over all the trials
    avgPupL = zeros(numPupSamples,1);
    avgPupR = zeros(numPupSamples,1);
    
    for i2 = 1:numFixed
        aux1 = aux.pupilL{idxFixed(i2)};
        aux2 = aux.pupilR{idxFixed(i2)};
        %aux1(1:5) = [];
        %aux2(1:5) = [];
        aux1Interp = aux1;
        aux2Interp = aux2;
        
        %check for blinks -- missing pupil data longer than 5 samples
        
        missIdxL = find(aux1<0);
        missIdxR = find(aux2<0);
        
        if(~isempty(missIdxL))
            aux11L = [1 find(diff(missIdxL)>1)'+1];
            startIdxL = missIdxL(aux11L);
            aux12L = [1 find(diff(missIdxL)>1)'];
            endIdxL = [missIdxL(aux12L(2:end)) ; missIdxL(end)];
            
            maxMissLenL = max(endIdxL - startIdxL);
            
            if(maxMissLenL > 10)
                %isFixedTrial(i1) = 0;
                %isBrokenTrial(i1) = 1;
                fprintf('Broken fixation found !! \n');
            else
                %now interpolate the blinks
                for i3 = 1:length(startIdxL)
                    if(endIdxL(i3)~=length(aux1) && startIdxL(i3)~=1)
                        preAux = aux1(startIdxL(i3)-1);
                        postAux = aux1(endIdxL(i3)+1);
                        interpAux = linspace(preAux,postAux,...
                            endIdxL(i3)-startIdxL(i3)+3);
                        aux1Interp(startIdxL(i3)-1:endIdxL(i3)+1) = interpAux;
                    else
                        if(endIdxL(i3)==length(aux1))
                            aux1Interp(startIdxL(i3)-1:end) = aux1(startIdxL(i3)-1);
                        end
                        if(startIdxL(i3)==1)
                            aux1Interp(startIdxL(i3):endIdxL(i3)+1) = aux1(endIdxL(i3)+1);
                        end
                    end
                end
            end
        end
        
        
        
        if(~isempty(missIdxR))
            aux11R = [1 find(diff(missIdxR)>1)'+1];
            startIdxR = missIdxR(aux11R);
            aux12R = [1 find(diff(missIdxR)>1)'];
            endIdxR = [missIdxR(aux12R(2:end)) ; missIdxR(end)];
            
            maxMissLenR = max(endIdxR - startIdxR);
            
            if(maxMissLenR > 10)
                %isFixedTrial(i1) = 0;
                %isBrokenTrial(i1) = 1;
                fprintf('Broken fixation found !! \n');
            else
                %now interpolate the blinks
                
                for i3 = 1:length(startIdxR)
                    if(endIdxR(i3)~=length(aux2) && startIdxR(i3)~=1)
                        preAux = aux2(startIdxR(i3)-1);
                        postAux = aux2(endIdxR(i3)+1);
                        interpAux = linspace(preAux,postAux,...
                            endIdxR(i3)-startIdxR(i3)+3);
                        aux2Interp(startIdxR(i3)-1:endIdxR(i3)+1) = interpAux;
                    else
                        if(endIdxR(i3)==length(aux2))
                            aux2Interp(startIdxR(i3)-1:end) = aux2(startIdxR(i3)-1);
                        end
                        
                        if(startIdxR(i3)==1)
                            aux2Interp(startIdxR(i3):endIdxR(i3)+1) = aux2(endIdxR(i3)+1);
                        end
                    end
                end
                
            end
        end
        
        
        
        %Filter the interpolated traces
        myFilt1.reset;
        prePad = ones(100,1)*aux1Interp(1);
        postPad = ones(100,1)*aux1Interp(end);
        auxFilt = [prePad; aux1Interp; postPad];     % DEBGU CHECK !!!!
        auxFilt = filtfilt(myFilt1.sosMatrix,myFilt1.ScaleValues,auxFilt);
        filtPupL = auxFilt(101:end-100);
        %filtPupL{i2} = filtfilt(myFilt1.sosMatrix,myFilt1.ScaleValues,aux1Interp);
        
        
        %Right Eye
        %myFilt = design(d);
        myFilt1.reset;
        prePad = ones(100,1)*aux2Interp(1);
        postPad = ones(100,1)*aux2Interp(end);
        auxFilt = [prePad; aux2Interp; postPad];     % DEBGU CHECK !!!!
        auxFilt = filtfilt(myFilt1.sosMatrix,myFilt1.ScaleValues,auxFilt);
        filtPupR = auxFilt(101:end-100);
        
        
        FilteredAlignedArray(i2,:) = (filtPupR(1:numPupSamples)+filtPupL(1:numPupSamples))*0.5;
        
    end
    
    %record the raw filtered pupil before doing anything
    auxRawPupil = FilteredAlignedArray;
    fixRawPupil = [fixRawPupil; FilteredAlignedArray];
    fixRawBaseline = [fixRawBaseline; mean(FilteredAlignedArray(:,1:3),2)];
    
    
    
    
    
    
    %DEBUG -- remove the linear trend in Baseline for each day
    auxSessionBase = mean(FilteredAlignedArray(:,1:3),2);
    X0 = ones(length(auxSessionBase),1);
    X1 =  aux.wagerOutcomeTimes';
    X1 = X1(isFixed);
    X1 = X1 - nanmean(X1);
    X = [X0 X1];
    auxSessionBaseAvg = mean(auxSessionBase);
    y = auxSessionBase - auxSessionBaseAvg;
    betaBaseTrend = robustfit(X,y,'bisquare',[],'off');
    trendBaseline = X*betaBaseTrend + auxSessionBaseAvg;
    auxTrendSubPupil = bsxfun(@minus,FilteredAlignedArray,trendBaseline);
    fixRawTrendSubPupil = [fixRawTrendSubPupil; auxTrendSubPupil];
    
    
    
    %DEBUG --- variables we need for doing the whitening in the main
    %regression
    prevTrialPupilRaw = zeros(size(FilteredAlignedArray));
    nextTrialPupilRaw = zeros(size(FilteredAlignedArray));
    prevTrialPupilRawTrendSub = zeros(size(FilteredAlignedArray));
    nextTrialPupilRawTrendSub = zeros(size(FilteredAlignedArray));
    nextTrialPupilMeanSub = zeros(size(FilteredAlignedArray));
    prevTrialPupilMeanSub = zeros(size(FilteredAlignedArray));
    prevTrialPupilMSZS = zeros(size(FilteredAlignedArray));
    nextTrialPupilMSZS = zeros(size(FilteredAlignedArray));
    auxAvgPupil = mean(FilteredAlignedArray,1);
    
    
    
    %2 - normalize each block (zscoring)
    %wagerStd = Std(wagerTrialIdx);
    %blkStart = find(diff([1 wagerStd(isFixed)]));
    %blkEnd = [find(diff(wagerStd(isFixed))) length(wagerStd(isFixed))];
    blkStart = zeros(4,1);
    blkEnd = zeros(4,1);
    blkStart(1) = 1;
    blkStart(2) = find(wagerTrialIdx(isFixed)>numSoundsPerBlock,1);
    blkStart(3) = find(wagerTrialIdx(isFixed)>2*numSoundsPerBlock,1);
    blkStart(4) = find(wagerTrialIdx(isFixed)>3*numSoundsPerBlock,1);
    
    blkEnd(4) = length(wagerTrialIdx(isFixed));
    blkEnd(1) = find(wagerTrialIdx(isFixed)>numSoundsPerBlock,1)-1;
    blkEnd(2) = find(wagerTrialIdx(isFixed)>2*numSoundsPerBlock,1)-1;
    blkEnd(3) = find(wagerTrialIdx(isFixed)>3*numSoundsPerBlock,1)-1;
    
    
    FilteredAlignedArrayMeanSub = FilteredAlignedArray;
    
    %baseLineWhitenedPupil = zeros(size(FilteredAlignedArray));
    for i2 = 1:length(blkStart)
        
        
        %DEBUG -- mean-subtract before z-scoring
        auxBlockMean = mean(FilteredAlignedArrayMeanSub(blkStart(i2):blkEnd(i2),:),1);
        aux2MS = bsxfun(@minus,FilteredAlignedArrayMeanSub(blkStart(i2):blkEnd(i2),:),...
            auxBlockMean);
        
        %aux2MS = reshape(aux2MS',1,size(aux2MS,1)*size(aux2MS,2));
        %mm1 = mean(aux2MS);
        %ss1 = std(aux2MS);
        %aux2MS = (aux2MS-mm1)/ss1;
        %aux2MS = reshape(aux2MS',numPupSamples,blkEnd(i2)-blkStart(i2)+1);
        %aux2MS = aux2MS';
        FilteredAlignedArrayMeanSub(blkStart(i2):blkEnd(i2),:) = aux2MS;
        
        %aux2 = FilteredAlignedArray(blkStart(i2):blkEnd(i2),:);
        aux2 = aux2MS;
        aux2 = reshape(aux2',1,size(aux2,1)*size(aux2,2));
        mm = mean(aux2);
        ss = std(aux2);
        aux2 = (aux2-mm)/ss;
        aux2 = reshape(aux2',numPupSamples,blkEnd(i2)-blkStart(i2)+1);
        aux2 = aux2';
        %aux2 = zscore(aux2,0,2);
        FilteredAlignedArray(blkStart(i2):blkEnd(i2),:) = aux2;
        
        
        aux2 = auxRawPupil(blkStart(i2):blkEnd(i2),:);
        %prevTrialPupil(blkStart(i2)+1:blkEnd(i2),:) = meanCentX(aux2(1:end-1,:));
        %prevTrialPupil(blkStart(i2)+1:blkEnd(i2),:) = bsxfun(@minus,aux2(1:end-1,:), ...
        %                                                     mean(aux2(1:end-1,:),2));
        prevTrialPupilRaw(blkStart(i2)+1:blkEnd(i2),:) = aux2(1:end-1,:);
        prevTrialPupilRaw(blkStart(i2),:) = zeros(1,size(aux2, 2));
        
        %nextTrialPupil(blkStart(i2)+1:blkEnd(i2),:) = bsxfun(@minus,aux2(2:end,:), ...
        %                                                     mean(aux2(2:end,:),2));
        nextTrialPupilRaw(blkStart(i2):blkEnd(i2)-1,:) = aux2(2:end,:);
        nextTrialPupilRaw(blkEnd(i2),:) = zeros(1,size(aux2,2));
        
        
        
        
        aux2 = auxTrendSubPupil(blkStart(i2):blkEnd(i2),:);
        prevTrialPupilRawTrendSub(blkStart(i2)+1:blkEnd(i2),:) = aux2(1:end-1,:);
        prevTrialPupilRawTrendSub(blkStart(i2),:) = zeros(1,size(aux2, 2));
        
        %nextTrialPupil(blkStart(i2)+1:blkEnd(i2),:) = bsxfun(@minus,aux2(2:end,:), ...
        %                                                     mean(aux2(2:end,:),2));
        nextTrialPupilRawTrendSub(blkStart(i2):blkEnd(i2)-1,:) = aux2(2:end,:);
        nextTrialPupilRawTrendSub(blkEnd(i2),:) = zeros(1,size(aux2,2));
        
        
        aux2 = aux2MS;
        prevTrialPupilMeanSub(blkStart(i2)+1:blkEnd(i2),:) = aux2(1:end-1,:);
        prevTrialPupilMeanSub(blkStart(i2),:) = zeros(1,size(aux2, 2));
        
        %nextTrialPupil(blkStart(i2)+1:blkEnd(i2),:) = bsxfun(@minus,aux2(2:end,:), ...
        %                                                     mean(aux2(2:end,:),2));
        nextTrialPupilMeanSub(blkStart(i2):blkEnd(i2)-1,:) = aux2(2:end,:);
        nextTrialPupilMeanSub(blkEnd(i2),:) = zeros(1,size(aux2,2));
        
        
        aux2 = FilteredAlignedArray(blkStart(i2):blkEnd(i2),:);
        prevTrialPupilMSZS(blkStart(i2)+1:blkEnd(i2),:) = aux2(1:end-1,:);
        prevTrialPupilMSZS(blkStart(i2),:) = zeros(1,size(aux2, 2));
        
        %nextTrialPupil(blkStart(i2)+1:blkEnd(i2),:) = bsxfun(@minus,aux2(2:end,:), ...
        %                                                     mean(aux2(2:end,:),2));
        nextTrialPupilMSZS(blkStart(i2):blkEnd(i2)-1,:) = aux2(2:end,:);
        nextTrialPupilMSZS(blkEnd(i2),:) = zeros(1,size(aux2,2));
        
        
        
        %run length STD calculation
        runLengthStd(blkStart(i2)) = nan;
        
    end
    
    fixRunLengthSTD = [fixRunLengthSTD; runLengthStd(isFixed)];
    
    fixPrevTrialRawPupil = [fixPrevTrialRawPupil; prevTrialPupilRaw];  %traces from trial before and after
    fixNextTrialRawPupil = [fixNextTrialRawPupil; nextTrialPupilRaw];
    fixPrevTrialRawTrendSubPupil = [fixPrevTrialRawTrendSubPupil; prevTrialPupilRawTrendSub ];  %traces from trial before and after
    fixNextTrialRawTrendSubPupil = [fixNextTrialRawTrendSubPupil; nextTrialPupilRawTrendSub ];
    fixPrevTrialMeanSubPupil = [fixPrevTrialMeanSubPupil; prevTrialPupilMeanSub];  %traces from trial before and after
    fixNextTrialMeanSubPupil = [fixNextTrialMeanSubPupil; nextTrialPupilMeanSub];
    fixPrevTrialMSZSPupil = [fixPrevTrialMSZSPupil; prevTrialPupilMSZS];
    fixNextTrialMSZSPupil = [fixNextTrialMSZSPupil; nextTrialPupilMSZS];
    
    
    %Josh's suggestion of removing mean trend, then smoothing
    fixMeanSubtrPupil = [fixMeanSubtrPupil; FilteredAlignedArrayMeanSub];
    %fixMsBaseline = [fixMsBaseline; mean(auxMeanSubtrPup(:,1:3),2)];
    
    fixFilteredPupil = [fixFilteredPupil; FilteredAlignedArray];
    
    
    
    %RUN LENGTH RELATED CALCULATIONS
    
    %PUPIL AS A FUNCTION OF RUN LENGTH
    % the difference from above is taht we include broken fixation trials
    
    
    blkStart = zeros(4,1);
    blkEnd = zeros(4,1);
    blkStart(1) = 1;
    blkStart(2) = find(wagerTrialIdx>numSoundsPerBlock,1);
    blkStart(3) = find(wagerTrialIdx>2*numSoundsPerBlock,1);
    blkStart(4) = find(wagerTrialIdx>3*numSoundsPerBlock,1);
    
    blkEnd(4) = length(wagerTrialIdx);
    blkEnd(1) = find(wagerTrialIdx>numSoundsPerBlock,1)-1;
    blkEnd(2) = find(wagerTrialIdx>2*numSoundsPerBlock,1)-1;
    blkEnd(3) = find(wagerTrialIdx>3*numSoundsPerBlock,1)-1;
    
    
    auxFixTrialNo = 1:blkEnd(4);
    auxFixTrialNo(blkStart(2):blkEnd(2)) = auxFixTrialNo(blkStart(2):blkEnd(2)) - ...
        blkEnd(1);
    auxFixTrialNo(blkStart(3):blkEnd(3)) = auxFixTrialNo(blkStart(3):blkEnd(3)) - ...
        blkEnd(2);
    auxFixTrialNo(blkStart(4):blkEnd(4)) = auxFixTrialNo(blkStart(4):blkEnd(4)) - ...
        blkEnd(3);
    
    
    fixTrialNo = [fixTrialNo auxFixTrialNo(isFixed)];
    
    avgRLperTAC = zeros(numTAC,1);
    stdRLperTAC = zeros(numTAC,1);
    medianRLperTAC = zeros(numTAC,1);
    
    %Calculate the average run-length as a function of TAC
    runLengths = -1*ones(size(wagerTAC));
    
    auxWagerTAC = wagerTAC;
    auxWagerTAC(blkStart) = -1;   %Don't count runLengths for block starts
    for i3 = 1:numTAC
        tacWagerIdx = find(auxWagerTAC==i3);
        tacWagerTrial = wagerTrialIdx(tacWagerIdx);
        tacPrevWagerTrial = wagerTrialIdx(max(tacWagerIdx-1,1));
        avgRLperTAC(i3) = mean(tacWagerTrial - tacPrevWagerTrial);
        stdRLperTAC(i3) = std(tacWagerTrial - tacPrevWagerTrial);
        medianRLperTAC(i3) = median(tacWagerTrial - tacPrevWagerTrial);
        runLengths(tacWagerIdx) = tacWagerTrial - tacPrevWagerTrial;
    end
    fixRunLengths = [fixRunLengths runLengths(isFixed)];
    
    
    %We need to do this here because some block starts have very large
    %tsince last probe values
    timeSinceLastProbe = aux.timeSinceLastProbe;
    timeSinceLastProbe(blkStart) = 2;
    
    
    
    
    % concatenate subject data
    fixOutcomes = [fixOutcomes wagerOutcomes(isFixed)];
    if(i1 > 106)
        fixDay1Sess = [fixDay1Sess ; 1*ones(length(wagerOutcomes(isFixed)),1)];
    else
        fixDay1Sess = [fixDay1Sess ; 0*ones(length(wagerOutcomes(isFixed)),1)];
    end
    probeOutcomes = [probeOutcomes wagerOutcomes];
    probeCentralOutcomes = [probeCentralOutcomes (wagerOutcomes > 20 & wagerOutcomes < 160)];
    if(i1 > 106)
        probeDay1Sess = [probeDay1Sess ; 1*ones(length(wagerOutcomes),1)];
    else
        probeDay1Sess = [probeDay1Sess ; 0*ones(length(wagerOutcomes),1)];
    end
    
    fixSubjID  = [fixSubjID ID(i1)*ones(1,sum(isFixed))];
    probeSubjID  = [probeSubjID ID(i1)*ones(1,length(wagerMeans))];
    fixMean = [fixMean wagerMeans(isFixed)];
    probeMean = [probeMean wagerMeans];
    fixStd = [fixStd Std(wagerTrialIdx(isFixed))];
    probeStd = [probeStd Std(wagerTrialIdx)];
    fixTAC = [fixTAC wagerTAC(isFixed)];
    probeTAC = [probeTAC wagerTAC];
    fixbCPTAC = [fixbCPTAC TA_bCP(wagerTrialIdx(isFixed))+1];
    probebCPTAC = [probebCPTAC TA_bCP(wagerTrialIdx)+1];
    fixBaseline = [fixBaseline; mean(FilteredAlignedArray(:,1:3),2)];
    fixIsFixed = [fixIsFixed isFixed];
    fixPrediction = [fixPrediction prediction(isFixed)];
    probePrediction = [probePrediction prediction];
    fixEstimate = [fixEstimate estimate(isFixed)];
    probeEstimate = [probeEstimate estimate];
    
    fixNormPredErr = [fixNormPredErr auxNormPredErr(isFixed)];
    fixNormPercpErr = [fixNormPercpErr auxNormPercpErr(isFixed)];
    
    fixWager = [fixWager wager(isFixed)];
    probeWager = [probeWager wager];
    
    bettingFrac = sum(wager(isFixed))/length(wager(isFixed));
    if ( bettingFrac > 0.15 & bettingFrac < 0.85)
        auxExt = zeros(size(wager(isFixed)));
    else
        auxExt = ones(size(wager(isFixed)));
    end
    fixExtremeWager = [fixExtremeWager auxExt];
    
    auxLastOutcomeTAC = TAC(wagerTrialIdx(isFixed)-1)+1;
    fixLastOutcomeTAC = [fixLastOutcomeTAC auxLastOutcomeTAC];
    
    % Session ID
    fixSessionID = [fixSessionID; i1*ones(length(wager(isFixed)),1)];
    probeSessionID = [probeSessionID; i1*ones(length(wager),1)];
    
    
    % Timing related variables ***************************************
    fixTsinceLastProbe = [fixTsinceLastProbe timeSinceLastProbe(isFixed)];
    disp(['Subj -> ' num2str(i1)]);
    [avgRLperTAC stdRLperTAC medianRLperTAC];
    wagerOutcomeTimes = aux.wagerOutcomeTimes;
    fixWagerOutcomeTimes = [fixWagerOutcomeTimes wagerOutcomeTimes(isFixed)];
    auxTimeToPredict = aux.timeToPredict;
    fixTimeToPredict = [fixTimeToPredict auxTimeToPredict(isFixed)];
    auxTimeToEstimate = aux.timeInEstimate;
    fixTimeToEstimate = [fixTimeToEstimate ; auxTimeToEstimate(isFixed)];
    
    aux_1 = aux.outcomeTimes(wagerTrialIdx(1:end-1)+1);
    aux_2 = aux.wagerOutcomeTimes(1:end-1) + 2.5 + auxTimeToEstimate(1:end-1)';
    auxTimeToBet = aux_1 - aux_2;
    fixTimeToBet = [fixTimeToBet; [auxTimeToBet(isFixed(1:end-1))' ; -1] ];
    
    aux_1 = aux.outcomeTimes(wagerTrialIdx(1:end-1)+1);
    aux_2 = aux.wagerOutcomeTimes(2:end) - auxTimeToPredict(2:end);
    fixInterProbeInterval = [fixInterProbeInterval; aux_2' - aux_1'];
    
    %DEBUG
    if(length(probeSessionID) ~= length(probeSubjID))
        disp('Something fishy!!')
        keyboard
    end
end


% MRN
% CREATE PREDICTIONS FROM FRUGFUN5:
HazardRate=.15; driftRate=0; likelihoodWeight=1;

% create an array of new blocks:
newBlock=false(length(allStd), 1);
newBlock(1)=true;
for i = 2:length(allStd)
    if allStd(i)~=allStd(i-1)||allSessID(i)~=allSessID(i-1)
        newBlock(i)=true;
    end
end

startBlock=find(newBlock);
endBlock=[startBlock(2:end)-1; length(newBlock)];
modPred=nan(length(allStd), 1);

% Loop through blocks and run model.
for i = 1:length(startBlock)
    sel=startBlock(i):endBlock(i);
    noise=unique(allStd(sel));
    [B, totSig, R, pCha, sigmaU] = frugFun5_wRange(allOutcome(sel), ...
        HazardRate, noise, driftRate, likelihoodWeight, false,...
        150, 1, [0, 180]);
    modPred(sel)=B(1:end-1);
end


% IF THE forceRunModelSims flag is set to true model sims are regenerated 
if forceRunModelSims
    NoFitting_trueHazard_trainingAudWidth
    print(['Model sims written to ' modelData_dir]);
end






% Now, calculate the variance as a function of PE, to get the weights for
% weighted regression
% as a sanity check, you should see similar results with ROBUSTFIT


% MODEL SELECTION NOTE : the model with a column of ones and the spatial
% bias term raised to power 1 has the lowest BIC -- so we use that here
uniqueStds = [10,20];
uniqueIDs = unique(ID);

predErr=probeOutcomes'-probePrediction';
percErr=probeOutcomes'-probeEstimate';
%predErr = fixNormPredErr';
%percErr = fixNormPercpErr';

absPE=abs(predErr);
abPEbins=prctile(absPE, 0:10:100);

priorWtVIF = [];

clear betas betaInt binStd peBinStd priorWtVIF
for i0 = 1:length(uniqueIDs)
    selID = fixSubjID == uniqueIDs(i0); %  & fixbCP < 7;  %the last part is for big change-points
    selIdx = selID;
    selSTD = fixStd(selIdx)';
    
    perceptualErr = percErr(selIdx)-nanmean(percErr(selIdx));
    predictionErr = predErr(selIdx) -nanmean(predErr(selIdx));
    
    %create the decision matrix
    X0 = ones(length(predictionErr),1);
    
    % columns corresponding to b1 - b6
    Xpred10 = zeros(length(predictionErr),numTAC);
    for i4 = 1:numTAC
        Xpred10(fixTAC(selIdx) == i4, i4) = 1;
    end
    Xpred10 = bsxfun(@times,Xpred10,predictionErr.*(selSTD==10));
    Xpred10 = bsxfun(@minus,Xpred10,nanmean(Xpred10,1));
    
    Xpred20 = zeros(length(predictionErr),numTAC);
    for i4 = 1:numTAC
        Xpred20(fixTAC(selIdx) == i4, i4) = 1;
    end
    Xpred20 = bsxfun(@times,Xpred20,predictionErr.*(selSTD==20));
    Xpred20 = bsxfun(@minus,Xpred20,nanmean(Xpred20,1));
    
    
    
    XtimePred10 = fixTimeToPredict(selIdx)'.*(selSTD==10);
    XtimePred10 = (XtimePred10 - mean(XtimePred10))./std(XtimePred10);
    
    XtimePred20 = fixTimeToPredict(selIdx)'.*(selSTD==20);
    XtimePred20 = (XtimePred20 - mean(XtimePred20))./std(XtimePred20);
    
    XtLastProbe10 = fixTsinceLastProbe(selIdx)'.*(selSTD==10);
    XtLastProbe10 = (XtLastProbe10 - mean(XtLastProbe10))/std(XtLastProbe10);
    
    XtLastProbe20 = fixTsinceLastProbe(selIdx)'.*(selSTD==20);
    XtLastProbe20 = (XtLastProbe20 - mean(XtLastProbe20))/std(XtLastProbe20);
    
    rawPred = fixPrediction(selIdx)';
    
    XSpace10 = (90 - rawPred).*(selSTD==10);
    XSpace10 = (XSpace10 - nanmean(XSpace10))./std(XSpace10);
    
    XSpace20 = (90 - rawPred).*(selSTD==20);
    XSpace20 = (XSpace20 - nanmean(XSpace20))./std(XSpace20);
    
    X = [X0 (Xpred10+Xpred20)  (XSpace10+XSpace20) ...
        (XtimePred10+XtimePred20) (XtLastProbe10+XtLastProbe20)];
    
    [B, BINT, res] = regress(perceptualErr, X);
    %[B stats] = robustfit(X,perceptualErr','bisquare',[],'off');
    
    betas(:,i0) =  B;
    
    
    %MRN: compute residuals for different absolute PE magnitudes
    for i = 2:length(abPEbins)
        peSel=absPE(selIdx)>abPEbins(i-1) & absPE(selIdx)<=abPEbins(i);
        peBinStd(i0,i-1)=nanstd(res(peSel));
    end
    
    
    VIF=nan(size(X, 2), 1);
    for i = 2:size(X, 2)        % don't calculate VIF for column of ones
        xy = zscore(X(:,i));
        xx=X(:, [1:i-1 i+1:size(X, 2)]);
        xx(:, 2:end)=zscore(xx(:, 2:end));
        [~,~,R,~,STATS] = regress(xy,xx);
        VIF(i)=1./(1-STATS(1));
    end
    priorWtVIF(i0,:) = VIF;
    
end


% Correcting for heteroscedasticity : Individual
% subject measurements seem noisy, so we'll use the average standard
% deviation for each subject to decide how much to weigh each error bin.

% And we can create an array of standard deviations


trialExpStd=nan(size(absPE));
for i = 1:length(abPEbins)-1
    trialExpStd(absPE>=abPEbins(i))=nanmean(peBinStd(:,i));
    allBinStd(i)=nanmean(peBinStd(:,i));
end

fixTrialExpStd = trialExpStd(logical(fixIsFixed));

% Plotting initialisation
getCbColors
discCols = lines(8);





% Calculate the STD and MAD for PE and percpErr

subjPercpMAD = nan(size(uniqueIDs));
subjPredMAD = nan(size(uniqueIDs));

subj_estErr_std = nan(size(uniqueIDs));
subj_predErr_std = nan(size(uniqueIDs));
subj_normPredErr_std = nan(size(uniqueIDs));  %normalise error by noise condition
subj_normEstErr_std = nan(size(uniqueIDs));   %normalise error by noise condition

subj_estErr_std_hi = nan(size(uniqueIDs));
subj_estErr_std_lo = nan(size(uniqueIDs));
subj_predErr_std_hi = nan(size(uniqueIDs));
subj_predErr_std_lo = nan(size(uniqueIDs));

subjPercpMAD_hi = nan(size(uniqueIDs));
subjPredMAD_hi = nan(size(uniqueIDs));

subjPercpMAD_lo = nan(size(uniqueIDs));
subjPredMAD_lo = nan(size(uniqueIDs));

for i1 = 1:length(uniqueID)
    
    selSubj = probeSubjID == i1;
    selSubj_pred = probeSubjID == i1 & probeTAC ~= 1; %discard CP trials for pred error
    selSubj_hi = probeSubjID == i1 & probeStd == 20;
    selSubj_lo = probeSubjID == i1 & probeStd == 10;
    
    selSubj_pred_hi = probeSubjID == i1 & probeStd == 20 & probeTAC ~= 1;
    selSubj_pred_lo = probeSubjID == i1 & probeStd == 10 & probeTAC ~= 1;
    
    
    subj_STD = probeStd(selSubj);
    subj_STD_pred = probeStd(selSubj_pred);
    
    subjPredErr = probeOutcomes(selSubj_pred)-probePrediction(selSubj_pred);
    auxMedPred = nanmedian(subjPredErr);
    subjPredMAD(i1) = nanmedian(abs(subjPredErr - auxMedPred));
    
    subjPercpErr = probeOutcomes(selSubj)-probeEstimate(selSubj);
    auxMedPercp = nanmedian(subjPercpErr);
    subjPercpMAD(i1) = nanmedian(abs(subjPercpErr - auxMedPercp));
    
    
    subj_estErr_std(i1) = std(subjPercpErr);
    subj_predErr_std(i1) = std(subjPredErr);
    subj_normPredErr_std(i1) = std(subjPredErr./subj_STD_pred);
    subj_normEstErr_std(i1) = std(subjPercpErr./subj_STD);
    
    
    subjPredErr_hi = probeOutcomes(selSubj_pred_hi)-probePrediction(selSubj_pred_hi);
    auxMedPred = nanmedian(subjPredErr_hi);
    subjPredMAD_hi(i1) = nanmedian(abs(subjPredErr_hi - auxMedPred));
    subj_predErr_std_hi(i1) = std(subjPredErr_hi);
    
    subjPercpErr_hi = probeOutcomes(selSubj_hi)-probeEstimate(selSubj_hi);
    auxMedPercp = nanmedian(subjPercpErr_hi);
    subjPercpMAD_hi(i1) = nanmedian(abs(subjPercpErr_hi - auxMedPercp));
    subj_estErr_std_hi(i1) = std(subjPercpErr_hi);
    
    subjPredErr_lo = probeOutcomes(selSubj_lo)-probePrediction(selSubj_lo);
    auxMedPred = nanmedian(subjPredErr_lo);
    subjPredMAD_lo(i1) = nanmedian(abs(subjPredErr_lo - auxMedPred));
    subj_predErr_std_hi(i1) = std(subjPredErr_lo);
    
    subjPercpErr_lo = probeOutcomes(selSubj_lo)-probeEstimate(selSubj_lo);
    auxMedPercp = nanmedian(subjPercpErr_lo);
    subjPercpMAD_lo(i1) = nanmedian(abs(subjPercpErr_lo - auxMedPercp));
    subj_estErr_std_lo(i1) = std(subjPercpErr_lo);
end



% Trraining stuff
std_subj_reports_training = [];
training_angles = unique(training_all_outcome);

for i0 = 1:length(unique(ID))
    
    sel_subj = training_all_ID == i0;
    
    aux_reports = training_all_estimate(sel_subj);
    aux_outcomes = training_all_outcome(sel_subj);
    std_subj_reports_training(i0) = std(aux_reports-aux_outcomes);
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               LOAD THE DATA NEEDED FROM NORMATIVE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~forceRunModelSims % if you've run the model sims, you should already have this stuff.
    load(modelData_file);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               FLAG THAT CONTROLS WHETHER OR NOT WE USE DAY 1 DATA FOR
%               ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


use_day1 = false;



%% *****************************************************************************
%                           TESTING VERSIONS OF THE LINEAR MODEL
%
% ******************************************************************************


% This file analyses various versions of the linear model and serves as a
% reference for their performance. Note, there is always one intercept term

%    v1 -- both noise conditions collapsed *without* a PE*Stab*RU term
%    v1a -- same as v1 but with a PE*noise term



%% ----------MODEL-BASED REGRESSION FOR THE DIFFERENT VERSIONS


useTotalUnc = true;  %if true then we don't mean center per noise block
uniqueStds = [10,20];
uniqueIDs = unique(ID);

calcVIF = true;

linMod_AIC  = [];
linMod_BIC = [];

linMod_v1_beta_PW = zeros(size(probeSubjID));
linMod_v1_sim_perceptual_err = zeros(size(probeSubjID));
priorWtInfResid_v1 = nan(size(probeSubjID))';


absPE=abs(probeOutcomes'-probePrediction');
abPEbins=prctile(absPE, 0:10:100);



% Load simulated data to get benchmark values for regression coefficients.
if ~forceRunModelSims
    load(modelData_file)
end

% probeModEstSim = simulated model perceptions.
simOptPercError=repmat(probeOutcomes', 1, size(probeModEstSim,2)) - probeModEstSim;


% calculate prior weights vs TAC for the model
clear   beta_v1  beta_v1a beta_v1b
clear   subj_VIF_v1  subj_VIF_v1a subj_VIF_v1b peBinStd peBinStd_model



varianceComp=nan(size(probeEstimate));
% biasComp(selID)=linMod_v1_sim_perceptual_err;
justBiasComp=nan(size(probeEstimate));
biasAndSpace=nan(size(probeEstimate));


% load ~/Desktop/MRN_tmp_workspace_8-4-16


beta_v1b_sim=nan(length(uniqueIDs), 9, size(simOptPercError, 2));


for i0 = 1:length(uniqueIDs)
    
    clear aic bic
    
    if(use_day1)
        selID = (probeSubjID == uniqueIDs(i0))';
    else
        selID = (probeSubjID == uniqueIDs(i0))' & probeDay1Sess == 0;   % -- nonCP trials :  & probeTAC ~= 1;
    end
    
    subjBets = probeWager(selID)';
    betFrac = sum(subjBets)/length(subjBets);
    if(betFrac > 0.1 & betFrac < 0.9)
        subjNormalBets(i0) = true;
    else
        subjNormalBets(i0) = false;
    end
    
    
    
    % PREDICTION ERROR
    mcPredErr = probeOutcomes(selID)' - probePrediction(selID)';
    mcPredErr = mcPredErr - nanmean(mcPredErr);
    
    mcPredErr_lo = (probeOutcomes(selID)' - probePrediction(selID)');
    mcPredErr_lo = mcPredErr_lo.*(probeStd(selID)' == 10);
    mcPredErr_lo = mcPredErr_lo - nanmean(mcPredErr_lo(mcPredErr_lo~=0));
    
    mcPredErr_hi = (probeOutcomes(selID)' - probePrediction(selID)');
    mcPredErr_hi = mcPredErr_hi.*(probeStd(selID)' == 20);
    mcPredErr_hi = mcPredErr_hi - nanmean(mcPredErr_hi(mcPredErr_hi~=0));
    
    % PERCEPTUAL ERROR
    mcPercpErr = probeOutcomes(selID)' - probeEstimate(selID)';
    mcPercpErr = mcPercpErr - nanmean(mcPercpErr);  % Do you want to m-c PcE?
    
    
    %Correcting perceptual error for bias (load the training data file
    %before this
    aux_estimate = probeEstimate(selID)';
    [~, aux_spatial_idx] = min(abs(bsxfun(@minus,probeOutcomes(selID)',mean_train_angles')),[],2);
    aux_estimate = aux_estimate + train_PE_subject_bias_binned(i0,aux_spatial_idx)';
    %mcPercpErr = probeOutcomes(selID)' - aux_estimate;
    %mcPercpErr = mcPercpErr - nanmean(mcPercpErr);
    mcSpatialBias = train_PE_subject_bias_binned(i0,aux_spatial_idx)';
    allTrainingSpatialBias(selID)=mcSpatialBias;
    mcSpatialBias = mcSpatialBias - nanmean(mcSpatialBias);
    
    % STABLITY TERM
    mcSubjCPterm = probeModelStabSubj(selID);
    mcSubjCPterm =  mcSubjCPterm - nanmean(mcSubjCPterm);
    
    mcSubjCPterm_lo = probeModelStabSubj(selID).*(probeStd(selID)' == 10) ;
    mcSubjCPterm_lo =  mcSubjCPterm_lo - nanmean(mcSubjCPterm_lo);
    
    mcSubjCPterm_hi = probeModelStabSubj(selID).*(probeStd(selID)' == 20) ;
    mcSubjCPterm_hi =  mcSubjCPterm_hi - nanmean(mcSubjCPterm_hi);

    % UNCERTAINTY TERM
    mcSubjUncTerm = probeModelUncSubj(selID); % This term is: audNoise.^2 ./ (totUncSubj.^2 + audNoise.^2);
    mcSubjUncTerm = mcSubjUncTerm - nanmean(mcSubjUncTerm);
    
    mcSubjUncTerm_lo = probeModelUncSubj(selID) .*(probeStd(selID)' == 10);
    ru_idx = mcSubjUncTerm_lo~=0;
    mcSubjUncTerm_lo(ru_idx) = mcSubjUncTerm_lo(ru_idx) - nanmean(mcSubjUncTerm_lo);
    
    mcSubjUncTerm_hi = probeModelUncSubj(selID) .*(probeStd(selID)' == 20);
    ru_idx = mcSubjUncTerm_hi~=0;
    mcSubjUncTerm_hi(ru_idx) = mcSubjUncTerm_hi(ru_idx) - nanmean(mcSubjUncTerm_hi);
    
    % mean-centered per noise condition
    mcSubjTrueUncTerm = mcSubjUncTerm_hi + mcSubjUncTerm_lo;
    
    if(useTotalUnc)
        mcSubjTrueUncTerm = mcSubjUncTerm;
    else
        mcSubjUncTerm = mcSubjTrueUncTerm;
    end
    
    
    mcWager = probeWager(selID)';
    mcWager = mcWager - nanmean(mcWager);
    
    % THIS is a signed error so stuff can get pulled toward center!!!
    mcSpaceTerm = (90 - probePrediction(selID)');
    mcSpaceTerm = mcSpaceTerm - nanmean(mcSpaceTerm);
    
    mcNoise = probeStd(selID)';
    mcNoise = mcNoise - nanmean(mcNoise);
    
    % THIS TERM IS ABSOLUTE!!!! so things near center can be less biased
    %modulation of prior usage depending on spatial location
    mcSpacePEmod = abs(90 - probePrediction(selID)');
    mcSpacePEmod = mcSpacePEmod - nanmean(mcSpacePEmod);
    
    Y = mcPercpErr;
    
    X1 = [ones(length(mcPredErr),1) mcPredErr  mcPredErr.*mcSubjCPterm ...
        mcPredErr.*mcSubjUncTerm, mcPredErr.*mcWager mcSpaceTerm mcSpatialBias];
    
    X1_a = [ones(length(mcPredErr),1) mcPredErr  mcPredErr.*mcSubjCPterm ...
        mcPredErr.*mcSubjUncTerm, mcPredErr.*mcWager mcNoise.*mcPredErr mcSpaceTerm mcSpatialBias];
    
    X1_b = [ones(length(mcPredErr),1) mcPredErr  mcPredErr.*mcSubjCPterm ...
        mcPredErr.*mcSubjUncTerm, mcPredErr.*mcWager mcNoise.*mcPredErr ...
        mcPredErr.*mcSpacePEmod mcSpaceTerm mcSpatialBias]; % mcSpaceTerm ];
    
    
    [b1,bint1,h1,p1,LL1] = regressW_mike(Y, trialExpStd(selID), X1);
    [b1a,bint1a,h1a,p1a,LL1a] = regressW_mike(Y, trialExpStd(selID), X1_a);
    [b1b,bint1b,h1b,p1b,LL1b] = regressW_mike(Y, trialExpStd(selID), X1_b);
    
    
    % MRN added 8/4/16 to check robustness of weighted regression
    % likelihoods against those computed using standard regression. 
    [~, ~,  R1]  = regress(Y,  X1);  
    [~, ~, R1a]  = regress(Y,  X1_a);
    [~, ~, R1b]  = regress(Y, X1_b);
    
    LL1_alt=sum(log(normpdf(R1, 0, nanstd(R1))));
    LL1a_alt=sum(log(normpdf(R1, 0, nanstd(R1a))));
    LL1b_alt=sum(log(normpdf(R1, 0, nanstd(R1b))));
   
    beta_v1(i0,:)  = b1;
    beta_v1a(i0,:) = b1a;
    beta_v1b(i0,:) = b1b;
    
    
    % MRN: compute regression coefficients for optimal generative model to
    % use as benchmarks
    for j = 1:size( simOptPercError,2)
        Y1 = simOptPercError(selID,j);
        [b1_sim,bint1_sim,h1_sim,p1_sim,LL1_sim] = regressW_mike(Y1, trialExpStd(selID), X1);
        beta_v1_sim(i0,:,j) = b1_sim;
    end
    
    
    %MRN: compute residuals for different absolute PE magnitudes
    res = (Y - X1*b1);
    for i = 2:length(abPEbins)
        peSel=absPE(selID)>abPEbins(i-1) & absPE(selID)<=abPEbins(i);
        peBinStd(i0,i-1)=nanstd(res(peSel));
    end
    
    % ------------- comaprison of linear model prior wt vs. behaviour -----
    %calculate pwt using beta coeffs estimated above
    linMod_v1_beta_PW(selID) = beta_v1(i0,2) + beta_v1(i0,3)*mcSubjCPterm + ...
        beta_v1(i0,4)*mcSubjUncTerm + beta_v1(i0,5)*mcWager;
    
    %calculate pwt using simulated percepts
    linMod_v1_sim_perceptual_err(selID) = X1*b1; % changed by MRN 7/11 to reflect winning model
    yHat_withoutSpatialTerms=X1_b(:,1:end-3)*b1b(1:end-3);
    
    
    % CALCULATE RESIDUALS FOR V1 model
    %Prior influence -- positive values indicate more PE influence than suspected
    R = (Y - X1*b1);
    priorWtInfResid_v1(selID)= R.*sign(mcPredErr);  %multiply by sign of raw pred err or mcPredErr
    
    % Bias/Variance decomposition:
    varianceComp(selID)=R;
   % biasComp(selID)=linMod_v1_sim_perceptual_err;
    justBiasComp(selID)= yHat_withoutSpatialTerms;
    biasAndSpace(selID)= linMod_v1_sim_perceptual_err(selID);
    

    %calculate AIC BIC
    numObs = length(Y);
    [aic(1),bic(1)] = aicbic(LL1,size(X1,2),numObs);
    [aic(2),bic(2)] = aicbic(LL1a,size(X1_a,2),numObs);
    [aic(3),bic(3)] = aicbic(LL1b,size(X1_b,2),numObs);
    
    % alternative AIC/BIC based on OLS regression:
    [aic_alt(1),bic_alt(1)] = aicbic(LL1_alt,size(X1,2),numObs);
    [aic_alt(2),bic_alt(2)] = aicbic(LL1a_alt,size(X1_a,2),numObs);
    [aic_alt(3),bic_alt(3)] = aicbic(LL1b_alt,size(X1_b,2),numObs);

    
    linMod_AIC(:,i0) =  aic;
    linMod_BIC(:,i0) = bic;
    
    linMod_AIC_alt(:,i0) =  aic_alt;
    linMod_BIC_alt(:,i0) = bic_alt;

    %Calculate the variance inflation factor for 1s < t < 1.5 s
    if(calcVIF)
        VIF_v1=nan(size(X1, 2), 1);
        VIF_v1a=nan(size(X1_a, 2), 1);
        VIF_v1b=nan(size(X1_b, 2), 1);
        
        for i = 2:size(X1, 2)        % don't calculate VIF for column of ones
            xy = zscore(X1(:,i));
            xx=X1(:, [1:i-1 i+1:size(X1, 2)]);
            xx(:, 2:end)=zscore(xx(:, 2:end));
            [~,~,R,~,STATS] = regress(xy,xx);
            VIF_v1(i)=1./(1-STATS(1));
        end
        subj_VIF_v1(i0,:) = VIF_v1;
        
        for i = 2:size(X1_a, 2)        % don't calculate VIF for column of ones
            xy = zscore(X1_a(:,i));
            xx=X1_a(:, [1:i-1 i+1:size(X1_a, 2)]);
            xx(:, 2:end)=zscore(xx(:, 2:end));
            [~,~,R,~,STATS] = regress(xy,xx);
            VIF_v1a(i)=1./(1-STATS(1));
        end
        
        subj_VIF_v1a(i0,:) = VIF_v1a;
        
        for i = 2:size(X1_b, 2)        % don't calculate VIF for column of ones
            xy = zscore(X1_b(:,i));
            xx=X1_b(:, [1:i-1 i+1:size(X1_b, 2)]);
            xx(:, 2:end)=zscore(xx(:, 2:end));
            [~,~,R,~,STATS] = regress(xy,xx);
            VIF_v1b(i)=1./(1-STATS(1));
        end
        
        subj_VIF_v1b(i0,:) = VIF_v1b;
    end
    
    
    
    
%     for iReps = 1:size(probeModEstSim_noisySubjective,2)
%          % PREDICTION ERROR
%         mcPredErr = probeOutcomes(selID)' - probeModPredSim_noisySubjective(selID,iReps);
%         mcPredErr = mcPredErr - nanmean(mcPredErr);
%         
%     end
end

avBeta_v1_sim=nanmean(beta_v1_sim, 3);
devBeta_v1_sim=nanstd(beta_v1_sim, [], 3);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Subject and model errors & Betting vs. SAC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



subjPredErrSAC_hi = [];
subjPredErrSAC_low = [];
subjPredReliability_SAC_low = [];
subjPredReliability_SAC_hi = [];
subjEstErrSAC_hi = [];
subjEstErrSAC_low = [];

subjPredErrSAC_wrtMean_hi = [];
subjPredErrSAC_wrtMean_low = [];
subjEstErrSAC_wrtMean_hi = [];
subjEstErrSAC_wrtMean_low = [];

subj_normPredErrSAC = [];


subjMeanPredErrSAC_hi = [];
subjMeanPredErrSAC_low = [];
subjMeanEstErrSAC_hi = [];
subjMeanEstErrSAC_low = [];

%optimal model predictions
modelPredErrSAC_hi = [];
modelPredErrSAC_low = [];


%Stability and uncertainty terms
modelStab_pure_SAC_low = [];
modelStab_pure_SAC_hi = [];
modelStab_subjective_SAC_low = [];
modelStab_subjective_SAC_hi = [];

modelUnc_pure_SAC_low = [];
modelUnc_pure_SAC_hi = [];
modelUnc_subjective_SAC_low = [];
modelUnc_subjective_SAC_hi = [];


stabTerm_SAC_hi = [];
stabTerm_SAC_low = [];
uncTerm_SAC_hi = [];
uncTerm_SAC_low = [];

%betting stuff
subj_bet_freq = [];
subjBetSAC_hi = [];
subjBetSAC_low = [];
bet_prop_hi = [];
bet_prop_low = [];
norm_bet_prop_hi = [];
norm_bet_prop_low = [];


uncTerm_true_SAC_hi = [];
uncTerm_true_SAC_low = [];




subjPredErr = probeOutcomes - probePrediction;
subjPredErr_wrtMean = probeMean - probePrediction;
subjEstErr = probeOutcomes - probeEstimate;
subjEstErr_wrtMean = probeMean - probeEstimate;

modPredErr = bsxfun(@minus,modPred(allIsProbe),probeOutcomes');
%modEstErr = bsxfun(@minus,probeModEstSim_pure,probeOutcomes');

uniqueStds = unique(fixStd);
numTAC = 6;


% Removing "irreducible uncertainty from RU"
aux_ru_hi = (1- probeModelUncSubj).*(probeStd==20)';
aux_ru_hi(aux_ru_hi~=0) = aux_ru_hi(aux_ru_hi~=0) - nanmean(aux_ru_hi(aux_ru_hi~=0));
aux_ru_lo = (1- probeModelUncSubj).*(probeStd==10)';
aux_ru_lo(aux_ru_lo~=0) = aux_ru_lo(aux_ru_lo~=0) - nanmean(aux_ru_lo(aux_ru_lo~=0));

true_ru = aux_ru_hi + aux_ru_lo;



% --------------- some variables for Figure 2 ----------------------
bin_start = -20;
window_size = 20;
window_step = 5;
bin_end = 200 - window_size;
bin_edges = bin_start:window_step:bin_end;


%corrected estimates
probeCorrectedEstimate = probeEstimate;

%       ----------------- PREDICTION ----------------------

median_subj_reports_pred = [];
median_subj_means_pred = [];

for i0 = 1:length(unique(ID))
    for i1 = 1:length(bin_edges-1)
        sel_subj = probeSubjID == i0 & probeDay1Sess' == 0 & probeOutcomes >= bin_edges(i1) & ...
            probeOutcomes < bin_edges(i1) + window_size & ...
            probeTAC ~= 1;
        
        aux_pred = probePrediction(sel_subj);
        aux_outcome = probeMean(sel_subj);
        
        median_subj_reports_pred(i0,i1) = median(aux_pred);
        median_subj_means_pred(i0,i1) = median(aux_outcome);
    end
end


median_subj_reports_est = [];
median_subj_outcomes_est = [];

for i0 = 1:length(unique(ID))
    for i1 = 1:length(bin_edges-1)
        sel_subj = probeSubjID == i0 & probeDay1Sess' == 0 & probeOutcomes >= bin_edges(i1) & ...
            probeOutcomes < bin_edges(i1) + window_size;
        aux_estimate = probeEstimate(sel_subj);
        aux_outcome = probeOutcomes(sel_subj);
        
        median_subj_reports_est(i0,i1) = median(aux_estimate);
        median_subj_outcomes_est(i0,i1) = median(aux_outcome);
        
        
        probeCorrectedEstimate(sel_subj) = probeCorrectedEstimate(sel_subj) + ...
            (median(aux_outcome) - median(aux_estimate));
    end
end




std_subj_reports_training = [];
training_angles = unique(training_all_outcome);

for i0 = 1:length(unique(ID))
    
    sel_subj = training_all_ID == i0;
    
    aux_reports = training_all_estimate(sel_subj);
    aux_outcomes = training_all_outcome(sel_subj);
    std_subj_reports_training(i0) = std(aux_reports-aux_outcomes);
    
end




clear subj_est_err_disp subj_pred_err_disp


for i1 = 1:length(uniqueID)
    
    subj_bet_freq(i1) = nanmean(probeWager(probeSubjID == i1));
    
    
    
    for i2 = 1:numTAC
        
        if(use_day1)
            selSAC_hi = probebCPTAC == i2 & probeSubjID == i1 & probeStd == 20;
            selSAC_low = probebCPTAC == i2 & probeSubjID == i1 & probeStd == 10;
        else
            selSAC_hi = probebCPTAC == i2 & probeSubjID == i1 & probeStd == 20 & probeDay1Sess' == 0 ;
            selSAC_low = probebCPTAC == i2 & probeSubjID == i1 & probeStd == 10 & probeDay1Sess' == 0 ;
        end
        
        
        
        subjPredErrSAC_hi(i2,i1) = nanstd(subjPredErr(selSAC_hi));
        subjPredErrSAC_low(i2,i1) = nanstd(subjPredErr(selSAC_low));
        subjEstErrSAC_hi(i2,i1) = nanstd(subjEstErr(selSAC_hi));
        subjEstErrSAC_low(i2,i1) = nanstd(subjEstErr(selSAC_low));
        
        subjPredErrSAC_wrtMean_hi(i2,i1) = nanstd(subjPredErr_wrtMean(selSAC_hi));
        subjPredErrSAC_wrtMean_low(i2,i1) = nanstd(subjPredErr_wrtMean(selSAC_low));
        subjEstErrSAC_wrtMean_hi(i2,i1) = nanstd(subjEstErr_wrtMean(selSAC_hi));
        subjEstErrSAC_wrtMean_low(i2,i1) = nanstd(subjEstErr_wrtMean(selSAC_low));
        
        
        
        subjMeanPredErrSAC_hi(i2,i1) = nanmean(subjPredErr(selSAC_hi));
        subjMeanPredErrSAC_low(i2,i1) = nanmean(subjPredErr(selSAC_low));
        subjMeanEstErrSAC_hi(i2,i1) = nanmean(subjEstErr(selSAC_hi));
        subjMeanEstErrSAC_low(i2,i1) = nanmean(subjEstErr(selSAC_low));
        
        
        modelPredErrSAC_hi(i2,i1) = nanmean(nanstd(modPredErr(selSAC_hi)));
        modelPredErrSAC_low(i2,i1) = nanmean(nanstd(modPredErr(selSAC_low)));
        
        
        %         modelStab_subjective_SAC_hi(i2,i1) = 1- nanmean(mean(probeModStabSim_subjective(selSAC_hi,:)));
        %         modelStab_subjective_SAC_low(i2,i1) = 1- nanmean(mean(probeModStabSim_subjective(selSAC_low,:)));
        %
        %         modelStab_pure_SAC_hi(i2,i1) = 1- nanmean(mean(probeModStabSim_pure(selSAC_hi,:)));
        %         modelStab_pure_SAC_low(i2,i1) = 1- nanmean(mean(probeModStabSim_pure(selSAC_low,:)));
        %
        %
        %         modelUnc_subjective_SAC_hi(i2,i1) = nanmean(mean(probeModUncSim_subjective(selSAC_hi,:)));
        %         modelUnc_subjective_SAC_low(i2,i1) = nanmean(mean(probeModUncSim_subjective(selSAC_low,:)));
        %
        %         modelUnc_pure_SAC_hi(i2,i1) = nanmean(mean(probeModUncSim_pure(selSAC_hi,:)));
        %         modelUnc_pure_SAC_low(i2,i1) = nanmean(mean(probeModUncSim_pure(selSAC_low,:)));
        %
        
        %betting
        
        bet_prop_hi(i2,i1) = nanmean(probeWager(selSAC_hi));
        bet_prop_low(i2,i1) = nanmean(probeWager(selSAC_low)) ;
        
        norm_bet_prop_hi(i2,i1) = bet_prop_hi(i2,i1) - subj_bet_freq(i1);
        norm_bet_prop_low(i2,i1) = bet_prop_low(i2,i1) - subj_bet_freq(i1);
        
        % Uncertainty and stability terms as a function of SAC
        stabTerm_SAC_hi(i2,i1) = nanmean(probeModelStabSubj(selSAC_hi));
        stabTerm_SAC_low(i2,i1) = nanmean(probeModelStabSubj(selSAC_low));
        
        uncTerm_SAC_hi(i2,i1) = nanmean(probeModelUncSubj(selSAC_hi));
        uncTerm_SAC_low(i2,i1) = nanmean(probeModelUncSubj(selSAC_low));
        
        uncTerm_true_SAC_hi(i2,i1) = nanmean(true_ru(selSAC_hi));
        uncTerm_true_SAC_low(i2,i1) = nanmean(true_ru(selSAC_low));
        
    end
    
    subjPredReliability_SAC_hi(:,i1) = ....
        subjPredErrSAC_hi(:,i1).^2 ./ ...
        ( subjPredErrSAC_hi(:,i1).^2 + std_subj_reports_training(i1).^2);
    
    
    subjPredReliability_SAC_low(:,i1) = ....
        subjPredErrSAC_low(:,i1).^2 ./ ...
        ( subjPredErrSAC_low(:,i1).^2 + std_subj_reports_training(i1).^2);
end

% Don't consider extreme betting subjects
sel_bet_subj = subj_bet_freq < 0.1 | subj_bet_freq > 0.9;

norm_bet_prop_hi(:,sel_bet_subj) = [];
norm_bet_prop_low(:,sel_bet_subj) = [];

%kkMakeErrFig




% CALCULATE THE THEORETICAL BOUNDS (and contrasts) FOR ESTIMATION ERRORS
cp_cont_wts = [1, ones(1, 5)./-5];
sac_cont_wts = [-2 -1 0 1 2];

theory_subjective_est_err_hi = bsxfun(@plus,1./subjPredErrSAC_hi.^2,1./std_subj_reports_training.^2);
theory_subjective_est_err_hi = sqrt(1./theory_subjective_est_err_hi);

theory_subjective_est_err_low = bsxfun(@plus,1./subjPredErrSAC_low.^2,1./std_subj_reports_training.^2);
theory_subjective_est_err_low = sqrt(1./theory_subjective_est_err_low);

theory_optimal_est_err_hi = bsxfun(@plus,1./modelPredErrSAC_hi.^2,1./std_subj_reports_training.^2);
theory_optimal_est_err_hi = sqrt(1./theory_optimal_est_err_hi);

theory_optimal_est_err_low = bsxfun(@plus,1./modelPredErrSAC_low.^2,1./std_subj_reports_training.^2);
theory_optimal_est_err_low = sqrt(1./theory_optimal_est_err_low);


aux_theory_subjective_EstErr_combined = (theory_subjective_est_err_low' + ...
    theory_subjective_est_err_hi')*0.5;
theory_estErr_subjective_cp_contrast = mean(bsxfun(@times, ...
    aux_theory_subjective_EstErr_combined, cp_cont_wts),2);
theory_estErr_subjective_sac_contrast = mean(bsxfun(@times, ...
    aux_theory_subjective_EstErr_combined(:,2:end), sac_cont_wts),2);


theory_estErr_subjective_cp_contrast_hi = mean(bsxfun(@times, ...
    theory_subjective_est_err_hi', cp_cont_wts),2);
theory_estErr_subjective_sac_contrast_hi = mean(bsxfun(@times, ...
    theory_subjective_est_err_hi(2:end,:)', sac_cont_wts),2);

theory_estErr_subjective_cp_contrast_low = mean(bsxfun(@times, ...
    theory_subjective_est_err_low', cp_cont_wts),2);
theory_estErr_subjective_sac_contrast_low = mean(bsxfun(@times, ...
    theory_subjective_est_err_low(2:end,:)', sac_cont_wts),2);

theory_estErr_subjective_noise_contrast = nanmean(theory_subjective_est_err_hi(2:end,:)' - ...
    theory_subjective_est_err_low(2:end,:)',2);



theory_bet_windowSize = 16;
theory_bet_frequency_hi = [];
theory_bet_frequency_low = [];
for i1 = 1:length(uniqueIDs)
    for i2 = 1:numTAC
        aux_mass = normcdf(theory_bet_windowSize/2,0,theory_subjective_est_err_low(i2,i1)/2) - ...
            normcdf(-theory_bet_windowSize/2,0,theory_subjective_est_err_low(i2,i1)/2);
        %aux_mass = aux_mass/(1-aux_mass);
        
        theory_bet_frequency_low(i2,i1) = aux_mass;
        
        aux_mass = normcdf(theory_bet_windowSize/2,0,theory_subjective_est_err_hi(i2,i1)/2) - ...
            normcdf(-theory_bet_windowSize/2,0,theory_subjective_est_err_hi(i2,i1)/2);
        theory_bet_frequency_hi(i2,i1) = aux_mass;
    end
end


theory_bet_frequency_low = bsxfun(@minus,theory_bet_frequency_low,mean(theory_bet_frequency_low));
theory_bet_frequency_hi = bsxfun(@minus,theory_bet_frequency_hi,mean(theory_bet_frequency_hi));

aux_theory_bet_combined = (theory_bet_frequency_low' + ...
    theory_bet_frequency_hi')*0.5;
aux_theory_bet_combined = aux_theory_bet_combined(~sel_bet_subj,:);
theory_bet_cp_contrast = mean(bsxfun(@times, ...
    aux_theory_bet_combined, cp_cont_wts),2);
theory_bet_sac_contrast = mean(bsxfun(@times, ...
    aux_theory_bet_combined(:,2:end), sac_cont_wts),2);
theory_bet_noise_contrast = nanmean(theory_bet_frequency_hi(2:end,~sel_bet_subj)' - ...
    theory_bet_frequency_low(2:end,~sel_bet_subj)',2);

theory_bet_cp_contrast_low = mean(bsxfun(@times, ...
    theory_bet_frequency_low', cp_cont_wts),2);
theory_bet_sac_contrast_low = mean(bsxfun(@times, ...
    theory_bet_frequency_low(2:end,:)', sac_cont_wts),2);

theory_bet_cp_contrast_hi = mean(bsxfun(@times, ...
    theory_bet_frequency_hi', cp_cont_wts),2);
theory_bet_sac_contrast_hi = mean(bsxfun(@times, ...
    theory_bet_frequency_hi(2:end,:)', sac_cont_wts),2);


% PREDICTION ERROR CONTRASTS FOR OPTIMAL MODEL
avNcp_cont_wts = [0, ones(1, 5)]./5;
cp_cont_wts = [1, ones(1, 5)./-5];
sac_cont_wts = [-2 -1 0 1 2];


aux_pe_model_combined = (modelPredErrSAC_hi' + modelPredErrSAC_hi')*0.5;
model_predErr_cp_contrast = mean(bsxfun(@times,aux_pe_model_combined, cp_cont_wts),2);
model_predErr_sac_contrast = mean(bsxfun(@times,aux_pe_model_combined(:,2:end), sac_cont_wts),2);
model_predErr_noise_contrast = nanmean(modelPredErrSAC_hi(2:end,:)' - ...
    modelPredErrSAC_low(2:end,:)',2);

model_predErr_cp_contrast_hi = mean(bsxfun(@times,modelPredErrSAC_hi', cp_cont_wts),2);
model_predErr_sac_contrast_hi = mean(bsxfun(@times,modelPredErrSAC_hi(2:end,:)', sac_cont_wts),2);

model_predErr_cp_contrast_low = mean(bsxfun(@times,modelPredErrSAC_low', cp_cont_wts),2);
model_predErr_sac_contrast_low = mean(bsxfun(@times,modelPredErrSAC_low(2:end,:)', sac_cont_wts),2);

% PREDICTION ERROR CONTRASTS FOR SUBJECTS



aux_pe_combined = (subjPredErrSAC_hi' + subjPredErrSAC_low')*0.5;

predErr_avNcp_contrast = sum(bsxfun(@times,aux_pe_combined, avNcp_cont_wts),2);
predErr_cp_contrast = mean(bsxfun(@times,aux_pe_combined, cp_cont_wts),2);
predErr_sac_contrast = mean(bsxfun(@times,aux_pe_combined(:,2:end), sac_cont_wts),2);
predErr_noise_contrast = nanmean(subjPredErrSAC_hi(2:end,:)' - ...
    subjPredErrSAC_low(2:end,:)',2);

predErr_cp_contrast_hi = mean(bsxfun(@times,subjPredErrSAC_hi', cp_cont_wts),2);
predErr_sac_contrast_hi = mean(bsxfun(@times,subjPredErrSAC_hi(2:end,:)', sac_cont_wts),2);

predErr_cp_contrast_lo = mean(bsxfun(@times,subjPredErrSAC_low', cp_cont_wts),2);
predErr_sac_contrast_lo = mean(bsxfun(@times,subjPredErrSAC_low(2:end,:)', sac_cont_wts),2);





%ESTIMATION ERROR CONTRASTS

cp_cont_wts = [1, ones(1, 5)./-5];
sac_cont_wts = [-2 -1 0 1 2];

aux_estE_combined = (subjEstErrSAC_hi' + subjEstErrSAC_low')*0.5;

estErr_cp_contrast = mean(bsxfun(@times,aux_estE_combined, cp_cont_wts),2);
estErr_sac_contrast = mean(bsxfun(@times,aux_estE_combined(:,2:end), sac_cont_wts),2);
estErr_noise_contrast = nanmean(subjEstErrSAC_hi(2:end,:)' - ...
    subjEstErrSAC_low(2:end,:)',2);

estErr_cp_contrast_hi = mean(bsxfun(@times,subjEstErrSAC_hi', cp_cont_wts),2);
estErr_sac_contrast_hi = mean(bsxfun(@times,subjEstErrSAC_hi(2:end,:)', sac_cont_wts),2);

estErr_cp_contrast_lo = mean(bsxfun(@times,subjEstErrSAC_low', cp_cont_wts),2);
estErr_sac_contrast_lo = mean(bsxfun(@times,subjEstErrSAC_low(2:end,:)', sac_cont_wts),2);


% BETTING CONTRASTS


cp_cont_wts = [1, ones(1, 5)./-5];
sac_cont_wts = [-2 -1 0 1 2];

aux_bet_combined = (norm_bet_prop_hi' + norm_bet_prop_low')*0.5;

bet_cp_contrast = mean(bsxfun(@times,aux_bet_combined, cp_cont_wts),2);
bet_sac_contrast = mean(bsxfun(@times,aux_bet_combined(:,2:end), sac_cont_wts),2);
bet_noise_contrast = nanmean(norm_bet_prop_hi(2:end,:)' - ...
    norm_bet_prop_low(2:end,:)',2);

bet_cp_contrast_hi = mean(bsxfun(@times,norm_bet_prop_hi', cp_cont_wts),2);
bet_sac_contrast_hi = mean(bsxfun(@times,norm_bet_prop_hi(2:end,:)', sac_cont_wts),2);

bet_cp_contrast_lo = mean(bsxfun(@times,norm_bet_prop_low', cp_cont_wts),2);
bet_sac_contrast_lo = mean(bsxfun(@times,norm_bet_prop_low(2:end,:)', sac_cont_wts),2);







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Compare model and subject prior weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note (Kamesh): probeModEstSim is the variable that has the model estimation
%error


uniqueStds = [10,20];
uniqueIDs = unique(ID);
%uniqueID([5]) = [];
nReps = size(probeModEstSim,2);

% calculate prior weights vs TAC for the model
clear  beta_model_pw beta_subj_pw

for i0 = 1:length(uniqueIDs)
    
    if(use_day1)
        selID = probeSubjID == uniqueIDs(i0);
    else
        selID = probeSubjID == uniqueIDs(i0) & probeDay1Sess' == 0 &  probebCPTAC < 7;
    end
    
    clear Bsubj stats
    subjPredErr = probeOutcomes' - probePrediction';
    subjPercErr = probeOutcomes' - probeEstimate';
    
    %DEBUG
    %selID = selID & abs(subjPercErr') > -1 & abs(subjPercErr') < 90;
    
    selSTD = probeStd(selID)';
    
    perceptualErr = subjPercErr(selID)-nanmean(subjPercErr(selID));
    predictionErr = subjPredErr(selID) -nanmean(subjPredErr(selID));
    
    %create the decision matrix
    X0 = ones(length(predictionErr),1);
    
    % columns corresponding to b1 - b6
    Xpred = zeros(length(predictionErr),2*numTAC);
    for i4 = 1:numTAC
        Xpred(probebCPTAC(selID) == i4, i4) = 1;
        Xpred(probebCPTAC(selID) == i4, i4+6) = 1;
    end
    Xpred(:,1:6) = bsxfun(@times,Xpred(:,1:6),predictionErr.*(selSTD==10));
    Xpred(:,7:12) = bsxfun(@times,Xpred(:,7:12),predictionErr.*(selSTD==20));
    Xpred = bsxfun(@minus,Xpred,nanmean(Xpred,1));
    
    % find indices of big change-point trials
    bCPIdx = probebCPTAC < 7;
    
    
    rawPred = probePrediction(selID)';
    
    % how far left of center are we?
    XSpace = (90 - rawPred);
    XSpace = (XSpace - nanmean(XSpace))./std(XSpace);
    
    
    Xsubj = [ones(size(Xpred,1),1) Xpred XSpace];
    
    %[Bsubj, stats] = robustfit(Xsubj,perceptualErr,[],[], 'off');
    
    [Bsubj, ~,~,LLsubj] = regressW_mike(perceptualErr, ...
        trialExpStd(selID),Xsubj);
    
    
    beta_subj_pw(:,i0) = Bsubj;
    
    % Now do the regression for the model
    clear Bmodel
    for iReps = 1:nReps
        
        
        predErr = probeOutcomes'- probeModPredSim(:,iReps);
        percErr = probeOutcomes'- probeModEstSim(:,iReps);
        modTrialExpStd = trialExpStd;
        perceptualErr = percErr(selID)-nanmean(percErr(selID));
        predictionErr = predErr(selID) -nanmean(predErr(selID));
        
        %create the decision matrix
        X0 = ones(length(predictionErr),1);
        
        % columns corresponding to b1 - b6
        Xpred = zeros(length(predictionErr),2*numTAC);
        for i4 = 1:numTAC
            Xpred(probebCPTAC(selID) == i4, i4) = 1;
            Xpred(probebCPTAC(selID) == i4, i4+6) = 1;
        end
        Xpred(:,1:6) = bsxfun(@times,Xpred(:,1:6),predictionErr.*(selSTD==10));
        Xpred(:,7:12) = bsxfun(@times,Xpred(:,7:12),predictionErr.*(selSTD==20));
        Xpred = bsxfun(@minus,Xpred,nanmean(Xpred,1));
        
        rawPred = modPred(selID);
        
        % how far left of center are we?
        XSpace = (90 - rawPred);
        XSpace = (XSpace - nanmean(XSpace))./std(XSpace10);
        
        Xmodel = [ones(size(Xpred,1),1) Xpred XSpace];
        %[Bmodel, stats] = robustfit(Xmodel,perceptualErr,[],[], 'off');
        
        [Bmodel, ~,~,LLmodel] = regressW_mike(perceptualErr, ...
            trialExpStd(selID),Xmodel);
        
        beta_model_pw(:,i0,iReps) = Bmodel;
        
    end
    
    
end



avg_beta_model_pw = nanmean(beta_model_pw,3);


clear modBetaInt_hi modBetaInt_lo modAvgBetas_hi modAvgBetas_lo

modBetaInt_hi = std(avg_beta_model_pw(8:13,:),[],2)./sqrt(i0);
modAvgBetas_hi = mean(avg_beta_model_pw(8:13,:),2);

modBetaInt_lo = std(avg_beta_model_pw(2:7,:),[],2)./sqrt(i0);
modAvgBetas_lo = mean(avg_beta_model_pw(2:7,:),2);

avg_beta_model_lo = avg_beta_model_pw(2:7,:);
avg_beta_model_hi = avg_beta_model_pw(8:13,:);

subjAvgBetas_lo =  mean(beta_subj_pw(2:7,:),2);
subjBetaInt_lo = std(beta_subj_pw(2:7,:),[],2)./sqrt(i0);

subjAvgBetas_hi = mean(beta_subj_pw(8:13,:),2);
subjBetaInt_hi = std(beta_subj_pw(8:13,:),[],2)./sqrt(i0);

beta_subj_lo = beta_subj_pw(2:7,:);
beta_subj_hi = beta_subj_pw(8:13,:);


% ------------ Get the theoretical prior bias using subjective and optimal
% priors

clear theoretical_subjective_priorBias_hi theoretical_subjective_priorBias_low ...
    theoretical_optimal_priorBias_hi theoretical_optimal_priorBias_low

aux_priorBias_hi = bsxfun(@plus,subjPredErrSAC_hi.^2,std_subj_reports_training.^2);
theoretical_subjective_priorBias_hi = bsxfun(@times,1./aux_priorBias_hi,...
    std_subj_reports_training.^2);

aux_priorBias_low = bsxfun(@plus,subjPredErrSAC_low.^2,std_subj_reports_training.^2);
theoretical_subjective_priorBias_low = bsxfun(@times,1./aux_priorBias_low,...
    std_subj_reports_training.^2);

theoryBias_subjective_Int_hi = std(theoretical_subjective_priorBias_hi,[],2)./sqrt(29);
theoryBias_subjective_avg_hi = mean(theoretical_subjective_priorBias_hi,2);

theoryBias_subjective_Int_low = std(theoretical_subjective_priorBias_low,[],2)./sqrt(29);
theoryBias_subjective_avg_low = mean(theoretical_subjective_priorBias_low,2);


aux_priorBias_hi = bsxfun(@plus,modelPredErrSAC_hi.^2,std_subj_reports_training.^2);
theoretical_optimal_priorBias_hi = bsxfun(@times,1./aux_priorBias_hi,...
    std_subj_reports_training.^2);

aux_priorBias_low = bsxfun(@plus,modelPredErrSAC_low.^2,std_subj_reports_training.^2);
theoretical_optimal_priorBias_low = bsxfun(@times,1./aux_priorBias_low,...
    std_subj_reports_training.^2);

theoryBias_optimal_Int_hi = std(theoretical_optimal_priorBias_hi,[],2)./sqrt(29);
theoryBias_optimal_avg_hi = mean(theoretical_optimal_priorBias_hi,2);

theoryBias_optimal_Int_low = std(theoretical_optimal_priorBias_low,[],2)./sqrt(29);
theoryBias_optimal_avg_low = mean(theoretical_optimal_priorBias_low,2);

% Compute the contrasts for use in figures


cp_cont_wts = [1, ones(1, 5)./-5];
sac_cont_wts = [-2 -1 0 1 2];
jit=.05;

% auxBetas_hi = beta_subj_hi(pltIdx(2:end),:)';
% auxBetas_lo = beta_subj_lo(pltIdx(2:end),:)';

aux_beta_combined =  (beta_subj_hi'+beta_subj_lo')*0.5;

pw_sac_contrast_combined = mean(bsxfun(@times,aux_beta_combined(:,2:end),sac_cont_wts),2);
pw_sac_contrast_hi = mean(bsxfun(@times,beta_subj_hi(2:end,:)',sac_cont_wts),2);
pw_sac_contrast_lo = mean(bsxfun(@times,beta_subj_lo(2:end,:)',sac_cont_wts),2);

pw_avNcp_contrast=sum(bsxfun(@times, aux_beta_combined, avNcp_cont_wts),2);

pw_cp_contrast_combined = mean(bsxfun(@times,aux_beta_combined,cp_cont_wts),2);
pw_cp_contrast_hi = mean(bsxfun(@times,beta_subj_hi',cp_cont_wts),2);
pw_cp_contrast_lo = mean(bsxfun(@times,beta_subj_lo',cp_cont_wts),2);

pw_noise_contrast = nanmean(beta_subj_hi' - beta_subj_lo',2);



% Now, calculate the contrasts for the theoretical prior bias
aux_optimal_priorBias_combined =  (theoretical_optimal_priorBias_low'+ ...
    theoretical_optimal_priorBias_hi')*0.5;
theory_optimal_pw_cp_contrast = mean(bsxfun(@times,aux_optimal_priorBias_combined,cp_cont_wts),2);
theory_optimal_pw_sac_contrast = mean(bsxfun(@times,aux_optimal_priorBias_combined(:,2:end),sac_cont_wts),2);
theory_optimal_pw_noise_contrast = nanmean(theoretical_optimal_priorBias_hi' - ...
    theoretical_optimal_priorBias_low',2);


aux_subjective_priorBias_combined =  (theoretical_subjective_priorBias_low'+ ...
    theoretical_subjective_priorBias_hi')*0.5;
theory_subjective_pw_cp_contrast = mean(bsxfun(@times,aux_subjective_priorBias_combined,cp_cont_wts),2);
theory_subjective_pw_sac_contrast = mean(bsxfun(@times,aux_subjective_priorBias_combined(:,2:end),sac_cont_wts),2);
theory_subjective_pw_noise_contrast = nanmean(theoretical_subjective_priorBias_hi' - ...
    theoretical_subjective_priorBias_low',2);

theory_subjective_pw_cp_contrast_hi = mean(bsxfun(@times,theoretical_subjective_priorBias_hi',cp_cont_wts),2);
theory_subjective_pw_sac_contrast_hi = mean(bsxfun(@times,theoretical_subjective_priorBias_hi(2:end,:)',sac_cont_wts),2);

theory_subjective_pw_cp_contrast_low = mean(bsxfun(@times,theoretical_subjective_priorBias_low',cp_cont_wts),2);
theory_subjective_pw_sac_contrast_low = mean(bsxfun(@times,theoretical_subjective_priorBias_low(2:end,:)',sac_cont_wts),2);



%% ---------------- Posterior Predictive Checking w/ Selected Model -----------
%                  (Comparing subject and model prior wt. vs. SAC)
%                        (NOISE CONDITIONS SPLIT)


% Calculate prior weight for the noise conditions collapsed
uniqueStds = [10,20];
uniqueIDs = unique(ID);

% calculate prior weights vs TAC for the model
clear  linModel_beta_pw linModel_sim_pw_split subj_pw_allTrials

for i0 = 1:length(uniqueIDs)
    
    if(use_day1)
        selID = probeSubjID == uniqueIDs(i0);
    else
        selID = probeSubjID == uniqueIDs(i0) & probeDay1Sess' == 0;
    end
    selSTD = probeStd(selID)';
    
    clear Bsubj stats BlinMod
    subjPredErr = probeOutcomes' - probePrediction';
    subjPercErr = probeOutcomes' - probeEstimate';
    
    perceptualErr = subjPercErr(selID)-nanmean(subjPercErr(selID));
    predictionErr = subjPredErr(selID) -nanmean(subjPredErr(selID));
    
    
    %create the decision matrix
    X0 = ones(length(predictionErr),1);

    % columns corresponding to b1 - b6
    Xpred = zeros(length(predictionErr),2*numTAC);
    for i4 = 1:numTAC
        Xpred(probeTAC(selID) == i4, i4) = 1;
        Xpred(probeTAC(selID) == i4, i4+6) = 1;
    end
    Xpred(:,1:6) = bsxfun(@times,Xpred(:,1:6),predictionErr.*(selSTD==10));
    Xpred(:,7:12) = bsxfun(@times,Xpred(:,7:12),predictionErr.*(selSTD==20));
    Xpred = bsxfun(@minus,Xpred,nanmean(Xpred,1));
    
    rawPred = probePrediction(selID)';
    
    % how far left of center are we?
    XSpace = (90 - rawPred);
    XSpace = (XSpace - nanmean(XSpace))./std(XSpace);
    
    Xsubj = [ones(size(Xpred,1),1) Xpred XSpace];
    
    %[Bsubj, stats] = robustfit(Xsubj,perceptualErr,[],[], 'off');
    
    [Bsubj, ~,~,LLsubj] = regressW_mike(perceptualErr, ...
        trialExpStd(selID),Xsubj);
    
    
    %do the regression w/ perceptual error simulated from linear model
    linMod_perceputal_err = linMod_v1_sim_perceptual_err(selID)';
    [BlinMod, ~,~,LLlinMod] = regressW_mike(linMod_perceputal_err, ...
        trialExpStd(selID),Xsubj);

    
    subj_pw_allTrials(:,i0) = Bsubj;
    linModel_sim_pw_split(:,i0) =  BlinMod;
    
    % run through perceptions simulated from the generative model and
    % compute the same thing:
    clear BgenMod
    for j = 1:size(probeModEstSim,2)
        genMod_perceputal_err= probeOutcomes(selID)' - probeModEstSim(selID,j);
        BgenMod(j,:) = regressW_mike(genMod_perceputal_err, trialExpStd(selID),Xsubj);
        
    end
        
    genModel_sim_pw_split(:,i0) =  nanmean(BgenMod);
    
end

subjAvgBetas_allTrials_lo =  mean(subj_pw_allTrials(2:7,:),2);
subjBetaInt_allTrials_lo = std(subj_pw_allTrials(2:7,:),[],2)./sqrt(i0);

subjAvgBetas_allTrials_hi = mean(subj_pw_allTrials(8:13,:),2);
subjBetaInt_allTrials_hi = std(subj_pw_allTrials(8:13,:),[],2)./sqrt(i0);


linearApprxAvgBetas_split_lo = mean(linModel_sim_pw_split(2:7,:),2);
linearApprxBetaInt_split_lo =  std(linModel_sim_pw_split(2:7,:),[],2)./sqrt(i0);

linearApprxAvgBetas_split_hi = mean(linModel_sim_pw_split(8:13,:),2);
linearApprxBetaInt_split_hi =  std(linModel_sim_pw_split(8:13,:),[],2)./sqrt(i0);


% MRN added this to measure perceptual bias in the generative model:

genModAvgBetas_split_lo = mean(genModel_sim_pw_split(2:7,:),2);
genModBetaInt_split_lo =  std(genModel_sim_pw_split(2:7,:),[],2)./sqrt(i0);

genModAvgBetas_split_hi = mean(genModel_sim_pw_split(8:13,:),2);
genModBetaInt_split_hi =  std(genModel_sim_pw_split(8:13,:),[],2)./sqrt(i0);











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            PUPIL ANALYSIS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% step 1: do analysis
% step 2: make figures
% step 3: cleanup code and ensure that all necessary variables are in main
% script. 


% things we need -- excluding first session:
% 1) selection for subjects and day1 data:
% 2) prediction, outcome, perception
% 3) pupil trace for fixated trials
  

  %% pupil analysis -- MRN
  % but i'm not sure how kamesh is dealing with low frequency stuff... 
  
  trialsToDrop  = 0;  % we may want to drop a few trials at the beginning of each block?
  numComps      = 3; % uncertainty effect requires low frequency effects to be absorbed...
  

% create spatialbias term for mixed effects regression model:  
  for i = 1:max(probeSubjID)
       sel=probeSubjID==i;
      [~, aux_spatial_idx] = min(abs(bsxfun(@minus,probeOutcomes(sel)',mean_train_angles')),[],2);
      mcSpatialBias(sel) = train_PE_subject_bias_binned(i,aux_spatial_idx)'; 
  end
  
  
  % Create variables selected for fixated trials that were not from session
  %             #1
  fixIsFixed = logical(fixIsFixed)';
  fixPriorWtInfResid = priorWtInfResid_v1(fixIsFixed);


  nofirst_PercErr        = fixOutcomes(~fixDay1Sess)'-fixEstimate(~fixDay1Sess)';
  nofirst_PredErr        = fixOutcomes(~fixDay1Sess)'-fixPrediction(~fixDay1Sess)';
  nofirst_Stab           = meanCentX(probeModelStabSubj(logical(probeIsFixed&~probeDay1Sess')));
  nofirst_Rel            = meanCentX(probeModelUncSubj(logical(probeIsFixed&~probeDay1Sess')));
  nofirst_Bet            = meanCentX(fixWager(~fixDay1Sess)');
  nofirst_Std_invCode    = fixStd(~fixDay1Sess)'.*-1 + 15;
  nofirst_Space          = fixPrediction(~fixDay1Sess)';
  nofirst_TimeToPred     = fixTimeToPredict(~fixDay1Sess)';
  % fixTsinceLastProb
  nofirst_CentErr        = 90 - fixPrediction(~fixDay1Sess)';
  nofirst_BlkNum         =  probeBlkNum(logical(probeIsFixed)&~probeDay1Sess');
  nofirst_SesNum         = probeSessionID(logical(probeIsFixed)&~probeDay1Sess');
  nofirst_TrialsWB       = fixTrialNo(~fixDay1Sess)';
  nofirst_tSinceLastProbe=fixTsinceLastProbe(~fixDay1Sess)';
  nofirst_trialNumber    =fixTrialNo(~fixDay1Sess)';
  nofirst_sesID          =fixSessionID(~fixDay1Sess)';
  nofirst_expStd         =trialExpStd(probeIsFixed&~probeDay1Sess');
  noFirst_spaceBias      = mcSpatialBias(probeIsFixed&~probeDay1Sess');
  
  noFirst_residPriorWt   =fixPriorWtInfResid(~fixDay1Sess);
  noFirst_subjID         =fixSubjID(~fixDay1Sess)';
  nofirst_TrialNo        =fixTrialNo(~fixDay1Sess)';


  for i =1:length(probeSessionID)
      if i==1||probeSessionID(i)~=probeSessionID(i-1)
          probeTrialWithinSes(i)=1;
      else
          probeTrialWithinSes(i)=probeTrialWithinSes(i-1)+1;
      end
  end  % trial within session as a potential confounding factor

 % create cosine series of regressors to soak of really low frequency
 % changes
  allCosComps=nan(length(probeSessionID), numComps);
  clear allCosComps
  for i =1:max(fixSessionID);
      sesTrials=fixSessionID==i;     
      ns=1:sum(sesTrials);
      clear cosComp
      for k =1:numComps;
          for n = ns
              cosComp(n,k)= cos((pi.*(2.*n - 1).*(k-1))./(2.*sum(sesTrials)));
          end
      end
      allCosComps(sesTrials,:)=cosComp;
  end
  
  nofirst_rawPupil=fixRawPupil(~fixDay1Sess,:);
  baseLine=nanmean(nofirst_rawPupil(:,1:3), 2);
  pTrialBaseline=[nan; baseLine(1:end-1)];
  nTrialBaseline=[baseLine(2:end); nan];

  nofirst_allCosComps=allCosComps(~fixDay1Sess,:);
  
  clear simp_baseP simp_baseF simp_baseR2 baseP baseF baseR2 baseCoeffs respCoeffs respResidual respCoeffs avg_subj_pup_traceMRN
  
  regLabels={'Relevance', 'Reliability', 'Bet', 'Residual', 'trialSinceLastProbe', 'prevTrialDiameter'}
  
  % NOTE: 6th regressor should always be previous trial diameter, this is
  % hard coded to be replaced as we loop through timepoints (with previous
  % trial data from the relevant timepoint). 
  
  
 avg_subj_pup_trace_KK = []; %evoked response
  
  

 load('respTime.mat');

  
  
  
  
  for i = 1:max(fixSubjID)
      
      % create cosine components:
      % question: should we include noise or just soak up low freq with cos componenets? , : fixStd(sel)',
      sel=noFirst_subjID==i &nofirst_TrialsWB>trialsToDrop;
      
      % store average pupil trace for individual differences analysis:
      avg_subj_pup_traceMRN(i,:)= nanmean(nofirst_rawPupil(sel,:));
      
      % looking at the evoked response (relative to baseline)
      avg_subj_pup_trace_KK(i,:) = nanmean(bsxfun(@minus,nofirst_rawPupil(sel,:),baseLine(sel)),1);
      
      xMat=zScoreX([nofirst_Stab(sel), nofirst_Rel(sel), nofirst_Bet(sel),  noFirst_residPriorWt(sel), nofirst_tSinceLastProbe(sel), pTrialBaseline(sel)]);
      numRawCoeffs=size(xMat, 2);
      selSes=unique(nofirst_sesID(sel));
      for j = 1:length(selSes)
          xMat=[xMat, nofirst_allCosComps(sel,:).* repmat(fixSessionID(sel)==selSes(j), 1, numComps)];
      end
      
      % compute stats on a null model for comparison:
      null_base_logLike(i)= sum(log(normpdf(zscore(baseLine(sel)), 0, 1)));  
      [null_base_BIC(i), null_base_AIC(i)]=computeBIC(null_base_logLike(i), 1, sum(sel));
      
      
      % run reduced baseline regression -- to get goodness of fit info. 
      [~,~,R, ~, STATS] =regress(zscore(baseLine(sel)), [ones(sum(sel),1),  xMat(:,1:4) ]);
      
      % store goodness of fit stuff:
      simp_base_logLike(i)=sum(log(normpdf(R, 0, nanstd(R))));
      [simp_base_BIC(i), simp_base_AIC(i)]=computeBIC(simp_base_logLike(i), 5, sum(sel));
      simp_baseR2(i)=STATS(1);
      simp_baseF(i)=STATS(2);
      simp_baseP(i)=STATS(3);

      % run baseline regression:
      [b,~,R, RINT, STATS] =regress(zscore(baseLine(sel)), xMat(:,:));
      baseCoeffs(i,:)=b(1:numRawCoeffs);
      
      % Store goodness of fit info for best model:      
      full_base_logLike(i)=sum(log(normpdf(R, 0, nanstd(R))));
      [full_base_BIC(i), full_base_AIC(i)]=computeBIC(full_base_logLike(i), size(xMat,2), sum(sel));
      baseR2(i)=STATS(1);
      baseF(i)=STATS(2);
      baseP(i)=STATS(3);     
      [B,BINT,baseResidual(sel)] = regress(baseLine(sel),xMat(:,5:end));
      
      
      
      % run other regressions:
      % Add baseline from CURRENT trial into the regression model:
      xMat(:,end+1)=baseLine(sel);
      % from an older analysis. 
      for j = 1:size(nofirst_rawPupil, 2)
          % replace previous baseline with previous diamater
          Y=zscore(nofirst_rawPupil(sel,j));
          xMat(:,6)= [nan; Y(1:end-1)]; % add pupil measurment from PREVIOUS trial into regrssion model. 
          b=regress(nofirst_rawPupil(sel,j), xMat(:,:));
          respCoeffs(i,j,:)=b(1:numRawCoeffs);
          
          if j ==respTime

              % compute stats on a null model for comparison:
              null_resp_logLike(i)= nansum(log(normpdf(zscore(nofirst_rawPupil(sel,j)), 0, 1)));
              [null_resp_BIC(i), null_resp_AIC(i)]=computeBIC(null_resp_logLike(i), 1, sum(sel));
              
              % Simple model (no nuisance variables)
              [~,~,R,~,STATS]=regress(zscore(nofirst_rawPupil(sel,j)), [ones(length(xMat),1), xMat(:,1:4)]);
               simp_resp_logLike(i)= nansum(log(normpdf(R, 0, nanstd(R))));
              [simp_resp_BIC(i), simp_resp_AIC(i)]=computeBIC(simp_resp_logLike(i), 1, sum(sel));
              simp_respR2(i)=STATS(1);
              simp_respF(i)=STATS(2);
              simp_respP(i)=STATS(3);
              
              % Full model (with nuisance variables)
              [~,~,R,~,STATS]=regress(zscore(nofirst_rawPupil(sel,j)), xMat(:,:));
              full_resp_logLike(i)= nansum(log(normpdf(R, 0, nanstd(R))));
              [full_resp_BIC(i), full_resp_AIC(i)]=computeBIC(full_resp_logLike(i), 1, sum(sel));
              respR2(i)=STATS(1);
              respF(i)=STATS(2);
              respP(i)=STATS(3);
               
          end
          
          
          
          [B,BINT, respResidual(sel,j)] = regress(Y,xMat(:,5:end));
      end
  end
  
  xTickLabs={'Base model', 'w Nuisance'}
  rel_baseline_AIC= [simp_base_AIC' - null_base_AIC',  full_base_AIC'- null_base_AIC'];
  rel_response_AIC= [simp_resp_AIC' - null_resp_AIC',  full_resp_AIC'- null_resp_AIC'];
  
  
%   hold on
%   lims=[min(rel_baseline_AIC(:)) ,   max(rel_baseline_AIC(:))];
%   plot(lims, lims, '--k')
%   plot(rel_baseline_AIC(:,1), rel_baseline_AIC(:,2), 'o')
%   ylabel('Relative AIC (full model)')
%   xlabel('Relative AIC (base model)')
%   title('Baseline regression goodness of fit')
%   xlim(lims)
%   ylim(lims)
%   
%   
%   rel_response_AIC= [simp_resp_AIC' - null_resp_AIC',  full_resp_AIC'- null_resp_AIC'];
%   hold on
%   lims=[min(rel_response_AIC(:)) ,   max(rel_response_AIC(:))];
%   plot(lims, lims, '--k')
%   plot(rel_response_AIC(:,1), rel_response_AIC(:,2), 'o')
%   ylabel('Relative AIC (full model)')
%   xlabel('Relative AIC (base model)')
%   title('Baseline regression goodness of fit')
%   xlim(lims)
%   ylim(lims)

  
    
  %make_figureS2;

  
  % Baseline stats are easy... just one ttest:
  [H_base,P_base,CI,STATS]=ttest(baseCoeffs)
  
  
  % display stats for baseline coefficients
  % Remember: bet statistics are computed by ommitting subjects that always
  % give same PDW response (>.9|<.1)
  X1Labels=regLabels;
  meanBeta=nanmean(baseCoeffs);
  steBeta =nanstd(baseCoeffs)./sqrt(length(baseCoeffs));
  for i = 1:length(X1Labels)
      disp(sprintf('pupil baseline coefficient: %s \n \n \t Mean/STE beta= %s / %s \n \t t-stat = %s \n \t p-value= %s', ...
          X1Labels{i}, num2str(meanBeta(i)), num2str(steBeta(i)), num2str(STATS.tstat(i)), num2str(P_base(i))))
  end

  [h, p]=ttest(respCoeffs);
  
  
  
  
  
  
  % MRN 1/12/17 --
  % OK, we also want some basic descriptive plots showing the average pupil
  % diameter at baseline and time of response 
  
  % Goal: Create 6 panels consisting of pupil baseline and peak response
  % binned and sorted by "reliability", "relevance" and "confidence". 
    
  % subtract out mean subject effects:

  avPupTrace=nanmean(avg_subj_pup_traceMRN);
  respTime = find(avPupTrace==max(avPupTrace), 1)
  
  %save ~/Dropbox/auditoryTaskManuscript/respTime.mat respTime
  
 
  
  
  % get real values, not mean centered:  
  Stab           = (probeModelStabSubj(logical(probeIsFixed&~probeDay1Sess')));
  Rel            = (probeModelUncSubj(logical(probeIsFixed&~probeDay1Sess')));
  Bet            = (fixWager(~fixDay1Sess)');

  nBins=5
  pTiles=0: (100./nBins):100
  
  % Next steps:
  % 1) understand why confidence is not showing up
  % 2) look at residual diameter after accounting for other simple things
  % 3) individual differences in diameter at each timepoint.  
  
  % loop through subjects
  
  MC_pup=nan(size(nofirst_rawPupil));
  MC_respResidual=nan(length(nofirst_rawPupil),1);
  clear betBinResponse stabBinBaseline stabBinResponse relBinBaseline relBinResponse medStab medRel
  
  for i = 1:max(fixSubjID)
      sel=noFirst_subjID==i &nofirst_TrialsWB>trialsToDrop;      
      stabPrctiles= prctile(Stab(sel),pTiles);
      relPrctiles=  prctile(Rel(sel),pTiles);

      % loop through bins of stability and reliability:
      MC_pup(sel,:)=meanCentX(nofirst_rawPupil(sel,:));
      [B,BINT,R]=regress(MC_pup(sel,respTime),  [ones(sum(sel),1),  nanmean(MC_pup(sel,1:3),2)]);
      MC_respResidual(sel)=R;

      for j = 1:(length(stabPrctiles)-1)
          stabBinBaseline(i,j)=nanmean(nanmean(MC_pup(sel&Stab>=stabPrctiles(j)&Stab<stabPrctiles(j+1) ,1:3)));
          stabBinResponse(i,j)=(nanmean(MC_respResidual(sel&Stab>=stabPrctiles(j)&Stab<stabPrctiles(j+1) )));
          medStab(i,j)=nanmean(Stab(sel&Stab>=stabPrctiles(j)&Stab<stabPrctiles(j+1)));
          relBinBaseline(i,j)=nanmean(nanmean(MC_pup(sel&Rel>=relPrctiles(j)&Rel<relPrctiles(j+1) ,1:3)));
          relBinResponse(i,j)=(nanmean(MC_respResidual(sel&Rel>=relPrctiles(j)&Rel<relPrctiles(j+1))));
          medRel(i,j)=nanmean(Rel(sel&Rel>=relPrctiles(j)&Rel<relPrctiles(j+1)));
      end

      subBetFrac(i)=nanmean(Bet(sel));
      betBinBaseline(i,1)=nanmean(nanmean(MC_pup(sel&Bet==0 ,1:3)));
      betBinBaseline(i,2)=nanmean(nanmean(MC_pup(sel&Bet==1 ,1:3)));
      
      betBinResponse(i,1)=(nanmean(MC_respResidual(sel&Bet==0)));
      betBinResponse(i,2)=(nanmean(MC_respResidual(sel&Bet==1)));
      
      
      subMeanBaseline(i)=nanmean(nanmean(nofirst_rawPupil(sel,1:3)));
      subMeanResponse(i)=nanmean(nanmean(nofirst_rawPupil(sel,respTime)));

  end
    
  ste_stabBaseline=nanstd(stabBinBaseline)./sqrt(size(stabBinBaseline, 1));
  ste_stabResponse=nanstd(stabBinResponse)./sqrt(size(stabBinResponse, 1));
 
  ste_relBaseline=nanstd(relBinBaseline)./sqrt(size(relBinBaseline, 1));
  ste_relResponse=nanstd(relBinResponse)./sqrt(size(relBinResponse, 1));

  ste_betBaseline=nanstd(betBinBaseline(subjNormalBets))./sqrt(sum(subjNormalBets));
  ste_betResponse=nanstd(betBinResponse(subjNormalBets))./sqrt(sum(subjNormalBets));
 
  
  
  % Create figure:
%   subplot(2, 3, 1)
%   hold on
%   plot(nanmean(medStab), nanmean(stabBinBaseline), 'o')
%   plot([nanmean(medStab); nanmean(medStab)], ...
%       [nanmean(stabBinBaseline)+ste_stabBaseline; nanmean(stabBinBaseline)-ste_stabBaseline],  '-b')
%   
%   ylabel('Baseline pupil')
%   subplot(2, 3, 4)
%   hold on
%   plot(nanmean(medStab), nanmean(stabBinResponse), 'o')
%   plot([nanmean(medStab); nanmean(medStab)], ...
%       [nanmean(stabBinResponse)+ste_stabResponse; nanmean(stabBinResponse)-ste_stabResponse],  '-b')
%   ylabel('Pupil response')
%   xlabel('Stability')
%   
%   subplot(2, 3, 2)
%   hold on
%   plot(nanmean(medRel), nanmean(relBinBaseline), 'o')
%   plot([nanmean(medRel); nanmean(medRel)], ...
%       [nanmean(relBinBaseline)+ste_relBaseline; nanmean(relBinBaseline)-ste_relBaseline],  '-b')
%   
%   subplot(2, 3, 5)
%   hold on
%   plot(nanmean(medRel), nanmean(relBinResponse), 'o')
%   plot([nanmean(medRel); nanmean(medRel)], ...
%       [nanmean(relBinResponse)+ste_relResponse; nanmean(relBinResponse)-ste_relResponse],  '-b')
%   xlabel('Reliability')
% 
%   subplot(2, 3, 3)
%   hold on
%   plot([0, 1], nanmean(betBinBaseline(subjNormalBets,:)), 'o')
%   plot([0, 1; 0, 1], ...
%       [nanmean(betBinBaseline(subjNormalBets,:))+ste_betBaseline; nanmean(betBinBaseline(subjNormalBets,:))-ste_betBaseline],  '-b')
%   xlim([-.5, 1.5])
%   
%   subplot(2, 3, 6)
%   hold on
%   plot([0, 1], nanmean(betBinResponse(subjNormalBets,:)), 'o')
%   plot([0, 1; 0, 1], ...
%       [nanmean(betBinResponse(subjNormalBets,:))+ste_betResponse; nanmean(betBinResponse(subjNormalBets,:))-ste_betResponse],  '-b')
%   xlabel('Confidence')
%   xlim([-.5, 1.5])
% 
%   %saveas(gcf, 'basicPupilAnalyses.eps', 'epsc2')
%   close all
  
  
  % Next up: Individual differences:
  % 1) average baseline and response per subject (response residuals?)
  % correlated with:
  
  % -- 1) 
  
  
  [~,~,subjMeanRespResid] = regress(subMeanResponse',[ones(length(subMeanBaseline), 1), subMeanBaseline'])

  
  % perceptual bias       = beta_v1(:,2)
  % stability driven bias = beta_v1(:,3)
  
%   subplot(2,2,1)
%   plot( beta_v1(:,2), subMeanBaseline', 'o')
%   ylabel('Baseline diameter')
%   [RHO,PVAL]=corr(beta_v1(:,2), subMeanBaseline')
%   
%   
%   subplot(2,2,2)
%   plot(beta_v1(:,3), subMeanBaseline',  'o')
%   [RHO,PVAL]=corr(beta_v1(:,3), subMeanBaseline')
%   
%   subplot(2,2,3)
%   plot(beta_v1(:,2), subjMeanRespResid',  'o')
%   [RHO,PVAL]=corr(beta_v1(:,2), subjMeanRespResid)
%   
%   ylabel('Pupil response')
%   xlabel('Perceptual bias')
%  
%   subplot(2,2,4)
%   plot(beta_v1(:,3), subjMeanRespResid',  'o')
%   xlabel('Relevance dependent bias')
%   [RHO,PVAL]=corr(beta_v1(:,3), subjMeanRespResid)
%   saveas(gcf, 'basicPupIndDiffAnalysis.eps', 'epsc2')
%   close all

  
  % rest of the stats are trickier...
  % permutation testing for significance (multiple comparisons correction
  maskTime=1:150; % don't include baseline measurements
  toPlot=[1, 2, 3, 4];
  clear allIsSig
  for i = 1:length(toPlot)
      k=toPlot(i);
      input.data=respCoeffs(:,maskTime,k);
      input.clt=.05;
      input.numReps=10000;
      input.corrType='size';
      output=clusterBasedPermutationTesting(input);
      unique(output.corrP)
      allIsSig(:,i)=output.corrP<.05;
      p(i)=min(output.corrP);
      maxClustSize(i)=max(output.clusterSize);
  end
  significantEpochs = logical(allIsSig);
  
  % sloppy, but hardcoding the timing for now to get window of time in
  % seconds:
  timeInterval=2.5./150;
  maxContigTime=timeInterval.*maxClustSize;
  
  X1Labels=regLabels;
  for i = 1:length(toPlot)
      disp(sprintf('Pupil response coefficient: %s \n \n \t thresh time = %s  \n \t p-value= %s \n \n', ...
          X1Labels{i}, num2str(maxContigTime(i)),  num2str(p(i))))
  end

  
  
  clear effectPeak
  % make plot:
  for i =1:length(toPlot)
      effectPeak(i)=find(p(1,:,toPlot(i))==min(p(1,:,toPlot(i))));
  end
  
  
  % create a LOSO version of effectPeak, such that the pupil-predicted
  % learning figure can be independent of the more standard pupil analyses.
  
  
  
  % now for individual differences: DOES behavioral variance relate to
  % pupil variance?
  
  sharedVar=(zscore(beta_v1(:,2)) + zscore(beta_v1(:,3)));
  uniqueVar=(zscore(beta_v1(:,2))- zscore(beta_v1(:,3)));
  
  
  [RHO,PVAL]=corr(beta_v1(:,2), beta_v1(:,3));
  
  
  
  xMat=[ones(size(sharedVar)), zscore(beta_v1(:,2)), zscore(beta_v1(:,3))];
  subBase=nanmean(avg_subj_pup_traceMRN(:,1:3), 2);

  % Compute correlations between pupil and behavior at baseline:
  [baseIndPupCorr, base_pupCorrInt] = regress(subBase, [ones(length(uniqueID),1) sharedVar uniqueVar,]);

    % Compute correlations between pupil and behavior at other timepoints, controlling for baseline:
  clear indPupCorr b_int
  for i1 = 1:size(avg_subj_pup_traceMRN,2)
      X = zScoreX([ones(length(uniqueID),1) sharedVar uniqueVar, subBase]); %betaRU  aux_1 aux_2
      y = zscore(avg_subj_pup_traceMRN(:,i1));
      [indPupCorr(i1,:) b_int(i1,:,:)] = regress(y,X);%ridge(y,X,0.01);
  end
  
 
  
  
  for i = 1:4
      isSig_uncorrected(:,i)=sign(b_int(:,i,1).*b_int(:,i,2))==1;
  end
  
  sigDuration=sum(isSig_uncorrected).*timeInterval;
  input.numReps=500;
  sigDuration_perm=nan(input.numReps, length(sigDuration));
  for i = 1:input.numReps
       % permute columns 2&3 of explanatory matrix:
       permOrder=randperm(length(sharedVar));
       X=zScoreX(X);
       X(:,2)=X(permOrder,2);
       X(:,3)=X(permOrder,3);
       
       % loop through all pupil measurements and compute correlation:
       for i1 = 1:size(avg_subj_pup_traceMRN,2)
           y = zscore(avg_subj_pup_traceMRN(:,i1));
           [~, b_int_perm(i1,:,:)] = regress(y,X);%ridge(y,X,0.01);
       end
       
       
       % compute significant epochs:
       for i1 = 1:4
           isSig_uncorrected_perm(:,i1)=sign(b_int_perm(:,i1,1).*b_int_perm(:,i1,2))==1;
       end
       
       % 
       sigDuration_perm(i,:)=sum(isSig_uncorrected_perm).*timeInterval;

  end
  
 sharedP= (1- sum(sigDuration(2)> sigDuration_perm(:,2))./input.numReps);
 uniqueP= (1- sum(sigDuration(3)> sigDuration_perm(:,3))./input.numReps);


 % SHARED effect for individual differences (overall perceptual bias &
 % relevance modulation of bias)

% time points that include significant effects for individual differences regressors: 
 sharedSigTimes=isSig_uncorrected(:,2).*  ( sharedP<.05 );
 uniqueSigTimes=isSig_uncorrected(:,3).*   (uniqueP<.05 );

 disp(sprintf('Shared effect of perceptual bias and relevance modulation thereof:  \n \n \t thresh time = %s  \n \t p-value= %s \n \n', ...
    num2str(sigDuration(2)),  num2str(sharedP)))
 
 disp(sprintf('Unique effects [PE - relevance*PE]:  \n \n \t thresh time = %s  \n \t p-value= %s \n \n', ...
    num2str(sigDuration(3)),  num2str(uniqueP)))

  %% Identify peak times of pupil effects... use these for pupil-predicted perceptual bias figures
  
  
  clear peakResponse
  for i = 1:size(respCoeffs, 3)
      ind=find(nanmean(respCoeffs(:,:,i))==max(nanmean(respCoeffs(:,:,i))));
      peakResponse(:,i)=respCoeffs(:,ind, i)
  end
  
  [r, p]=corr([beta_v1, peakResponse])
  
  
  [I,J]=find(p<.05)
  
  clear B BINT
  for i =1: size( respCoeffs(:,:, 1), 2)
      [B(i,:),BINT(i,:,:)] =  regress(respCoeffs(:,i, 1), xMat);
  end
  
  hold on
  plot([0, 150], [0, 0], '--k')
  
  clear a
  for i = 2:size(BINT, 2);
      a(i)=plot(B(:,i), '-', 'color', cbColors(i+1,:));
      plot(BINT(:,i,1), '--', 'color', cbColors(i+1,:));
      plot(BINT(:,i,2), '--', 'color', cbColors(i+1,:));
  end
  legend(a(2:end), 'basePriorUse', 'stabilityAdjustment')
  title('correlation with CP pup effect')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            pupil-predicted bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

             
% Can we use pupil diameter to predict perceptual bias? Within subjects?
% across subjects?


% MRN changed this so that we use the time of peak pupil response, rather
% than a time selected according to our previous analysis. 

time            = linspace(0,150.0/60.0,150);
peakTimes       = ones(4,1).*time(respTime);           %time(effectPeak);
useResTimepoints= true; % if true, we will compute things on residual. 
                        % If false, predictions will be based on all other
                        % timepoints.
  
% now lets decomponse individual differences and within subject differencs:
subID_dummy=[];
subSoundResponse=zeros(size(nofirst_PredErr));
subBaseResponse=zeros(size(nofirst_PredErr));
subConfResponse=zeros(size(nofirst_PredErr));
subResResponse=zeros(size(nofirst_PredErr));
mcSoundResponse=zeros(size(nofirst_PredErr));
mcBaseResponse=zeros(size(nofirst_PredErr));
mcConfResponse=zeros(size(nofirst_PredErr));
mcResResponse=zeros(size(nofirst_PredErr));




for i = 1:max(noFirst_subjID)
    sel=i==noFirst_subjID;

    % Get average pupil diameter at various timepoints:
    meanSoundResponse(i)=avg_subj_pup_traceMRN(i,effectPeak(1));
    meanBaseResponse(i)=nanmean(avg_subj_pup_traceMRN(i,1:3));
    meanConfResponse(i)=avg_subj_pup_traceMRN(i,effectPeak(3));
    meanResResponse(i)=avg_subj_pup_traceMRN(i,effectPeak(4));
         
 
    % Extract residuals from peak timepoints
    subID_dummy(:,end+1)=sel';
    mcSoundResponse(sel)   = zScoreX(respResidual(sel, effectPeak(1)));
    mcBaseResponse(sel)    = zScoreX(baseResidual(sel)');
    mcConfResponse(sel)    = zScoreX(respResidual(sel, effectPeak(3)));
    mcResResponse(sel)     = zScoreX(respResidual(sel, effectPeak(4)));

end

mcSoundResponse(~isfinite(mcSoundResponse))=0;
mcBaseResponse(~isfinite(mcBaseResponse))=0;
mcConfResponse(~isfinite(mcConfResponse))=0;
mcResResponse(~isfinite(mcResResponse))=0;


% For average responses, regress out baseline diameter:
xes=[ones(size(meanBaseResponse))', meanBaseResponse']
[B,BINT,R_soundResponse] = regress(meanSoundResponse',xes)
[B,BINT,R_confResponse] = regress(meanConfResponse',xes)
[B,BINT,R_resResponse] = regress(meanResResponse',xes)

% Place subject average pupil response, after accounting for baseline, in
% containers for regression:
for i = 1:max(noFirst_subjID)
    sel=i==noFirst_subjID;
    subSoundResponse(sel)  = R_soundResponse(i);
    subBaseResponse(sel)   = meanBaseResponse(i);
    subConfResponse(sel)   = R_confResponse(i);
    subResResponse(sel)    = R_resResponse(i);
end


% line  1) basic terms
%       2) trial by trial terms
%       3) individual differences terms
%       4) interactions

% NOTE: kamesh changed code such that:
% allTrainingSpatialBias is now only computed for non day 1 probe trials. 

% normalize pupil measures per subject to put all subjects on the same
% footing and de-emphasize outliers. 

% Goals: 
% 1) create exact version of model-based analysis used from previous figure
% 2) add pupil terms to that model
% 3) compute likelihood ratio stats and AIC
% 4) recreate figure with new model

% mcSpaceTerm = centerError.


% Random subject effects:

subSpecPE_terms=subID_dummy.*repmat(nofirst_PredErr, 1,  size(subID_dummy,2));
sepRelTerms=subID_dummy.*repmat(nofirst_PredErr.*nofirst_Rel, 1,  size(subID_dummy,2));
sepBetTerms=subID_dummy.*repmat(nofirst_PredErr.*nofirst_Bet, 1,  size(subID_dummy,2));
sepStabTerms=subID_dummy.*repmat(nofirst_PredErr.*nofirst_Stab, 1,  size(subID_dummy,2));
sepCETerms=subID_dummy.*repmat(nofirst_CentErr, 1,  size(subID_dummy,2));
sepSpatialBiasTerms=subID_dummy.*repmat(noFirst_spaceBias, 1,  size(subID_dummy,2));

baseTerms    =[ones(size(nofirst_PredErr)),  nofirst_PredErr,noFirst_spaceBias, nofirst_CentErr, noFirst_spaceBias] ;
fixedModTerms=[baseTerms, nofirst_Stab, nofirst_Rel, nofirst_Bet] ;
mixedModTerms=[subID_dummy, subSpecPE_terms, sepStabTerms, sepRelTerms, sepBetTerms,sepCETerms, sepSpatialBiasTerms];

pupTerms=[nofirst_PredErr.*mcSoundResponse, nofirst_PredErr.*mcBaseResponse, nofirst_PredErr.*mcConfResponse, ...
    nofirst_PredErr.*zscore(subSoundResponse), nofirst_PredErr.*zscore(subBaseResponse), nofirst_PredErr.*zscore(subConfResponse), ...
    ] ;

simplePupTerms=[nofirst_PredErr.*mcResResponse, nofirst_PredErr.*zscore(subResResponse), nofirst_PredErr.*  zscore(mcResResponse.*zscore(subResResponse))] ;
sel=nofirst_TrialNo>trialsToDrop;

[B,BINT, r] = regress(nofirst_PercErr(sel),baseTerms(sel,:));
LL_fixedMod=sum(log(normpdf(r, 0, nanstd(r))))

[B,BINT, r] = regress(nofirst_PercErr(sel), [baseTerms(sel,:), pupTerms(sel,:)]);
LL_fixedPlusPup=sum(log(normpdf(r, 0, nanstd(r))))

[B,BINT, r] = regress(nofirst_PercErr(sel), [baseTerms(sel,:), simplePupTerms(sel,:)]);
LL_fixedPlusSparsePup=sum(log(normpdf(r, 0, nanstd(r))))

[B,BINT, r] = regress(nofirst_PercErr(sel), [fixedModTerms(sel,:)]);
LL_fixedWModBasedTerms=sum(log(normpdf(r, 0, nanstd(r))))

[B,BINT, r] = regress(nofirst_PercErr(sel), [fixedModTerms(sel,:), pupTerms(sel,:)]);
LL_fixedWModBasedTermsPlusPup=sum(log(normpdf(r, 0, nanstd(r))))

[B,BINT, r] = regress(nofirst_PercErr(sel), [fixedModTerms(sel,:), simplePupTerms(sel,:)]);
LL_fixedWModBasedTermsPlusSparsePup=sum(log(normpdf(r, 0, nanstd(r))))

[B,BINT, r1] = regress(nofirst_PercErr(sel), [mixedModTerms(sel,:)]);
LL_randModBasedTerms=sum(log(normpdf(r1, 0, nanstd(r1))))

[B,BINT, r] = regress(nofirst_PercErr(sel), [mixedModTerms(sel,:), pupTerms(sel,1:3)]);
LL_randModBasedTerms_wPup=sum(log(normpdf(r, 0, nanstd(r))))

[B,BINT, r2] = regress(nofirst_PercErr(sel), [mixedModTerms(sel,:), simplePupTerms(sel,[1])]);
LL_randModBasedTerms_wSparsePup=sum(log(normpdf(r2, 0, nanstd(r2))))


% likelihood ratio tests for adding pupil terms to each model:
[h,pValue,stat,cValue] = lratiotest(LL_fixedPlusPup,LL_fixedMod,6)
[h,pValue,stat,cValue] = lratiotest(LL_fixedWModBasedTermsPlusPup,LL_fixedWModBasedTerms,6)
[h,pValue,stat,cValue] = lratiotest(LL_randModBasedTerms_wPup,LL_randModBasedTerms,3)


% likelihood ratio tests for adding residual pupil terms to each model:
[h,pValue1,stat1,cValue1] = lratiotest(LL_fixedPlusSparsePup,LL_fixedMod,3)
[h,pValue2,stat2,cValue2] = lratiotest(LL_fixedWModBasedTermsPlusSparsePup,LL_fixedWModBasedTerms,3)
[h,pValue3,stat3,cValue3] = lratiotest(LL_randModBasedTerms_wSparsePup,LL_randModBasedTerms,1)

LR_isSig=[pValue1, pValue2, pValue3]<.05;
LR_isReallySig=[pValue1, pValue2, pValue3]<.0001;;


LikeRatio_stats=[stat1, cValue1, pValue1; stat2, cValue2 pValue2; stat3, cValue3, pValue3 ]





% Get stats for coefficients in simplest regression model:
coefStats=regstats(nofirst_PercErr(sel),[baseTerms(sel,:), simplePupTerms(sel,:)])
ind=size(baseTerms(sel,:), 2)+2;

% rows: within, between, interaction
% columns: beta, t-stat, pvalue:
coeffStats=[coefStats.tstat.beta(ind:end), coefStats.tstat.t(ind:end), coefStats.tstat.pval(ind:end) ]




% Get weights from weighted regression:
if ~useResTimepoints
    [Bw BINTw] = regressW_mike(nofirst_PercErr(sel),nofirst_expStd(sel), [baseTerms(sel,:), pupTerms(sel,:)]);
    noPE_mat=[normSoundResponse, normBaseline, normBetResponse, ...
        zscore(subSoundResponse),zscore(subBaseResponse), zscore(subConfResponse), ...
        ] ;
    
    [BIC1, AIC1]=computeBIC([LL_fixedMod, LL_fixedPlusPup,  LL_fixedWModBasedTerms, LL_fixedWModBasedTermsPlusPup, LL_randModBasedTerms, LL_randModBasedTerms_wPup ], ...
        [size(baseTerms,2), size(baseTerms,2)+3, size(fixedModTerms, 2), size(fixedModTerms, 2)+3, size(mixedModTerms, 2), size(mixedModTerms, 2)+1], ...
        [repmat(sum(sel), 1, 6)] )
    
else
    
    % this is SO slow:
    %[Bw BINTw] = regressW_mike(nofirst_PercErr(sel),nofirst_expStd(sel), [baseTerms(sel,:), simplePupTerms(sel,:)]);
    [Bw BINTw] = regress(nofirst_PercErr(sel), [baseTerms(sel,:), simplePupTerms(sel,:)])    
    noPE_mat=[mcResResponse, zscore(subResResponse),  zscore(mcResResponse.*zscore(subResResponse))] ;

    
    
        statsForCoeffs = regstats(nofirst_PercErr(sel), [baseTerms(sel,:), simplePupTerms(sel,:)])    

    
    
    [BIC1, AIC1]=computeBIC([LL_fixedMod, LL_fixedPlusSparsePup,  LL_fixedWModBasedTerms, LL_fixedWModBasedTermsPlusSparsePup, LL_randModBasedTerms, LL_randModBasedTerms_wSparsePup ], ...
        [size(baseTerms,2), size(baseTerms,2)+2, size(fixedModTerms, 2), size(fixedModTerms, 2)+2, size(mixedModTerms, 2), size(mixedModTerms, 2)+1], ...
        [repmat(sum(sel), 1, 6)] )

end


X1Labels={'within', 'between', 'interaction'};
nineFiveHalfWidth=Bw-BINTw(:,1);
% NOTE: regstats adds a column of ones... so the coefficients are
% incremented by 1. 
coefNumbers=size(baseTerms(sel,:), 2)+1  : size(baseTerms(sel,:), 2)+size(simplePupTerms,2);
for i = coefNumbers
disp(sprintf('Base model coefficient values: %s \n \n \t Mean/95CI beta= %s / %s \n \t t-stat = %s \n \t p-value= %s', ...
    X1Labels{coefNumbers==i}, num2str(Bw(i)), num2str(nineFiveHalfWidth(i)), num2str(statsForCoeffs.tstat.t(i+1)), statsForCoeffs.tstat.pval(i+1)))
end




% predictor list: 
% 1) intercept
% 2) prediction error
% 3) training bias
% 4) center error
% 5) CP-response 
% 6) rel (baseline) response
% 7) bet (end) response
% 8)  CP- sub response 
% 9)  rel (baseline) sub response
% 10) bet (end) sub response

% compute predicted prior weights based on pupil response:


% Bw(2) % fixed effect of prediction error on perceptual error:

pupPredicted_bias= noPE_mat(sel,:)*Bw((size(baseTerms,2)+1):end);
hiBin=pupPredicted_bias>prctile(pupPredicted_bias, 75);
loBin=pupPredicted_bias<prctile(pupPredicted_bias, 25);

allPredBias=[nanmean(pupPredicted_bias(hiBin)), ...
    nanmean(pupPredicted_bias(~hiBin& ~loBin)), ...
    nanmean(pupPredicted_bias(loBin))]; 
    
selPredErr=nofirst_PredErr(sel);
selPercErr=nofirst_PercErr(sel);


%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               REVIEWER REQUESTED ANALYSES (KAMESH JAN 2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% How many trials were rejected due to blinks?

subj_broken_fix = zeros(length(uniqueIDs),1);
for i0 = 1:length(uniqueIDs)
    selID = probeSubjID == uniqueIDs(i0);
    aux_is_fixed = probeIsFixed(selID) == 0;
    subj_broken_fix(i0) = sum(aux_is_fixed)*1.0/length(aux_is_fixed);
    
end



% What were the median RTs for estimation and bets

subj_medianRT_est = zeros(length(uniqueIDs),1);
for i0 = 1:length(uniqueIDs)
    selID = fixSubjID == uniqueIDs(i0);
    aux_estTime = fixTimeToEstimate(selID);
    subj_medianRT_est(i0) = median(aux_estTime);
    
end

subj_medianRT_bet = zeros(length(uniqueIDs),1);
for i0 = 1:length(uniqueIDs)
    selID = fixSubjID == uniqueIDs(i0);
    aux_betTime = fixTimeToBet(selID);
    subj_medianRT_bet(i0) = median(aux_betTime(aux_betTime>0));
    
end


subj_medianRT_predict = zeros(length(uniqueIDs),1);
for i0 = 1:length(uniqueIDs)
    selID = fixSubjID == uniqueIDs(i0);
    aux_predTime = fixTimeToPredict(selID);
    subj_medianRT_predict(i0) = median(aux_predTime(aux_predTime>0));
    
end

