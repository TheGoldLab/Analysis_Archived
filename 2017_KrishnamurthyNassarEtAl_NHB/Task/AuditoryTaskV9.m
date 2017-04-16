

% This is an extension of the Auditory Task with the provisions for
% pupilometry measurements. This version of the task uses SnowDots ver 1.0
% and does not rely on the Psychophysics toolbox.

% This version of the task uses the Tobii eyetracker from Joe Kable's lab
% for better gaze tracking and more accurate pupil reports.

% Task types:
% 1- 'fullTask' - full task (predict, estimate, feedback)
% 2- 'estimateTask' - no prediction, just estimate, feedback (w/ CP structure)
% 3- 'trainNoFeedback' - meas. percp. accuracy w/o CP structure. No feedback
% 4- 'trainFeedback' - meas. percp. accuracy w/o CP structure. w/ feedback

% Typical order of doing these tasks:
% Day/Session 1 : subject does 4,2 & 1 (in that order)
% Day/Session 2 : subject does 3 & 1/2
% Day/Session 3 : subject does 3 & 2/1 (which was not done on Day2)

% -- Kamesh Krishnamurthy (kameshkk@gmail.com) U.Penn. May 2014


clear classes; disp('clearing workspace!');

%preTaskScreen;

% CRT/Local screen and fullscreen/window
options.fullScreen = true;
options.usingCRT = true;

%HRIR number used for the subject Use the best HRIR for eacmglCloseh subject
options.hrir = 1023;

%subject ID
options.subject = 'MK';

%task type & session number. Task type is variant of the task we are
%running on the subject
options.session = 'D4S3';
options.taskType = 'fullTask';

%if this is passive viewing 
%of cluster sounds, set this to true. You needn't make any other changes.
passive = false;  

%SET THE WAGER WINDOW SIZE
options.wagerWindowSize = 16;


% minimum number of outcomes between wager trials
options.minWagerRunLength = 8;     % ~5-6 seconds for the pupil to relax between wagers

%specify standard deviations for each block
%There will be one block for each std-dev.

stdevs =  [10 10 20 20];
rndIdx = randperm(length(stdevs));
stdevs = stdevs(rndIdx);
%stdevs =  [10];


if(strcmp(options.taskType,'fullTask') && passive)
    stdevs = [10];
end

%options.blockStds = stdevs(randperm(length(stdevs)));
options.blockStds = stdevs;

%specify the hazard rate (one element for each block)
hazardRate = 0.15;
options.blockHazards = [hazardRate hazardRate hazardRate hazardRate];
%options.blockHazards = [hazardRate];

% Create a vector of sampling probabilities
%pWagerSample = (1 - hazardRate).^(5:-1:0);
%pWagerSample(end-1) = 1;
%pWagerSample = [0.44 0.55 0.65 0.78 1.0 1.0];
pWagerSample = [0.55 0.64 0.67 0.78 0.95 1.0];
options.pWagerSample = pWagerSample;

%options.blockHazards = [0.15 0.15];

if(strcmp(options.taskType,'fullTask') && passive)
    options.blockHazards = [0.15];
end

%sample hazard
options.sampleHazard = 0.15;

%how many trials before a change-point
options.safetyTrials = 0; %5

%create list for sampling wager trials
numSamples = 50;
allAux=[];
for i = 1:numSamples;
aux = [1:6 1:6];
idx1 = randperm(length(aux));
allAux=[allAux aux(idx1)];
end
%we don't want a wager trial right at the start
allAux = [4 allAux];
if(strcmp(options.taskType,'fullTask') && passive)
    allAux = [50 50 50 50];
end

 
options.wagerTrialList = allAux;


%How many trials per block?
if((strcmp(options.taskType, 'trainNoFeedback') || ...
        strcmp(options.taskType, 'trainFeedback')))
    %how many training iterations?
    options.trainIterations = 6;     % no. of iterations for likelihood meas.
    options.trialsPerBlock = 195/15; %Do not change this
else
    strAux = options.session;
    if(strcmpi(strAux(1:2),'D1'))
       options.trialsPerBlock = 500;  % 500 on day1 and 600 on other days 
    else
       options.trialsPerBlock = 600;  % 500 on day1 and 600 on other days 
    end
    if(strcmp(options.taskType,'fullTask') && passive)
        options.trialsPerBlock = 80;
    end
end

%seed for the random number generator
%options.randSeed = 23232;

%rng shuffle;
%aux = rng;
%options.randSeed = aux.Seed;
%fprintf('Random number seed --> %d\n',aux.Seed);

a=clock;
options.randSeed=a(6);

% Mouse or keyboard
options.usingMouse = true; 


%path to data log file
logPath = '~/Dropbox/AuditoryTask/AuditoryTaskV9/Data/';
options.fileName = [logPath  options.subject ...
    options.session options.taskType '(' date() ').mat'];

if(length(options.blockStds) ~= length(options.blockHazards))
    error('blockStds and blockHazards must be of the same size');
end


%Use this for the version of the task w/ the wager
[tree list] = configureAuditoryTaskV9([], [], options);
%Use this for the version of the task w/o the wager
% [tree list] = configureAuditoryTaskV7a([], [], options);



tree.run();

