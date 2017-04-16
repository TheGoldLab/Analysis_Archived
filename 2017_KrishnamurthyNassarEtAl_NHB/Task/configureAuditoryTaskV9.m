function [tree, list] = configureAuditoryTaskV9(logic, avConfig, options)
% File for initializing the Auditory Task. This uses the AuditoryLogic
% and AuditoryAV templates which extend Ben's templates for Matt's Helicopter task. These
% templates are meant to collect common code for a broad range of
% predictive inference tasks.
%
% logic -- object that contains most of the mathematical settings and
% variables to store the state of the task
% avConf  -- object that contains the audio-visual settings and objects
% --Kamesh 25/10/12

sc=dotsTheScreen.theObject;
if(options.fullScreen && ~options.usingCRT)
    sc.reset('displayIndex', 1);
elseif(options.usingCRT)
    sc.reset('displayIndex', 2);
end

if nargin < 1 || isempty(logic)
    logic = AuditoryLogic9(options.randSeed);
    logic.blockStds = options.blockStds;
    logic.blockHazards = options.blockHazards;
    logic.trialsPerBlock = options.trialsPerBlock;
    logic.nBlocks = length(options.blockStds);
    logic.safetyTrials = options.safetyTrials;
    logic.sampleHazard = options.sampleHazard;
    logic.pWagerSample = options.pWagerSample;
    logic.minWagerRunLength = options.minWagerRunLength;
end
    

if nargin < 2 || isempty(avConfig)
    avConfig = AuditoryAV9();
    avConfig.wagerWindowSize = options.wagerWindowSize;
    avConfig.logic = logic;
    avConfig.hrirID = options.hrir;
    avConfig.dataPath = ['~/Data/WavFiles/subj' num2str(options.hrir) '/'];
end


% Configuring inputs -- you might want to use the mouse by default
if nargin < 3 || isempty(options) || ~isstruct(options)
    options.isKeyboardTrigger = false;
    options.triggerKeyName = 'KeyboardT';
    options.triggerMaxWait = 5*60;
    options.isPositionMapping = false;
end

if ~isfield(options, 'keyboardKeys')
    options.keyboardKeys.left = 'KeyboardLeftArrow';
    options.keyboardKeys.right = 'KeyboardRightArrow';
    options.keyboardKeys.commit = 'KeyboardSpacebar';
    options.keyboardKeys.info = 'KeyboardI';
    options.keyboardKeys.abort = 'KeyboardQ';
end


%set the wager trial sampling list
avConfig.logic.wagerTrialList = options.wagerTrialList;


% Configuring the eyetracker
if ~isfield(options, 'eyeTracker')
    options.eyeTracker.isEyeTracking = true;
    options.eyeTracker.sampleFrequency = 60;
    %fixation window for the Tobii eyetracker (relative coordinates)
    options.eyeTracker.fixWindow = [0.3 0.3 0.6 0.6];
end


% Data logging
topsDataLog.flushAllData();
t = topsDataLog.theDataLog;
t.fileWithPath = options.fileName;
% record some global info about the task
if(~(strcmp(options.taskType, 'trainNoFeedback') || ...
        strcmp(options.taskType, 'trainFeedback')))
    t.logDataInGroup(options.trialsPerBlock,'trialsPerBlock');
    t.logDataInGroup(options.blockStds,'blockStds');
    t.logDataInGroup(options.blockHazards,'blockHazards');
else
    t.logDataInGroup(options.trainIterations,'trainIterations');
end
t.logDataInGroup(options.hrir,'hrirNum');
t.logDataInGroup(options.subject,'subject');
t.logDataInGroup(options.taskType,'taskType');
t.logDataInGroup(options.session,'session');
t.logDataInGroup(options.randSeed,'randSeed');
t.logDataInGroup(options.wagerTrialList,'wagerTrialList');




%% TOPS containers for data and objects
list = topsGroupedList();
list{'logic'}{'object'} = logic;
list{'audio-visual'}{'object'} = avConfig;
list{'options'}{'struct'} = options;

%% -----------------------------INPUT------------------------------------
% Try to use a mouse and use a keyboard if not available.
usingMouse = options.usingMouse;

if usingMouse
    % These settings are very specific for the setup I'm using in the
    % psychophysics lab -- my laptop connected to the monitor and a
    % microsoft external mouse

    opt.Manufacturer = 'Microsoft';
    m = dotsReadableHIDMouse(opt);
    %m = dotsReadableHIDMouse;
    m.isExclusive = 1;
    m.isAutoRead = 1;
    
    m.flushData;
    m.initialize();
    
    
    % undefine any default events
    IDs = m.getComponentIDs();
    for ii = 1:numel(IDs)
        m.undefineEvent(IDs(ii));
    end
    %Define a mouse button press event
    m.defineEvent(3, 'left', 0, 0, true);
    m.defineEvent(4, 'right', 0, 0, true);
    ui = m;
    getInput = @getMouseInput;
    %store the mouse separately in case we need to use it
    list{'input'}{'mouse'} = m;
    
else
    
    
    kb = dotsReadableHIDKeyboard();
    
    %DEBUG -- for psychophysics lab keyboard
    trials = 0;
    maxTrials = 500;
    while(~strcmp(kb.deviceInfo.Product, 'Apple Keyboard') ...
                                            && trials < maxTrials)
                                        
         kb =  dotsReadableHIDKeyboard();                                                            
    end
    
    
    
    IDs = kb.getComponentIDs();
    for ii = 1:numel(IDs)
        kb.undefineEvent(IDs(ii));
    end
    isDown = strcmp({kb.components.name}, 'KeyboardDownArrow');
    downKey = kb.components(isDown);
    kb.defineEvent(downKey.ID,'pressed',0,0,true);
    isRight = strcmp({kb.components.name}, 'KeyboardRightArrow');
    rightKey = kb.components(isRight);
    kb.defineEvent(rightKey.ID,'right',0,0,true);
    isLeft = strcmp({kb.components.name}, 'KeyboardLeftArrow');
    leftKey = kb.components(isLeft);
    kb.defineEvent(leftKey.ID,'left',0,0,true);
    isUp = strcmp({kb.components.name}, 'KeyboardUpArrow');
    upKey = kb.components(isUp);
    kb.defineEvent(upKey.ID,'up',0,0,true);
    
    ui = kb;
    getInput = @getKbInput;
    %store the keyboard separately in case we need to use it
    list{'input'}{'keyboard'} = kb;
end


%DEBUG (for wager) -- testing only. Remove in psychophysics comp

kb = dotsReadableHIDKeyboard();
IDs = kb.getComponentIDs();
for ii = 1:numel(IDs)
    kb.undefineEvent(IDs(ii));
end
isDown = strcmp({kb.components.name}, 'KeyboardDownArrow');
downKey = kb.components(isDown);
kb.defineEvent(downKey.ID,'pressed',0,0,true);
isRight = strcmp({kb.components.name}, 'KeyboardRightArrow');
rightKey = kb.components(isRight);
kb.defineEvent(rightKey.ID,'right',0,0,true);
isLeft = strcmp({kb.components.name}, 'KeyboardLeftArrow');
leftKey = kb.components(isLeft);
kb.defineEvent(leftKey.ID,'left',0,0,true);
isUp = strcmp({kb.components.name}, 'KeyboardUpArrow');
upKey = kb.components(isUp);
kb.defineEvent(upKey.ID,'up',0,0,true);

%store the keyboard separately in case we need to use it
list{'input'}{'keyboard'} = kb;



avConfig.ui = ui;
list{'input'}{'controller'} = ui;
%list{'input'}{'mapping'} = uiMap;

%% Set up the eyetracker   (Tobii 60 XL for V9)

if options.eyeTracker.isEyeTracking
    %disp('Initializing Tobii Eye Tracker\n');
    %tetio_init();
    % Set to tracker ID to the product ID of the tracker you want to connect to.
    %trackerId = 'XL060-31500295.local.';
    %fprintf('Connecting to tracker "%s"...\n', trackerId);
    %tetio_connectTracker(trackerId)
	
    %currentFrameRate = tetio_getFrameRate;
    %fprintf('Eye Tracker Frame Rate: %d Hz.\n', currentFrameRate);
    
else
    options.eyeTracker.isEyeTracking = false;
end

%% Outline the structure of the exeriment with topsRunnable objects
% visualize the structure with tree.gui()
% run the experiment with tree.run()

% Define nodes of the topsTree tree->session->block->trial
% "tree" is the start point for the whole experiment
tree = topsTreeNode('open/close screen');
tree.iterations = 1;
tree.startFevalable = {@initialize, list};
tree.finishFevalable = {@terminate, list};

% "instructions" is a branch of the tree with an instructional slide show
instructions = topsTreeNode('instructions');
instructions.iterations = 1;
tree.addChild(instructions);

% "session" is a branch of the tree with the task itself
session = topsTreeNode('session');
session.iterations = logic.nBlocks;
session.startFevalable = {@startSession, logic};
session.finishFevalable = {@finishSession, logic};
tree.addChild(session);

% define a 'Block' node : child of session.
block = topsTreeNode('block');
block.iterations = logic.trialsPerBlock;
block.startFevalable = {@startBlock, list};
block.finishFevalable = {@finishBlock, list};
session.addChild(block);


% 'trial' node
% ConcurrentComposite is a collection of topsConcurrent objects which
% can run concurrently. 'trial' runs by running its children.
trial = topsConcurrentComposite('trial');
block.addChild(trial);
% StateMachine defines the flow control within a single trial
trialStates = topsStateMachine('trial states');
trial.addChild(trialStates);


%Store the 'tree' object in our topsList
list{'outline'}{'tree'} = tree;


%% Organize the flow through each trial
% The trial state machine will respond to user input commands
% and control timing.


switch options.taskType
    case 'fullTask'

        avConfig.taskType = 1;   %objects displayed depend on task type
        
        % In this version you fixate and then the window comes up and then
        % the sound plays
        
        if options.eyeTracker.isEyeTracking
            states = { ...
            'name'         'next'      'timeout'      'entry'                         'exit'                    'input'; ...
            'predict'      'fixate1'      0           {@setupPrediction avConfig}     {}                        {getInput avConfig};...
            'fixate1'      'outcome'             20          {@setupPacedFixation avConfig}  {@setupFixation avConfig} {@checkPacedFixation avConfig};...
            'outcome'      'checkNextSt'  0           {@setupOutcome avConfig}        {@pause 0.1}              {@isNextStateFixate avConfig};...
            'fixate2'      'estimate'     2.5         {}                              {@getEyeData avConfig}    {};...
            'estimate'     'wager'        0           {}                              {}                        {getInput avConfig};...
            'wager'        'success'      0           {@setupWager avConfig}          {}                        {@getWager avConfig};...
            'success'      ''             0.0         {@doSuccess avConfig}           {}                        {}; ...
            };
            
        else
             error('This task needs the Tobii 60 XL eyetracker. Revert to a previous version otherwise');
        end
        
        
%         if options.eyeTracker.isEyeTracking
%             states = { ...
%             'name'         'next'      'timeout'      'entry'                         'exit'                    'input'; ...
%             'predict'      'fixate1'      0           {@setupPrediction avConfig}     {}                        {getInput avConfig};...
%             'fixate1'      'outcome'      50          {@setupFixation avConfig}       {'tetio_readGazeData'}    {@checkPacedFixation avConfig};...
%             'outcome'      'checkNextSt'  0           {@setupOutcome avConfig}        {@pause 0.2}              {@isNextStateFixate avConfig};...
%             'fixate2'      'estimate'     2.5         {@setupFixation avConfig}       {@getEyeData avConfig}    {};...
%             'estimate'     'wager'        0           {}                              {}                        {getInput avConfig};...
%             'wager'        'success'      0           {@setupWager avConfig}          {}                        {@getWager avConfig};...
%             'success'      ''             0.0         {@doSuccess avConfig}           {}                        {}; ...
%             };
%             
%         else
%              error('This task needs the Tobii 60 XL eyetracker. Revert to a previous version otherwise');
%         end
        
        %note : isNextStateFixate checks if this is a wager trial.
        %for wager trials fixate2 is entered otherwise 'estimate' is
        %enetered
        
        
    case 'estimateTask'
        
        %In this version of the task, there is no prediction; otherwise,
        %everthing else is identical
        
        avConfig.taskType = 2;
        
        if options.eyeTracker.isEyeTracking
            states = { ...
            'name'         'next'      'timeout'      'entry'                         'exit'                    'input'; ...
            'fixate1'      'outcome'      20          {@setupPacedFixation avConfig}  {@setupFixation avConfig} {@checkPacedFixation avConfig};...
            'outcome'      'checkNextSt'  0           {@setupOutcome avConfig}        {@pause 0.1}              {@isNextStateFixate avConfig};...
            'fixate2'      'estimate'     2.5         {}                              {@getEyeData avConfig}    {};...
            'estimate'     'wager'        0           {}                              {}                        {getInput avConfig};...
            'wager'        'success'      0           {@setupWager avConfig}          {}                        {@getWager avConfig};...
            'success'      ''             0.0         {@doSuccess avConfig}           {}                        {}; ...
            };
            
        else
             error('This task needs the Tobii 60 XL eyetracker. Revert to a previous version otherwise');
        end
        
        
    case 'trainFeedback'
        if(~ options.eyeTracker.isEyeTracking)
           error('this version of the task need Eye tracking for training sessions'); 
        end
        avConfig.taskType = 4;
        logic.isChangePointTask = false;
        logic.nBlocks = options.trainIterations;
        session.iterations = logic.nBlocks;
        fixedOutcomes = zeros(logic.nBlocks, options.trialsPerBlock);
        aux = 0:15:180;
        for i1 = 1:logic.nBlocks
            idx = randperm(options.trialsPerBlock);
            fixedOutcomes(i1,:) = aux(idx);
        end
        fixedOutcomes = fixedOutcomes';
        logic.fixedOutcomes = fixedOutcomes(:);
        states = { ...
            'name'          'next'      'timeout'      'entry'                      'exit'                          'input'; ...
            'fixate1'      'outcome'      20          {@setupPacedFixation avConfig} {@setupFixation avConfig}    {@checkPacedFixation avConfig};...
            'outcome'       'fixate2'     0           {@setupOutcome avConfig}      {}                            {}; ...
            'fixate2'       'checkNextSt' 2.5         {}                            {@checkTrainingFixation avConfig}        {};...
            'checkNextSt'   'estimate'    0           {}                            {}                            {@isNextStateFailure avConfig}
            'estimate'      'success'     0           {}                            {}                            {getInput avConfig};...
            'success'       ''            0.8         {@doSuccess avConfig}         {}                            {}; ...
            'badFixation'   'failure'     0           {}                            {}                            {};...
            'failure'       ''            0           {@doFailure avConfig}         {@pause 0.2}                  {};...
            };
        %@waitForUser avConfig  -- removed this from the initial state
    case 'trainNoFeedback'
        avConfig.taskType = 3;
        logic.nBlocks = options.trainIterations;
        session.iterations = logic.nBlocks;
        logic.isChangePointTask = false;
        fixedOutcomes = zeros(logic.nBlocks, options.trialsPerBlock);
        aux = 0:15:180;
        for i1 = 1:logic.nBlocks
            idx = randperm(options.trialsPerBlock);
            fixedOutcomes(i1,:) = aux(idx);
        end
        fixedOutcomes = fixedOutcomes';
        logic.fixedOutcomes = fixedOutcomes(:);
        states = { ...
            'name'          'next'      'timeout'      'entry'                      'exit'                          'input'; ...
            'start'        'outcome'      0.8           {}                              {}                         {};...
            'outcome'       'estimate'    0           {@setupOutcome avConfig}      {@pause 1.1}                            {}; ...
            'estimate'      ''            0           {}                            {}                            {getInput avConfig};...
            };
    otherwise
        error('Enter a valid task type');
end


trialStates.addMultipleStates(states);
trialStates.startFevalable = {@startTrial list};
trialStates.finishFevalable = {@finishTrial list};

avConfig.list = list;




%% Helper functions
function startTrial(list)
logic = list{'logic'}{'object'};
logic.startTrial();  %This will generate outcome for the trial

av = list{'audio-visual'}{'object'};
av.prepareSoundFile();  %load the sound file before trial begins

ui = list{'input'}{'controller'};
ui.flushData();

% ************ DEBUG (DBGTOB)-- for Tobii eye tracker ******************

% place holder for eye tracker buffer flush before the start of the trial

%************** DEBUG END *********************************************


function finishTrial(list)
logic = list{'logic'}{'object'};
logic.finishTrial();

function startBlock(list)
%Draw the background and wait for user to proceed
av = list{'audio-visual'}{'object'};
av.waitForBlockStart();
av.numBrokenFix = 0;
av.numWagerTrial = 0;

%record time from the start of the block
av.blockStartTime = 0;
av.blockStartTime = tic;

%DEBUG -- let us start tracking at the beginning of every block
% and stop at the end.
tetio_startTracking;


logic = list{'logic'}{'object'};
logic.startBlock();

function finishBlock(list)

%place holder for finish blcok function
%DEBUG -- let us start tracking at the beginning of every block
% and stop at the end.
tetio_stopTracking;

av = list{'audio-visual'}{'object'};
logic = list{'logic'}{'object'};


function initialize(list)
%place holder for tree initialization
av = list{'audio-visual'}{'object'};
av.initialize();

function terminate(list)
%place holder for tree termination
av = list{'audio-visual'}{'object'};
av.terminate();
topsDataLog.writeDataFile();
opt = list{'options'}{'struct'};
fileName = opt.fileName;
data = topsDataLog.getSortedDataStruct;
save(fileName,'data','-append');

%DEBUG (Scoring Kamesh)
%display the total score at the end of the task
idx2 = strcmp('wagerTrialIdx',{data.group});
wagerTrialIdx = data(idx2);
wagerTrialIdx = cell2mat({wagerTrialIdx.item});
    
%need to correct for block index since trial no. gets set to
%one at the start of each block
idxTP = strcmp('trialsPerBlock',{data.group});
numTrialsPerBlock = data(idxTP);
numTrialsPerBlock = cell2mat({numTrialsPerBlock.item});
blockStartIdx = find(diff([500 wagerTrialIdx 1])<0);

idxScore = strcmp('score',{data.group});
score = data(idxScore);
score = cell2mat({score.item});
finalScores = score(blockStartIdx(2:end)-1);
fprintf('block scores are ---> \n');
fprintf('%d\n',finalScores);


