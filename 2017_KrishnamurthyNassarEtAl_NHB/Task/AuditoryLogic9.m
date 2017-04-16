classdef AuditoryLogic9 < handle
    % This class defines numeric and logical behaviors for the Auditory
    % task. The class structure is almost identical to PredInfAV which
    % contains the behaviors for Matt's predictive inference tasks.
    %
    % --------------------------------------------------------------
    % This has been modified  for the auditory task.
    % -Convention for angles : 180 (left) and 0 (right)
    % --Kamesh 25/10/2012
    
    % 9th April 2012:
    % Modifying the task so that only a few randomly selected trials are
    % real "trials", i.e. the subject has to do something; on other trials
    % the subject simply listens to the current outcome. The true trials
    % are selected randomly based on a fixed hazard (sampleHazard)
    
    properties
        % a name to identify this session
        name = '';
        
        % a time to identify this session
        time = 0;
        
        %which type of task is being run?
        % Task types:
        % 1- 'fullTask' - full task (predict, estimate, feedback)
        % 2- 'estimateTask' - no prediction, just estimate, feedback (w/ CP structure)
        % 3- 'trainNoFeedback' - meas. percp. accuracy w/o CP structure. No feedback
        % 4- 'trainFeedback' - meas. percp. accuracy w/o CP structure. w/ feedback
        taskType;
        
        %does the task have a change-point structure?
        isChangePointTask = true;
        
        % the number of blocks to run
        nBlocks = 1;
        
        % a standard deviation for each block
        blockStds = 5;
        
        % a hazard rate for each block
        blockHazards = 0.15;
        
        % a hazard for sampling the true trials
        sampleHazard = 0.25;
        
        % number of zero-hazard trials following a change trial
        safetyTrials = 1;
        
        % minimum run length of outcomes between wager trials
        minWagerRunLength = 8;
        wagerRunLength;
        
        % number of trials within each block
        trialsPerBlock = 10;
        
        % whether to shuffle the order of blockStds and blockHazards
        isBlockShuffle = false;
        
        % a list of outcomes to choose sequentially, instead of randomly
        fixedOutcomes = [];
        
        % seed for initializing random number generators
        randSeed = 0;
        
        %MOD FOR AUDITORY
        % maximum value for randomly chosen outcomes
        maxOutcome = 180;
        minOutcome = 0;
        
        % whether to reset subject's prediction each trial
        %MOD FOR AUDITORY -- the subject's prediction is reset
        isPredictionReset = true;
        
        % whether to limit subject's prediction each trial
        %MOD FOR AUDITORY -- the subject's prediction won't be limited
        isPredictionLimited = false;
        
        % whether the subject is currently allowed to make a prediction
        isPredictionActive = true;
        
        % arbitrary grand count of subject's "score" per block
        blockScore = 0;
        
        % arbitrary per-trial data to include in getStatus()
        trialData = [];
        
        % which element of fixedOutcomes was last picked
        fixedOutcomeIndex;
        
        % running count of good trials in each block
        blockCompletedTrials;
        
        % running count of good and bad trials in each block
        blockTotalTrials;
        
        %list of TACs that will be sampled as wager trials for each run
        wagerTrialList;
        wagerTrialNo = 1; %how many wager trials have been sampled so far?
        
        %are we sampling all TACs for wager?
        wagerTAC = 3;
        %wagerSampProb = 1/9*ones(1,9);
        
        pWagerSample;
        
        
    end
    
    properties (SetAccess = protected)
        
        % index of the current block
        currentBlock;
        
        % sequence of block numbers for choosing blockStds and blockHazards
        blockSequence;
        
        % whether the current trial is a change trial
        isChangeTrial;
        
        % running count of trials since the last change trial (1-based)
        steadyTrials;
        
        % remaining safety trials since the last change trial
        remainingSafety;
        
        % change hazard rate used in the current trial
        currentHazard;
        
        % std of the outcome distribution used in the current trial
        currentStd;
        
        % mean of the outcome distribution used in the current trial
        currentMean;
        
        % mean of the outcome distribution used in the previous trial
        lastMean;
        
        % pick from the outcome distribution for the current trial
        currentOutcome;
        
        % pick from the outcome distribution from the previous trial
        lastOutcome;
        
        % subject's pick for the current trial
        currentPrediction;
        
        % subject's pick from the previous trial
        lastPrediction;
        
        % whether the current trial is good and complete
        isGoodTrial;
        
        % whether the current trial is a true trial (post-decision wager)
        isWagerTrial;
        
        
        
        %don't have two wager trials within the same run
        wagersThisRun;
        
        % trial number within the block
        trialNo;
        
        
    end
    
    methods
        
        % Constructor takes no arguments.
        function self = AuditoryLogic9(randSeed)
            if nargin
                self.randSeed = randSeed;
            end
            self.startSession();
        end
        
        % Initialize random number functions from randSeed.
        function initializeRandomGenerators(self)
            % In the latest version of the task (V5 and later) the RNG
            % stuff is handled in AuditoryTask. I'm keeping this here just
            % in case we decide to use randseed.
            %rand('seed', self.randSeed);
            %randn('seed', self.randSeed);
            %rng(self.randSeed, 'twister');
        end
        
        % Set up for the first trial of the first block.
        function startSession(self)
            self.initializeRandomGenerators();
            
            % pick straight or shuffled sequence for block parameters
            if self.isBlockShuffle
                self.blockSequence = randperm(self.nBlocks);
            else
                self.blockSequence = 1:self.nBlocks;
            end
            
            % reset session totals
            self.fixedOutcomeIndex = 0;
            self.currentBlock = 0;
            self.blockCompletedTrials = 0;
            self.blockTotalTrials = 0;
            
            self.wagerTrialNo = 1;  %set the wagerTrial to 1
            
        end
        
        % Set up for a new block.
        function startBlock(self)
            % increment the block and choose block parameters
            self.currentBlock = self.currentBlock + 1;
            sequenceIndex = self.blockSequence(self.currentBlock);
            if(self.isChangePointTask)
                self.currentStd = self.blockStds(sequenceIndex);
                self.currentHazard = self.blockHazards(sequenceIndex);
            end
            
            % pick the first mean outcome at random
            if(self.isChangePointTask)
                self.pickChangeTrial();
            end
            
            % reset trial-by-trial state
            self.isPredictionActive = true;
            self.blockCompletedTrials = 0;
            self.blockTotalTrials = 0;
            %for change-point task
            self.isChangeTrial = false;
            self.steadyTrials = 0;
            self.remainingSafety = self.safetyTrials;
            self.currentOutcome = nan;
            self.currentPrediction = 0;
            self.lastMean = nan;
            self.lastOutcome = nan;
            self.lastPrediction = nan;
            
            self.trialNo = 0;
            
            %Initialize score to 0
            self.blockScore = 0;
            
            self.wagerRunLength = 1;
            
            %randomize the wagerTAC
            %self.wagerTrialNo = 1;   %DEBUG (Kamesh) this is to handle sampling with smaller blocks
            self.wagerTAC = self.wagerTrialList(self.wagerTrialNo);
            self.wagersThisRun = 0;
        end
        
        % Set up for a new trial.
        function startTrial(self)
            
            
            if(self.isChangePointTask)
                % remember several values from the previous trial
                self.lastMean = self.currentMean;
                self.lastOutcome = self.currentOutcome;
                self.lastPrediction = self.currentPrediction;
                
                
                % should this trial be a change trial?
                if self.remainingSafety > 0
                    % no, its still a safety trial
                    self.isChangeTrial = false;
                    self.remainingSafety = self.remainingSafety - 1;
                    self.steadyTrials = self.steadyTrials + 1;
                    
                else
                    
                    % maybe, toss a coin to find out
                    self.isChangeTrial = rand(1) < self.currentHazard;
                    if self.isChangeTrial
                        self.pickChangeTrial();  %pick a new mean
                        self.steadyTrials = 1;
                        self.remainingSafety = self.safetyTrials - 1;
                        
                        
                        %log the change-point trials --- 
                        topsDataLog.logDataInGroup(self.blockTotalTrials+1,...
                                              'cpTrialsNumber');
                        topsDataLog.logDataInGroup(1,'isCpTrial');
                        
                        
                        %begin new run for wagers
                        self.wagersThisRun = 0;
                        
                        %is this a wager trial?
                        if(self.wagersThisRun==0 && ...
                                 self.wagerTAC == self.steadyTrials && ...
                                  self.wagerRunLength > self.minWagerRunLength)
                             
                            %Now decide whether you actually want to sample this trial
                            % toss a coin with prob. pWagerSample(i) 
                            isSampled = rand(1) < self.pWagerSample(self.wagerTAC);
                            if(isSampled)
                                self.isWagerTrial = true;
                                fprintf('************** WAGER TRIAL ***************\n');
                                self.wagersThisRun = self.wagersThisRun + 1;
                                %log the trial number
                                topsDataLog.logDataInGroup(self.blockTotalTrials+1,...
                                              'wagerTrialIdx');           
                                %log what TAC is this wager trial
                                topsDataLog.logDataInGroup(self.wagerTAC,...
                                              'wagerTAC');
                                          
                                topsDataLog.logDataInGroup(self.wagerRunLength,...
                                              'wagerRunLength');
                                          
                                %pick wagerTAC from a randomized list
                                self.wagerTAC = self.wagerTrialList(self.wagerTrialNo);
                                self.wagerTrialNo = self.wagerTrialNo + 1;
                                self.wagerRunLength = 1;
                            else
                                self.isWagerTrial = false;
                            end    
                            
                        else
                            self.isWagerTrial = false;
                        end
                    else
                        self.steadyTrials = self.steadyTrials + 1;
    
                        %log the change-point trials
                        topsDataLog.logDataInGroup(0,'isCpTrial');
                        
                        %is this a wager trial?
                        if(self.wagersThisRun==0 && ...
                                 self.wagerTAC == self.steadyTrials && ...
                                 self.wagerRunLength > self.minWagerRunLength)
                             
                            %Now decide whether you actually want to sample this trial
                            % toss a coin with prob. pWagerSample(i) 
                            isSampled = rand(1) < self.pWagerSample(self.wagerTAC);
                            if(isSampled)
                                self.isWagerTrial = true;
                                fprintf('************** WAGER TRIAL ***************\n');                        
                                self.wagersThisRun = self.wagersThisRun + 1;
                            
                                %log the trial number
                                topsDataLog.logDataInGroup(self.blockTotalTrials+1,...
                                              'wagerTrialIdx');
                                                                
                                %log what TAC is this wager trial
                                topsDataLog.logDataInGroup(self.wagerTAC,...
                                              'wagerTAC');
                                          
                                topsDataLog.logDataInGroup(self.wagerRunLength,...
                                              'wagerRunLength');
                            %if(self.wagerTAC <= 7)
                            %self.wagerTAC = self.wagerTAC + 1;
                            %else
                            %    self.wagerTAC = 1;
                            %end
                            
                                %pick wagerTAC from a randomized list
                                self.wagerTAC = self.wagerTrialList(self.wagerTrialNo);
                                self.wagerTrialNo = self.wagerTrialNo + 1;
                                self.wagerRunLength = 1;
                            else
                                self.isWagerTrial = false;
                            end
                  
                        else
                            self.isWagerTrial = false;
                        end
                    end
                end
                
                
            end
            
            
            % generate a new outcome
            if isempty(self.fixedOutcomes)
                % pick a random outcome
                mu = self.currentMean;
                sigma = self.currentStd;
                
                %MOD FOR AUDITORY
                outcome = round(normrnd(mu, sigma));
                while outcome < self.minOutcome || outcome > self.maxOutcome
                    outcome = round(normrnd(mu, sigma));
                end
                
            else
                % pick the next fixed outcome
                self.fixedOutcomeIndex = self.fixedOutcomeIndex + 1;
                outcome = self.fixedOutcomes(self.fixedOutcomeIndex);
            end
            self.currentOutcome = outcome;
            
            %DEBUG
            if(self.isChangePointTask)
                fprintf('mean is %f\n',mu);
            else
                fprintf('outcome Index is %f\n',self.fixedOutcomeIndex);
            end
            fprintf('outcome is %f\n',self.currentOutcome);
            fprintf('trial number %f\n',self.blockTotalTrials);
            
            %Log data
            if(self.isChangePointTask)
                topsDataLog.logDataInGroup(self.currentMean,'mean');
                topsDataLog.logDataInGroup(self.currentStd,'std');
                topsDataLog.logDataInGroup(self.currentBlock,'block number');
                topsDataLog.logDataInGroup(self.currentHazard,'hazard');
            end
            % default to bad trial, unless set to good
            %self.isGoodTrial = false;
            
            
            self.wagerRunLength = self.wagerRunLength + 1;
        end
        
        
        
        % Assign a subject prediction for the current trial.
        function setPrediction(self, prediction)
            
            self.currentPrediction = prediction;
            %self.computeBehaviorParameters();
        end
        
        
        % Access the prediction which was reset or chosen by the subject.
        function prediction = getPrediction(self)
            prediction = self.currentPrediction;
        end
        
        %Access the current outcome
        function o = getCurrentOutcome(self)
            o = self.currentOutcome;
        end
        
        function setGoodTrial(self, isGoodTrial)
            self.isGoodTrial = isGoodTrial;
        end
        
        % Finish up the current trial.
        function finishTrial(self)
            
            %if self.isGoodTrial
            self.blockCompletedTrials = self.blockCompletedTrials + 1;
            %end
            self.blockTotalTrials = self.blockTotalTrials + 1;
            
            % if prediction is not active, revert to last prediction
            if ~self.isPredictionActive
                self.currentPrediction = self.lastPrediction;
            end
            
        end
        
        %finish up session and record data
        function finishSession(self)
            %do something useful if needed
        end
        
        % Update values for a change trial.
        function pickChangeTrial(self)
            self.currentMean = round(rand(1) .* self.maxOutcome);
            self.remainingSafety = self.safetyTrials;
            self.steadyTrials = 1;
        end
        
        
        
        % Summarize the current status of the session in a struct.
        function status = getStatus(self)
            props = properties(self);
            values = cell(size(props));
            for ii = 1:numel(props)
                values{ii} = self.(props{ii});
            end
            status = cell2struct(values, props);
        end
        
    end
end