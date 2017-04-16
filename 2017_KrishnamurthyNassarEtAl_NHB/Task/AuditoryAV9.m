classdef AuditoryAV9 < handle
    % MODIFIED FOR THE POST-DECISION WAGER TASK
    % --Kamesh 25/6/13
    
    properties
        
        blockStartTime;
        
        % hrir ID picked
        hrirID = 1005;
        
        %which type of task is being run?
        % Task types:
        % 1- 'fullTask' - full task (predict, estimate, feedback)
        % 2- 'estimateTask' - no prediction, just estimate, feedback (w/ CP structure)
        % 3- 'trainNoFeedback' - meas. percp. accuracy w/o CP structure. No feedback
        % 4- 'trainFeedback' - meas. percp. accuracy w/o CP structure. w/ feedback
        taskType;
        
        % a PredInfLogic object to work with
        logic;
        
        % a gray color that might be described as "light"
        lightGray = [1 1 1]*0.25;
        
        % a gray color that might be described as "medium"
        mediumGray = [1 1 1]*0.5;
        
        % a gray color that might be described as "dark"
        darkGray = [1 1 1]*0;
        
        arcColor =  [0 128 64]*0.5; % bluish
        
        estColor = [1 0 0.6];
        
        predColor = [1 0.5625 0];
        
        outcomeColor = [1 0 0];
        
        %Colors for isoluminant checkerboard
        %isoColor1 and isoColor2 are also used for the texture
        isoColor1 = [20 20 20]*1.0/255;
        isoColor2 = [30 30 30]*1.0/255; %40
        
        isoColor3 = [34 34 34]*1.0/255;
        isoColor4 = [30 30 30]*1.1/255;
        isoColor5 = [30 30 30]*1.1/255;
        isoColor6 = [13 13 13]*1.2/255;
        
        % degrees visual angle width of fixation point
        fixationSize = 1;
        
        % limits for drawing the arcedit AuditoryAV9
        arcLeftLimit = 200;
        arcRightLimit = -20;
        %inner and outer diameter for the arc
        arcRInner = 4; %2.5;
        arcROuter = 4.1;%2.57;
        arcRPrevCues1 = 5;%3.5;  %arc for displaying prev. pred., outcome, est.
        arcRPrevCues2 = 5.3;%3.8;
        
        
        %estimate made by subject
        estimate;
        
        %prediction made by subject
        prediction;
        
        
        %For post-decision wagering
        wArc;  %arc for wager window
        wText; %text at the bottom indicating 'HIGH' or 'LOW' wager
        wArcInner;  %2.3
        wArcOuter;  %2.45
        
        % fontsize for angles on ticks. This is different from the default
        % text 'fontsize' property. Two choices for different screen sizes
        tickFontSize1 = 20;
        tickFontSize2 = 32;
        
        fontSize1 = 36;
        fontSize2 = 64;
        fontSize3 = 30;
        
        % diameter for drawable targets
        estSides = 3;
        predSides = 4;
        prevPredSides = 4;
        dotSides = 30;
        dotHeight = 0.3;
        dotWidth = 0.3;
        prevDotHeight = 0.3;
        prevDotWidth = 0.3;
        
        %checker board settings
        checkerH = 10;
        checkerW = 10;
        
        % should the task run fullscreen
        isFullScreen = false;
        
        %fixation box dimensions
        fixX = 5;
        fixY = 5;
        
        %path to the sound file
        dataPath  = ['~/Data/WavFiles/subj' num2str(1005) '/'];
        
        
        %data input device
        ui;
        
        % how long to let the subject sit idle
        tIdle = 30;
        
        % how long to indicate that its time to predict
        tPredict = 0;
        
        % how long to let the subject update predictions
        tUpdate = 30;
        
        % how long to indicate subject commitment
        tCommit = 1;
        
        % how long to indicate the trial outcome
        tOutcome = 1;
        
        % how long to indicate the trial "delta" error
        tDelta = 1;
        
        % how long to indicate trial success
        tSuccess = 1;
        
        % how long to indicate trial failure
        tFailure = 1;
        
        % width of the task area in degrees visual angle
        width = 30;
        
        % height of the task area in degrees visual angle
        height = 20;
        
        % a string instructing the subject to wait, please
        pleaseWaitString = 'Please wait.';
        
        % a string instructing the subject to press a button
        pleasePressString = 'Press when ready.';
        
        % font size for all texts
        fontSize = 24;
        
        
        %EYE TRACKING RELATED
        %ASL eye tracker object
        eyeTracker;
        isEyeTracking;
        pupilMeasFcn;   %function that record pupil dia. when the sound is playing
        tLocal;
        
        %for self-paced fixation
        fixDuration;
        lastFixTime;
        pacedEtData;
        trialEtData;
        isFirstFixation;
        numBrokenFix;
        numWagerTrial;
        
        %tops classification object for fixation
        fixRegion;
        
        %auxiliary flag that keeps track of broken fixation during outcome
        %fixation
        isFixed;
        
        %tops groupedList object to access snowdots data structures
        list;
        
        %Is this a regular trial or a wager trial?
        isWagerTrial;
        
        %wager window size in degrees
        wagerWindowSize = 50;
        
        %low and high wager amoutns
        wagerLow = 5;
        wagerHigh = 15;
        
        %total score text
        scoreText;
        
        
        %low-latency sound player from Psychtoolbox
        pahandle;
        
        
        
        
        
    end
    
    
    properties(SetAccess = protected)
        
        % whether doing graphics remotely (true) or locally (false)
        isClient = false;
        
        %Have one ensemble which holds all drawables objects and
        %reference these objects by indexing into this ensemble.
        %Same idea for playables etc.
        
        % ensemble object for drawable objects
        drawables;
        
        %texture object for checkerboard texture
        checkTexture;
        
        %inter-block text string
        interBlockText;
        
        % fixation point (drawables index)
        fixation;
        
        % Arc (drawables index)
        arc;
        
        % Ticks (drawables index)
        ticks;
        
        % cursor indicating what the subject should do (drawables index)
        cursor;
        
        % estimator dot & line (drawables index) -- used to indicate posterior
        estimator;
        estLine;
        prevEstLine;
        prevEstText;
        
        % predictor dot & line (drawables index) -- used to indicate prior
        predictor;
        predLine;
        prevPredText;
        prevPredLine;
        prevPredText1;
        prevPredLine1;
        prevPredAngle;
        
        % center point (drawables index)
        center;
        
        %string which shows input angle
        angleText;
        
        
        % outcome dot (drawables index)  -- used to indicate actual outcome
        outcome;
        prevOutcomeLine;
        prevOutcomeText;
        
        
        % show last Nprev outcomes
        nPrev = 5;
        prevOutcomes = [];
        prevOutcomeTicks;
        rPrevOutTicks = 2.4; %Should be smaller than arcRInner
        
        
        
        %fixation window outline(drawables index)
        fixWindow;
        
        % ensemble object for playable objects
        playables;
        
        % sound to play with the trial outcome (playables index)
        outcomeSound;
        
        % display translation/centering related variables
        dispWidth;
        dispHeight;
        xShift;
        yShift;
        
        %are you predicting or estimating
        isPredicting;
        
        
        
        
        
    end
    
    methods
        
        %Constructor
        function self = AuditoryAV9(isClient)
            
            disp('--> constructor called');
            if nargin < 1 || isempty(isClient)
                isClient = false;
            end
            self.isClient = isClient;
            self.isEyeTracking = false;
            
            %Default text font size
            self.fontSize = 48;
            
            % ensemble for grouping drawables
            self.drawables = dotsEnsembleUtilities.makeEnsemble( ...
                'drawables', self.isClient);
            
            % ensemble for grouping playables
            self.playables = dotsEnsembleUtilities.makeEnsemble( ...
                'playables', self.isClient);
            
            % tell the ensemble how to open and close a window
            self.drawables.addCall({@dotsTheScreen.openWindow}, 'open');
            self.drawables.addCall({@dotsTheScreen.closeWindow}, 'close');
            % 'open' and 'close' won't run whenever the ensemble is run;
            % you need to invoke them explicitly
            self.drawables.setActiveByName(false, 'open');
            self.drawables.setActiveByName(false, 'close');
            
            
        end
        
        % Set up audio-visual resources as needed.
        function initialize(self)
            disp('Initializing->');
            
            %clear eye tracking data containers
            self.trialEtData = {};
            self.pacedEtData = {};
            self.isFirstFixation = true;
            
            
            %Wager window radius
            self.wArcInner = 0.92*self.arcRInner;  %2.3
            self.wArcOuter = 0.95*self.arcRInner;
            
            % create a 'playables' object for the outcome
            self.playables = dotsPlayableFile();
            self.playables.isBlocking = 0;
            
            self.prevPredAngle = pi/2;
            
            
            %----------------- Drawables objects ----------------------
            %texture object
            checkTexture1 = dotsDrawableTextures();
            checkTexture1.textureMakerFevalable = {@kameshTextureMaker,...
                self.checkerH,self.checkerW,...
                [],[],self.isoColor1,self.isoColor2};
            
            
            self.checkTexture = self.drawables.addObject(checkTexture1);
            
            %Text for inter-block screen
            interBlockText1 = dotsDrawableText();
            interBlockText1.color = self.isoColor3;
            interBlockText1.fontSize = self.fontSize;
            interBlockText1.string = 'Proceed When Ready!';
            self.interBlockText =  self.drawables.addObject(interBlockText1);
            
            %strings for prediction and estimation angles
            angleText1 = dotsDrawableText();
            angleText1.color = self.isoColor3;
            angleText1.fontSize = self.fontSize1;
            angleText1.isBold = true;
            angleText1.string = '';
            angleText1.x = 0; angleText1.y = -0.5;
            self.angleText = self.drawables.addObject(angleText1);
            
            
            
            % fixation target
            fix = dotsDrawableTargets();
            fix.colors = self.mediumGray;
            fix.width = self.fixationSize;
            fix.height = self.fixationSize;
            self.fixation = self.drawables.addObject(fix);
            
            % estimator dot
            est = dotsDrawableTargets();
            est.nSides = self.estSides;
            est.width = self.dotWidth;
            est.height = self.dotHeight;
            %est.colors = self.estColor;
            est.colors =  self.isoColor3; %self.isoColor2*1.1/255.0;
            self.estimator = self.drawables.addObject(est);
            
            %estimator line
            estLine1 = dotsDrawableLines();
            estLine1.isSmooth = true;
            estLine1.pixelSize = 2;
            estLine1.xFrom = 0; estLine1.yFrom = 0;
            estLine1.colors = est.colors;
            self.estLine = self.drawables.addObject(estLine1);
            
            
            %DEBUG -- show previous quantities
            %previous estimate text
            prevEstText1 = dotsDrawableText();
            
            %colorDEBUG
            prevEstText1.color = self.isoColor6;
            
            prevEstText1.fontSize = self.fontSize3;
            prevEstText1.string = 'E';
            prevEstText1.isBold = true;
            prevEstText1.isVisible = false;
            prevEstText1.x = 0; prevEstText1.y = 0;
            self.prevEstText = self.drawables.addObject(prevEstText1);
            
            %previous estimate guide line
            prevEstLine1 = dotsDrawableLines();
            prevEstLine1.pixelSize = 3;
            prevEstLine1.isSmooth = true;
            prevEstLine1.xFrom = 0; prevEstLine1.yFrom = 0;
            prevEstLine1.colors = est.colors;
            prevEstLine1.isVisible = false;
            self.prevEstLine = self.drawables.addObject(prevEstLine1);
            
            % predictor dot
            pred = dotsDrawableTargets();
            pred.height = self.dotHeight;
            pred.width = self.dotWidth;
            pred.nSides = self.predSides;
            %pred.colors = self.predColor;
            pred.colors = self.isoColor3; %self.isoColor2*1.1/255.0;
            self.predictor = self.drawables.addObject(pred);
            
            %DEBUG -- showing previous quantities
            %previous prediction text 1
            prevPredText1 = dotsDrawableText();
            
            %colorDEBUG
            prevPredText1.color = self.isoColor6;
            
            prevPredText1.fontSize = self.fontSize3;
            prevPredText1.string = 'P';
            prevPredText1.isBold = true;
            prevPredText1.x = 0; prevPredText1.y = 0;
            self.prevPredText = self.drawables.addObject(prevPredText1);
            
            %previous prediction guide line 1
            prevPredLine1 = dotsDrawableLines();
            prevPredLine1.pixelSize = 3;
            prevPredLine1.isSmooth = true;
            prevPredLine1.xFrom = 0; prevPredLine1.yFrom = 0;
            prevPredLine1.colors = self.isoColor4;
            self.prevPredLine = self.drawables.addObject(prevPredLine1);
            
            %previous prediction text 2
            prevPredText1 = dotsDrawableText();
            
            %colorDEBUG
            prevPredText1.color =  self.isoColor6;        %self.isoColor6;
            
            prevPredText1.fontSize = self.fontSize3;
            prevPredText1.string = 'P';
            prevPredText1.isBold = true;
            prevPredText1.x = 0; prevPredText1.y = 0;
            self.prevPredText1 = self.drawables.addObject(prevPredText1);
            
            %previous prediction guide line 2
            prevPredLine1 = dotsDrawableLines();
            prevPredLine1.pixelSize = 3;
            prevPredLine1.isSmooth = true;
            prevPredLine1.xFrom = 0; prevPredLine1.yFrom = 0;
            prevPredLine1.colors = self.isoColor6;
            self.prevPredLine1 = self.drawables.addObject(prevPredLine1);
            
            %predictor line
            predLine1 = dotsDrawableLines();
            predLine1.pixelSize = 2;
            predLine1.isSmooth = true;
            predLine1.xFrom = 0; predLine1.yFrom = 0;
            predLine1.colors = pred.colors;
            self.predLine = self.drawables.addObject(predLine1);
            
            % outcome dot
            outcom1 = dotsDrawableTargets();
            outcom1.nSides = self.dotSides;
            outcom1.height = self.dotHeight;
            outcom1.width = self.dotWidth;
            %outcom1.colors = self.outcomeColor;
            outcom1.colors = self.isoColor3; %self.isoColor2*1.1/255.0;
            self.outcome = self.drawables.addObject(outcom1);
            
            % previous outcome text
            prevOutcomeText1 = dotsDrawableText();
            
            %colorDEBUG
            prevOutcomeText1.color = self.isoColor6;      %self.isoColor6;
            
            prevOutcomeText1.fontSize = self.fontSize3;
            prevOutcomeText1.string = 'O';
            prevOutcomeText1.isBold = true;
            prevOutcomeText1.x = 0; prevOutcomeText1.y = -1.5;
            self.prevOutcomeText = self.drawables.addObject(prevOutcomeText1);
            
            % previous outcome guide line
            prevOutcomeLine1 = dotsDrawableLines();
            prevOutcomeLine1.pixelSize = 2;
            prevOutcomeLine1.isSmooth = true;
            prevOutcomeLine1.xFrom = 0; prevOutcomeLine1.yFrom = 0;
            prevOutcomeLine1.colors = self.isoColor6;
            self.prevOutcomeLine = self.drawables.addObject(prevOutcomeLine1);
            
            %fixation window outline
            fix1 = dotsDrawableLines();
            fix1.colors = self.isoColor2;
            
            % How big should the fixation box be?
            % convert from relative coordinates to visual degrees
            
            % *************** DBGTOB (need to convert correctly)
            %self.fixX = 5;
            %self.fixY = 5;
            
            fix1.x = [-self.fixX -self.fixX -self.fixX self.fixX self.fixX self.fixX...
                self.fixX -self.fixX];
            fix1.y = [-self.fixY self.fixY  self.fixY self.fixY self.fixY...
                -self.fixY -self.fixY -self.fixY];% + self.arcRInner/2; %vertical offset
            
            fix1.colors = self.isoColor2*1.1/255;
            fix1.isSmooth = true;
            fix1.isVisible = false;
            fix1.pixelSize = 3;
            self.fixWindow = self.drawables.addObject(fix1);
            
            
            
            %Post-decision wager stuff
            %wager window arc
            wArc1 = dotsDrawableArcs();
            wArc1.colors = self.isoColor3*1.05;
            wArc1.startAngle = 80;
            wArc1.sweepAngle = self.wagerWindowSize;
            wArc1.isSmooth = false;
            wArc1.nPieces = 10;
            wArc1.rInner = self.wArcInner;
            wArc1.rOuter = self.wArcOuter;
            self.wArc = self.drawables.addObject(wArc1);
            
            %wager text
            wText1 = dotsDrawableText();
            wText1.color = self.isoColor3;
            wText1.fontSize = self.fontSize1;
            wText1.isBold = true;
            wText1.string = '';
            wText1.x = 0; wText1.y = -0.5;
            self.wText = self.drawables.addObject(wText1);
            
            
            %Total score text
            sText1 = dotsDrawableText();
            sText1.color = self.isoColor3;
            sText1.fontSize = self.fontSize1;
            sText1.isBold = true;
            sText1.string = '';
            sText1.x = 0; sText1.y = -3.5; %-0.5;
            sText1.isVisible = true;
            self.scoreText = self.drawables.addObject(sText1);
            
            %             % DEBUG (Kamesh) Showing sequence of previous outcomes
            %             aux1 = dotsDrawableLines();
            %             aux1.pixelSize = 4;
            %             aux1.isSmooth = true;
            %             aux1.colors = self.isoColor4;
            %             xFrom = zeros(self.nPrev,1); xTo = xFrom;
            %             yFrom = zeros(self.nPrev,1); yTo = yFrom;
            %             self.prevOutcomeTicks;
            %             for i1 = 1:self.nPrev
            %                 self.prevOutcomes(end+1) = 90 + 5*randn;
            %
            %             end
            
            
            %Arc and ticks
            arc1 = dotsDrawableArcs();
            %arc1.colors = self.arcColor;
            arc1.colors = self.isoColor3; %self.isoColor2*1.1/255.0;
            arc1.startAngle = self.arcRightLimit;
            arc1.sweepAngle = self.arcLeftLimit - self.arcRightLimit;
            arc1.isSmooth = true;
            arc1.nPieces = 100;
            arc1.rInner = self.arcRInner;
            arc1.rOuter = self.arcROuter;
            self.arc = self.drawables.addObject(arc1);
            
            ticks1 = dotsDrawableLines();
            ticks1.pixelSize = 3;
            ticks1.isSmooth = true;
            ticks1.colors = arc1.colors;
            nTicks = 23;
            step = (self.arcLeftLimit - self.arcRightLimit)/(nTicks-1);
            rDot = (arc1.rInner+arc1.rOuter)/2;
            xTo = zeros(1,nTicks); xFrom = xTo;
            yFrom = zeros(1,nTicks); yTo = yFrom;
            for i1 = 1:nTicks
                xFrom(i1) = rDot*cos(((i1-1)*step+self.arcRightLimit)*pi/180);
                yFrom(i1) = rDot*sin(((i1-1)*step+self.arcRightLimit)*pi/180);
                if(mod(i1-1,2)==0)
                    xTo(i1) = (rDot+0.4)*cos(((i1-1)*step+self.arcRightLimit)*pi/180);
                    yTo(i1) = (rDot+0.4)*sin(((i1-1)*step+self.arcRightLimit)*pi/180);
                else
                    xTo(i1) = (rDot+0.2)*cos(((i1-1)*step+self.arcRightLimit)*pi/180);
                    yTo(i1) = (rDot+0.2)*sin(((i1-1)*step+self.arcRightLimit)*pi/180);
                end
            end
            ticks1.xTo = xTo; ticks1.xFrom = xFrom;
            ticks1.yTo = yTo; ticks1.yFrom = yFrom;
            self.ticks = self.drawables.addObject(ticks1);
            
            %Cursor
            curs = dotsDrawableTargets();
            curs.yCenter = 0.1*self.height;
            curs.colors = self.mediumGray;
            curs.width = self.fixationSize/2;
            curs.height = self.fixationSize/2;
            self.cursor = self.drawables.addObject(curs);
            
            
            
            % automate the task of drawing task graphics
            %[self.fixation, self.cursor, ...
            inds = [ self.checkTexture, self.predictor, self.outcome, ...
                self.arc, self.estimator, ...
                self.ticks, self.estLine, self.predLine, self.fixWindow, ...
                self.angleText];
            
            %For post-decision wager
            inds = [inds, self.wArc, self.wText, self.scoreText];
            
            %include previous estimate
            inds = [inds, self.prevEstLine, self.prevEstText];
            
            
            self.drawables.automateObjectMethod( ...
                'drawTask', @dotsDrawable.drawFrame, {}, inds, ...
                true);
            
            
            % open a drawing window and let objects react to it
            self.drawables.callByName('open');
            self.drawables.callObjectMethod(@prepareToDrawInWindow);
            
            %DEBGU -- trying to play faster sounds!! (60KHz sampling)
            % Configure the low-latency sound player from the Psychtoolbox
            InitializePsychSound;
            self.pahandle = PsychPortAudio('Open', [], [], 1, 44100, 2);
            
        end
        
        % Clean up audio-visual resources as needed.
        function terminate(self)
            self.drawables.callByName('close');
        end
        
        % Return concurrent objects, like ensembles, if any.
        function concurrents = getConcurrents(self)
            % TO CHECK -- where is this used? (Kamesh)
            concurrents = {self.drawables, self.playables};
        end
        
        
        
        
        % ------- state transition functions will be rewritten for the
        %-------- Auditory task (Kamesh)
        
        % Here's a description:
        % setupPrediction : function that sets up the screen for the prediction
        %                   task
        % setupOutcome : function that prepares the playable object and sets up the
        %                screen for outcome
        % doOutcome : play the sound file
        % setupEstimate : set up the screen for estimation
        % updateEstimate : function that records the subject's latest estimate
        % doSuccess : will show prediction, estimation and outcome
        % getInput will update the prediciton/estimation
        
        function setupPrediction(self)
            
            if(self.logic.isWagerTrial)
                %Put this outside the IF statement for bookkeeping
                self.isPredicting = true;
                
                %make the sound blocking
                playables = self.playables;
                playables.isBlocking = 0;
                
                %hide the text strings
                angleText = self.drawables.getObject(self.angleText);
                angleText.isVisible = false;
                
                % hide the estimator
                est = self.drawables.getObject(self.estimator);
                estLine = self.drawables.getObject(self.estLine);
                est.isVisible = true;
                estLine.isVisible = false;
                
                %hide outcome
                out = self.drawables.getObject(self.outcome);
                out.isVisible = false;
                
                %Show the outcome cues
                prevOutText = self.drawables.getObject(self.prevOutcomeText);
                prevOutLine = self.drawables.getObject(self.prevOutcomeLine);
                prevOutLine.isVisible = true;
                prevOutText.isVisible = true;
                
                
                % center the predictor
                pred = self.drawables.getObject(self.predictor);
                predLine = self.drawables.getObject(self.predLine);
                pred.isVisible = true;
                predLine.isVisible = true;
                pred.xCenter = 0;
                pred.yCenter = 0;
                predLine.xTo = pred.xCenter;
                predLine.yTo = pred.yCenter;
                
                %reveal total score text
                sText = self.drawables.getObject(self.scoreText);
                sText.isVisible = true;
                
                
                self.drawables.callByName('drawTask');
                self.drawables.getObject(self.checkTexture);
                
                topsDataLog.logDataInGroup(toc(self.blockStartTime),...
                    'tEnterPrediction');
            end
            
            
        end
        
        %for post-decision wager
        function setupWager(self)
            
            if(self.logic.isWagerTrial)
                
                % hide the predictor
                pred = self.drawables.getObject(self.predictor);
                predLine = self.drawables.getObject(self.predLine);
                pred.isVisible = false;
                predLine.isVisible = false;
                
                %hide outcome
                out = self.drawables.getObject(self.outcome);
                outcome = self.logic.getCurrentOutcome();
                rDot = (self.arcRInner + self.arcROuter)/2;
                out.xCenter = rDot*cos(outcome*pi/180);
                out.yCenter = rDot*sin(outcome*pi/180);
                out.isVisible = false;
                
                %hide the angle text strings
                angleText = self.drawables.getObject(self.angleText);
                angleText.isVisible = false;
                
                %show the wager arc
                wArc = self.drawables.getObject(self.wArc);
                wArc.isVisible = true;
                if(self.arcRightLimit > self.estimate - self.wagerWindowSize/2)
                    wArc.startAngle = self.arcRightLimit;
                    wArc.sweepAngle = self.estimate+self.wagerWindowSize/2 ...
                        - self.arcRightLimit;
                elseif(self.arcLeftLimit < self.estimate + self.wagerWindowSize/2)
                    wArc.startAngle = self.estimate - self.wagerWindowSize/2;
                    wArc.sweepAngle = self.arcLeftLimit - ...
                        (self.estimate - self.wagerWindowSize/2);
                else
                    wArc.startAngle = self.estimate - self.wagerWindowSize/2;
                    wArc.sweepAngle = self.wagerWindowSize;
                end
                
                
                %show the wager text strings
                wText = self.drawables.getObject(self.wText);
                wText.string = 'wager';
                wText.isVisible = true;
                
                % Hide the fixation window
                fixWindow = self.drawables.getObject(self.fixWindow);
                fixWindow.isVisible = false;
                
                %Show estimate guide line
                prevEstText = self.drawables.getObject(self.prevEstText);
                prevEstLine = self.drawables.getObject(self.prevEstLine);
                prevEstLine.isVisible = true;
                prevEstText.isVisible = true;
                
                
                %reveal total score text
                sText = self.drawables.getObject(self.scoreText);
                sText.isVisible = true;
                
                self.drawables.callByName('drawTask');
                
                topsDataLog.logDataInGroup(toc(self.blockStartTime),...
                    'tEnterWager');
                
            end
        end
        
        
        
        
        function setupOutcome(self)
            
            %keyboard
            
            % hide the predictor
            pred = self.drawables.getObject(self.predictor);
            predLine = self.drawables.getObject(self.predLine);
            pred.isVisible = false;
            predLine.isVisible = false;
            
            %hide outcome
            out = self.drawables.getObject(self.outcome);
            outcome = self.logic.getCurrentOutcome();
            rDot = (self.arcRInner + self.arcROuter)/2;
            out.xCenter = rDot*cos(outcome*pi/180);
            out.yCenter = rDot*sin(outcome*pi/180);
            out.isVisible = false;
            
            
            %hide the text strings
            angleText = self.drawables.getObject(self.angleText);
            angleText.isVisible = false;
            
            %post-decision wager stuff
            %hide the wager string
            wText = self.drawables.getObject(self.wText);
            wText.isVisible = false;
            %hide the wager arc
            wArc = self.drawables.getObject(self.wArc);
            wArc.isVisible = false;
            
            
            %hide the outcome cues
            rCue = self.arcRPrevCues2*1.1;
            prevOutText = self.drawables.getObject(self.prevOutcomeText);
            prevOutLine = self.drawables.getObject(self.prevOutcomeLine);
            prevOutText.x = rCue*cos(outcome*pi/180);
            prevOutText.y = rCue*sin(outcome*pi/180);
            prevOutText.rotation = outcome - 90;
            prevOutLine.xFrom = rDot*cos(outcome*pi/180);
            prevOutLine.xTo = 0.9*rCue*cos(outcome*pi/180);
            prevOutLine.yFrom = rDot*sin(outcome*pi/180);
            prevOutLine.yTo = 0.9*rCue*sin(outcome*pi/180);
            prevOutLine.isVisible = false;
            prevOutText.isVisible = false;
            
            %hide the estimator guide lines
            prevEstText = self.drawables.getObject(self.prevEstText);
            prevEstLine = self.drawables.getObject(self.prevEstLine);
            prevEstLine.isVisible = false;
            prevEstText.isVisible = false;
            
            %center the estimator
            est = self.drawables.getObject(self.estimator);
            estLine = self.drawables.getObject(self.estLine);
            est.isVisible = true;
            estLine.isVisible = true;
            est.xCenter = 0;
            est.yCenter = 0;
            estLine.xTo = est.xCenter;
            estLine.yTo = est.yCenter;
            
            
            %reveal total score text
            sText = self.drawables.getObject(self.scoreText);
            sText.isVisible = true;
            
            
            self.drawables.callByName('drawTask');
            
            % Now you are doone predicting
            self.isPredicting = false;
            
            
            
            % DEBUG -- try to play the sound file while entering the state
            % ********* this is to account for audioplayer latency
            %
            %playables = self.playables;
            
            %ISSUE THE PLAY COMMAND TO THE AUDIO PLAYER
            %playables.play();
            
            %USe the psychtoolbox low-latency player
            pahandle = self.pahandle;
            t1 = PsychPortAudio('Start', pahandle,  1, 0);
            
            
            
            
            %DEBUG
            %topsDataLog.logDataInGroup(topsClock(), 'debugEnterOutcome');
            %log the outcome at the right time(just before playing)
            topsDataLog.logDataInGroup(self.logic.getCurrentOutcome(),...
                'outcome');
            
            %EYE TRACKER --- get the local time when the play command was
            %issued to audioplayer
            
            topsDataLog.logDataInGroup(toc(self.blockStartTime),...
                'tOutcomePlay');
            
            %self.tLocal = tetio_localTimeNow;
            
            topsDataLog.logDataInGroup(tetio_localTimeNow,...
                'tLocalOutcome');
            
            
            %DEBUG
            %topsDataLog.logDataInGroup(topsClock(), 'debugLeaveOutcome');
        end
        
        
        function prepareSoundFile(self)
            
            %prepare the soundfile
            %             outcome = self.logic.getCurrentOutcome();
            %             playables = self.playables;
            %             auxOutcome = mod(outcome+270,360);
            %             fileName = [self.dataPath 'subj' num2str(self.hrirID) ...
            %                 '_az' num2str(auxOutcome) '_v1.wav'];
            %             %fileName = 'test.wav';
            %             playables.fileName = fileName;
            %             playables.prepareToPlay();
            
            
            %Use the Psychtoolbox low-latency audio player
            outcome = self.logic.getCurrentOutcome();
            pahandle = self.pahandle;
            auxOutcome = mod(outcome+270,360);
            fileName = [self.dataPath 'subj' num2str(self.hrirID) ...
                '_az' num2str(auxOutcome) '_v1.wav'];
            [s,Fs] = wavread(fileName);
            PsychPortAudio('FillBuffer', pahandle, s');
            
        end
        
        
        
        function doSuccess(self)
            
            %DEBUG (no visual cue during non-wager trial)
            if(self.logic.isWagerTrial)
                
                prevEstText = self.drawables.getObject(self.prevEstText);
                prevEstLine = self.drawables.getObject(self.prevEstLine);
                prevEstLine.isVisible = false;
                prevEstText.isVisible = false;
                
                if(self.logic.isWagerTrial)
                    %reveal the wager arc
                    wArc = self.drawables.getObject(self.wArc);
                    wArc.isVisible = true;
                    wArc.sweepAngle = self.wagerWindowSize;
                    
                    
                    %hide the wager text
                    wText = self.drawables.getObject(self.wText);
                    wText.isVisible = false;
                    
                    %reveal the estimate only on wager trials(from V7)
                    prevEstLine.isVisible = true;
                    prevEstText.isVisible = true;
                end
                
                
                %reveal the score always -- not just wager trials
                %reveal total score text
                sText = self.drawables.getObject(self.scoreText);
                sText.isVisible = true;
                
                
                
                %reveal outcome
                out = self.drawables.getObject(self.outcome);
                out.isVisible = true;
                
                %reveal outcome cues
                prevOutText = self.drawables.getObject(self.prevOutcomeText);
                prevOutLine = self.drawables.getObject(self.prevOutcomeLine);
                prevOutLine.isVisible = true;
                prevOutText.isVisible = true;
                
                
                %self.drawables.callByName('drawTask');
                s = dotsTheScreen.theObject;
                
                %texture
                tx = self.drawables.getObject(self.checkTexture);
                tx.draw;
                
                arc = self.drawables.getObject(self.arc);
                ticks = self.drawables.getObject(self.ticks);
                ticks.draw; arc.draw;
                if(self.taskType == 1)
                    %pred.draw; predLine.draw;
                end
                
                out.draw;
                
                if(self.logic.isWagerTrial)
                    wArc.draw; sText.draw;
                    prevEstLine.draw();
                    prevEstText.draw();
                end
                %angleText.draw; %outcomeText.draw;
                
                %DEBUG -- show previous quantities
                prevOutText.draw();
                prevOutLine.draw();
                s.nextFrame();
                
                
                %DEBUG (including the pause here instead of
                %configureAuditoryTask
                pause(0.3);
                
            else
                %no visual feedback for non-wager trials (unless training)
                %pause(0.1);
                %if(self.taskType == 4) %for trainFeedback
                
                %version with visual feedback on non-wager trials
                if(true)
                    
                    
                    %DEBUG -- give a short gap so that the display and the
                    %sound match
                    %pause(0.1);
                    
                    %reveal outcome
                    out = self.drawables.getObject(self.outcome);
                    out.isVisible = true;
                    
                    %reveal outcome cues
                    prevOutText = self.drawables.getObject(self.prevOutcomeText);
                    prevOutLine = self.drawables.getObject(self.prevOutcomeLine);
                    prevOutLine.isVisible = true;
                    prevOutText.isVisible = true;
                    
                    
                    %self.drawables.callByName('drawTask');
                    s = dotsTheScreen.theObject;
                    
                    %texture
                    tx = self.drawables.getObject(self.checkTexture);
                    tx.draw;
                    
                    % keep the estimator triangle and score text
                    % visible between trials
                    est = self.drawables.getObject(self.estimator);
                    sText = self.drawables.getObject(self.scoreText);
                    est.draw;
                    sText.draw;
                    
                    arc = self.drawables.getObject(self.arc);
                    ticks = self.drawables.getObject(self.ticks);
                    ticks.draw; arc.draw;
                    
                    out.draw;
                    
                    %DEBUG -- show previous quantities
                    prevOutText.draw();
                    prevOutLine.draw();
                    
                    s.nextFrame();
                    
                    %how long should the cue be up in non wager trials?
                    pause(0.22);
                    
                else
                    pause(0.05);
                end
            end
            
        end
        
        %Read input from the mouse
        function state = getMouseInput(self)
            
            % get input only on wager trials (for pred/est)
            if(self.taskType == 3 || self.taskType == 4 ...
                    || logical(self.logic.isWagerTrial) )
                %reveal the angle text string
                angleText = self.drawables.getObject(self.angleText);
                angleText.isVisible = true;
                
                %texture
                tx = self.drawables.getObject(self.checkTexture);
                
                m = self.ui;
                m.flushData;
                s = dotsTheScreen.theObject;
                arc = self.drawables.getObject(self.arc);
                rDot = (arc.rInner+arc.rOuter)/2;
                
                %DEBUG -- show previous quantities
                prevOutText = self.drawables.getObject(self.prevOutcomeText);
                prevOutLine = self.drawables.getObject(self.prevOutcomeLine);
                prevOutLine.isVisible = false;
                prevOutText.isVisible = false;
                prevEstText = self.drawables.getObject(self.prevEstText);
                prevEstLine = self.drawables.getObject(self.prevEstLine);
                
                sText = self.drawables.getObject(self.scoreText);
                
                if(self.isPredicting)
                    %reintroducing the prediction stage in V7
                    rCue = self.arcRPrevCues1;
                    dot = self.drawables.getObject(self.predictor);
                    l = self.drawables.getObject(self.predLine);
                    
                    %DEBUG -- show the cues
                    cueText = self.drawables.getObject(self.prevPredText);
                    cueText.isVisible = true;
                    cueLine = self.drawables.getObject(self.prevPredLine);
                    cueLine.isVisible = true;
                    
                    cueLine1.isVisible = true;
                    cueText1.isVisible = true;
                    
                    %DEBUG -- not sure why this is here
                    prevOutLine.isVisible = false;
                    prevOutText.isVisible = false;
                    
                else
                    %DEBUG -- show the guide lines
                    rCue = self.arcRPrevCues2;
                    cueText = self.drawables.getObject(self.prevEstText);
                    cueText.isVisible = true;
                    cueLine = self.drawables.getObject(self.prevEstLine);
                    cueLine.isVisible = true;
                    
                    %DEBUG -- not sure why this is here
                    prevOutLine.isVisible = false;
                    prevOutText.isVisible = false;
                    
                    dot = self.drawables.getObject(self.estimator);
                    l = self.drawables.getObject(self.estLine);
                end
                ticks = self.drawables.getObject(self.ticks);
                scaleFac = s.pixelsPerDegree;
                mXprev = m.x/scaleFac;
                mYprev = m.y/scaleFac;
                sensitivityFac =   0.6*0.9; %1.5*0.9; -- might want to lower this for motor error
                %theta = self.prevPredAngle; %start the pointer at the previous
                %predition ?
                theta = pi/2;
                state = [];
                while(isempty(state))
                    m.read();
                    state = m.getNextEvent();
                    mXcurr = m.x/scaleFac; mYcurr = -m.y/scaleFac;
                    dTheta = sensitivityFac*[mXcurr-mXprev mYcurr-mYprev]*...
                        [-dot.yCenter dot.xCenter]'/(rDot^2);
                    %[ theta dTheta mXcurr mYcurr m.isAvailable]
                    if(theta + dTheta >= (self.arcRightLimit*pi/180) && ...
                            theta + dTheta <= (self.arcLeftLimit*pi/180))
                        theta = theta + dTheta;
                    end
                    
                    %update the angleText
                    angleText.string =  sprintf('%3.0f',theta*180/pi);
                    dot.xCenter = rDot*cos(theta);
                    dot.yCenter = rDot*sin(theta);
                    l.xTo = dot.xCenter;
                    l.yTo = dot.yCenter;
                    
                    %DEBUG -- show previous quantities
                    %update the cue text and line
                    cueText.x = rCue*cos(theta);
                    cueText.y = rCue*sin(theta);
                    cueText.rotation = 1*theta*180/pi - 90;
                    cueLine.xFrom = rDot*cos(theta);
                    cueLine.xTo = 0.9*rCue*cos(theta);
                    cueLine.yFrom = rDot*sin(theta);
                    cueLine.yTo = 0.9*rCue*sin(theta);
                    
                    %All the drawing happens here
                    %texture
                    tx.draw;
                    
                    ticks.draw; arc.draw; dot.draw; l.draw; angleText.draw;
                    sText.draw;
                    cueLine.draw;
                    cueText.draw;
                    
                    s.nextFrame();
                    mXprev = mXcurr;
                    mYprev = mYcurr;
                end
                
                
                %log data
                if(self.isPredicting)
                    topsDataLog.logDataInGroup(theta,'prediction');
                    self.prediction = theta*180/pi;
                    
                else
                    topsDataLog.logDataInGroup(theta,'percept');
                    self.estimate = theta*180/pi;
                    % TO REMOVE (this is done in logETData) record good trial -- doing it here instead of 'doSuccess'
                    %topsDataLog.logDataInGroup(1,'isFixedTrial');
                end
                
            else
                state = '';
            end  % isWagerTrial
            
        end
        
        %Read input from keyboard
        function state  = getKbInput(self)
            
            %texture
            tx = self.drawables.getObject(self.checkTexture);
            
            kb = self.ui;
            kb.flushData;
            s = dotsTheScreen.theObject;
            arc = self.drawables.getObject(self.arc);
            rDot = (arc.rInner+arc.rOuter)/2;
            if(self.isPredicting)
                dot = self.drawables.getObject(self.predictor);
                l = self.drawables.getObject(self.predLine);
            else
                dot = self.drawables.getObject(self.estimator);
                l = self.drawables.getObject(self.estLine);
            end
            tick = self.drawables.getObject(self.ticks);
            scaleFac = s.pixelsPerDegree;
            thetaStep = 2e-2;
            kb.read();
            state = kb.getNextEvent();
            while(isempty(state))
                kb.read();state = kb.getNextEvent();
            end
            theta = pi/2;
            if(strcmp(state,'right'))
                theta = 0;
                dot.xCenter = rDot*cos(theta);
                dot.yCenter = rDot*sin(theta);
                l.xTo = dot.xCenter;
                l.yTo = dot.yCenter;
                
                %texture
                tx.draw;
                
                tick.draw;arc.draw; dot.draw; l.draw;
                s.nextFrame();
            elseif(strcmp(state,'left'))
                theta = pi;
                dot.xCenter = rDot*cos(theta);
                dot.yCenter = rDot*sin(theta);
                l.xTo = dot.xCenter;
                l.yTo = dot.yCenter;
                
                %texture
                tx.draw;
                
                tick.draw;arc.draw; dot.draw; l.draw;
                s.nextFrame();
            elseif(strcmp(state,'up'))
                theta = pi/2;
                dot.xCenter = rDot*cos(theta);
                dot.yCenter = rDot*sin(theta);
                l.xTo = dot.xCenter;
                l.yTo = dot.yCenter;
                
                %texture
                tx.draw;
                
                tick.draw;arc.draw; dot.draw; l.draw;
                s.nextFrame();
            end
            while(~strcmp(state,'pressed'))
                kb.read();
                state=kb.getNextEvent();
                if(strcmp(state,'right') || ...
                        dotsReadableHIDKeyboard.isEventHappening(kb,'right'))
                    if(theta - thetaStep >= 0 && theta - thetaStep <= pi)
                        theta = theta - thetaStep;
                        dot.xCenter = rDot*cos(theta);
                        dot.yCenter = rDot*sin(theta);
                        l.xTo = dot.xCenter;
                        l.yTo = dot.yCenter;
                        
                        %texture
                        tx.draw;
                        
                        tick.draw;arc.draw; dot.draw; l.draw;
                        s.nextFrame();
                    end
                elseif(strcmp(state,'left') || ...
                        dotsReadableHIDKeyboard.isEventHappening(kb,'left'))
                    if(theta + thetaStep >= 0 && theta + thetaStep <= pi)
                        theta = theta + thetaStep;
                        dot.xCenter = rDot*cos(theta);
                        dot.yCenter = rDot*sin(theta);
                        l.xTo = dot.xCenter;
                        l.yTo = dot.yCenter;
                        
                        %texture
                        tx.draw;
                        
                        tick.draw;arc.draw; dot.draw; l.draw;
                        s.nextFrame();
                    end
                end
            end
            
            %log data
            if(self.isPredicting)
                topsDataLog.logDataInGroup(theta,'prediction');
            else
                topsDataLog.logDataInGroup(theta,'percept');
                % TO REMOVE (this is done in logETData) record good trial -- doing it here instead of 'doSuccess'
                %topsDataLog.logDataInGroup(1,'isFixedTrial');
            end
        end
        
        % wait for the user to press some relevant key
        function state = waitForUser(self)
            ui = self.ui;
            ui.flushData;
            state = '';
            while(isempty(state))
                ui.read();
                state = ui.getNextEvent();
            end
        end
        
        
        %This is for post-decision wager -- get high or low wager
        %(TO BE CHANGED):
        % I want to implement this with mouse right & left clicks
        % but dotsReadableHIDMouse only detects one button in my Mac
        % trackpad. So I'm doing this with the keyboard for now, but we
        % should use the mouse in the psychophysics computer.
        
        function state = getWager(self)
            
            if(self.logic.isWagerTrial)
                
                
                %reveal the wager text string
                wText = self.drawables.getObject(self.wText);
                wText.isVisible = true;
                %wText.string = 'enter wager'; %Done earlier in setupWager
                
                %reveal the wager windo
                wArc = self.drawables.getObject(self.wArc);
                wArc.isVisible = true;
                
                %Estimate made by subject
                prevEstText = self.drawables.getObject(self.prevEstText);
                prevEstLine = self.drawables.getObject(self.prevEstLine);
                prevEstLine.isVisible = true;
                prevEstText.isVisible = true;
                
                %Update the score for this wager inside this function
                score = self.logic.blockScore;
                outcome = self.logic.getCurrentOutcome();
                wagerWindow = self.wagerWindowSize/2;
                winWager = (outcome - wagerWindow <= self.estimate) && ...
                    (self.estimate <= outcome + wagerWindow);
                topsDataLog.logDataInGroup(self.logic.blockScore,'score');
                
                %ticks
                ticks = self.drawables.getObject(self.ticks);
                
                %texture
                tx = self.drawables.getObject(self.checkTexture);
                
                
                kb = self.list{'input'}{'keyboard'};
                kb.flushData;
                s = dotsTheScreen.theObject;
                arc = self.drawables.getObject(self.arc);
                rDot = (arc.rInner+arc.rOuter)/2;
                
                state ='';
                
                while(~(strcmp(state,'left') || strcmp(state,'right')))
                    kb.read();
                    state = kb.getNextEvent();
                    
                    if(strcmp(state,'right') || ...
                            dotsReadableHIDKeyboard.isEventHappening(kb,'right'))
                        %Wager high
                        wText.string = 'HIGH';
                        %log the wager
                        topsDataLog.logDataInGroup(1,'wager');
                        %update score
                        if(winWager)
                            score = score + self.wagerHigh;
                        else
                            
                            %DEBUGscore
                            score = score - (self.wagerHigh*0.4 +2);
                            
                        end
                    elseif(strcmp(state,'left') || ...
                            dotsReadableHIDKeyboard.isEventHappening(kb,'left'))
                        %Wager Low
                        wText.string = 'LOW';
                        %log the wager
                        topsDataLog.logDataInGroup(0,'wager');
                        %update score
                        if(winWager)
                            score = score + self.wagerLow;
                        else
                            
                            %DEBUGscore
                            score = score - (self.wagerLow*0.4+1);
                            
                        end
                    end
                    
                    %update the score text and the total score
                    self.logic.blockScore = score;
                    sText = self.drawables.getObject(self.scoreText);
                    sText.string = num2str(score);
                    
                    %DEBUG WINDOW SIZE
                    wArc.sweepAngle = self.wagerWindowSize;
                    
                    tx.draw;
                    ticks.draw; arc.draw; wText.draw;
                    wArc.draw;
                    prevEstText.draw; prevEstLine.draw;
                    s.nextFrame();
                end
                
                tx.draw;
                ticks.draw; arc.draw; wText.draw;
                wArc.draw;
                prevEstText.draw; prevEstLine.draw;
                sText.draw();
                s.nextFrame();
                
            else
                state = '';
            end
        end
        
        
        % Give a message to the subject.
        function doMessage(self, message)
            disp('doMessage->');
            if nargin > 1
                disp(message);
            end
        end
        
        
        %Indicate failure trial
        function doFailure(self)
            disp('Failure trial!');
            
            %hide the fixation window;
            fixWindow = self.drawables.getObject(self.fixWindow);
            fixWindow.isVisible = false;
            self.drawables.callByName('drawTask');
            
            pause(0.7);
            aux = wavread('failure.wav');
            sound(aux);
            
            %make the subject repeat the failed outcome
            self.logic.fixedOutcomeIndex =  self.logic.fixedOutcomeIndex-1;
            self.logic.blockCompletedTrials = self.logic.blockCompletedTrials - 1;
            self.logic.blockTotalTrials = self.logic.blockTotalTrials - 1;
            
            %need to do this so that the 'block' tree node runs its child
            %'trial' the correct number of times
            tree = self.list{'outline'}{'tree'};
            n = tree.children{2}.children{1}.iterations;  %rhs is block.iterations
            tree.children{2}.children{1}.iterations = n+1;
            %fprintf('block iteration %d\n',tree.children{2}.children{1}.iterationCount);
            
            %TO REMOVE -- don't do this here. We record this in logETData
            %topsDataLog.logDataInGroup(0,'isFixedTrial');
            
        end
        
        
        %Draws the interblock screen and waits for user to proceed
        function waitForBlockStart(self)
            %texture
            tx = self.drawables.getObject(self.checkTexture);
            %text
            ibText = self.drawables.getObject(self.interBlockText);
            ibText.isVisible = true;
            s = dotsTheScreen.theObject;
            tx.draw; ibText.draw;
            s.nextFrame();
            ui = self.ui;
            ui.flushData;
            state = '';
            while(isempty(state))
                ui.read();
                state = ui.getNextEvent();
            end
            ibText.isVisible = false;
        end
        
        
        %*********************EYE TRACKING ****************************
        
        
        %set up display for fixation
        function setupPacedFixation(self)
            
            
            
            if(self.logic.isWagerTrial)
                
                
                %Store pupil from previous probe till now
                [leftEye, rightEye, timeStamp, trigSignal] = tetio_readGazeData();
                topsDataLog.logDataInGroup({leftEye, rightEye, timeStamp, trigSignal},...
                    'pupilSinceLastProbe');
                
                %DEBUG -- let us start tracking at the beginning of every block
                % and stop at the end. DO NOT DO IT HERE
                %tetio_startTracking;
                
                
                
                %this is for 'checkPacedFixation'
                self.fixDuration = 0;
                
                %set the isFixed flag to true; it will be flipped if fixation
                %is broken
                self.isFixed = true;
                
                %hide the text strings
                angleText = self.drawables.getObject(self.angleText);
                angleText.isVisible = false;
                outcomeText = self.drawables.getObject(self.prevOutcomeText);
                outcomeText.isVisible = false;
                
                
                % hide the predictor
                pred = self.drawables.getObject(self.predictor);
                predLine = self.drawables.getObject(self.predLine);
                pred.isVisible = false;
                predLine.isVisible = false;
                
                %hide outcome
                out = self.drawables.getObject(self.outcome);
                outcome = self.logic.getCurrentOutcome();
                rDot = (self.arcRInner + self.arcROuter)/2;
                out.xCenter = rDot*cos(outcome*pi/180);
                out.yCenter = rDot*sin(outcome*pi/180);
                out.isVisible = false;
                
                
                %center the estimator
                est = self.drawables.getObject(self.estimator);
                estLine = self.drawables.getObject(self.estLine);
                est.isVisible = true;
                estLine.isVisible = true;
                est.xCenter = 0;
                est.yCenter = 0; %self.arcRInner/2;
                estLine.xTo = est.xCenter;
                estLine.yTo = est.yCenter;
                
                %hide the wager arc
                wArc = self.drawables.getObject(self.wArc);
                wArc.isVisible = false;
                
                %hide the estimator guide lines
                prevEstText = self.drawables.getObject(self.prevEstText);
                prevEstLine = self.drawables.getObject(self.prevEstLine);
                prevEstLine.isVisible = false;
                prevEstText.isVisible = false;
                
                self.drawables.callByName('drawTask');
                
                %clear data carried over from prediction
                tetio_readGazeData;
                
                %Don't reveal window until end of paced fixation
                
                %                 %DEBUG -- give a short gap between predict and fixation
                %                 pause(0.8);
                %
                %                 %reveal the fixation window after the short pause;
                %                 fixWindow = self.drawables.getObject(self.fixWindow);
                %                 fixWindow.isVisible = true;
                %                 fixWindow.colors = self.isoColor2*1.1/255;
                %
                %                 self.drawables.callByName('drawTask');
                
                %record the times
                topsDataLog.logDataInGroup(toc(self.blockStartTime),...
                    'tEnterPacedFixation');
                
                
                topsDataLog.logDataInGroup(tetio_localTimeNow,...
                    'tLocalEnterPacedFixation');
                
                
                
            else
                fixWindow = self.drawables.getObject(self.fixWindow);
                fixWindow.isVisible = false;
            end
            
            if(self.taskType == 3 || self.taskType == 4)
                
                %DEBUG -- let us start tracking at the beginning of every block
                % and stop at the end. DO NOT DO IT HERE
                %tetio_startTracking;
                
                
                %this is for 'checkPacedFixation'
                self.fixDuration = 0;
                
                %set the isFixed flag to true; it will be flipped if fixation
                %is broken
                self.isFixed = true;
                
                %hide the text strings
                angleText = self.drawables.getObject(self.angleText);
                angleText.isVisible = false;
                outcomeText = self.drawables.getObject(self.prevOutcomeText);
                outcomeText.isVisible = false;
                
                
                % hide the predictor
                pred = self.drawables.getObject(self.predictor);
                predLine = self.drawables.getObject(self.predLine);
                pred.isVisible = false;
                predLine.isVisible = false;
                
                %hide outcome
                out = self.drawables.getObject(self.outcome);
                outcome = self.logic.getCurrentOutcome();
                rDot = (self.arcRInner + self.arcROuter)/2;
                out.xCenter = rDot*cos(outcome*pi/180);
                out.yCenter = rDot*sin(outcome*pi/180);
                out.isVisible = false;
                
                
                %center the estimator
                est = self.drawables.getObject(self.estimator);
                estLine = self.drawables.getObject(self.estLine);
                est.isVisible = true;
                estLine.isVisible = true;
                est.xCenter = 0;
                est.yCenter = 0; %self.arcRInner/2;
                estLine.xTo = est.xCenter;
                estLine.yTo = est.yCenter;
                
                %hide the wager arc
                wArc = self.drawables.getObject(self.wArc);
                wArc.isVisible = false;
                
                %hide the estimator guide lines
                prevEstText = self.drawables.getObject(self.prevEstText);
                prevEstLine = self.drawables.getObject(self.prevEstLine);
                prevEstLine.isVisible = false;
                prevEstText.isVisible = false;
                
                
                fixWindow = self.drawables.getObject(self.fixWindow);
                fixWindow.isVisible = false;
                
                self.drawables.callByName('drawTask');
                
                %For training tasks, clear data from before self-paced
                %fixation
                tetio_readGazeData;
            end
            % DEBUG -- how long does this take?
            %topsDataLog.logDataInGroup(topsClock(), 'debugLeaveFixation');
        end
        
        
        
        
        %set up display for fixation
        function setupFixation(self)
            
            if(self.logic.isWagerTrial)
                
                %clear data that might be carried over from paced fixation
                tetio_readGazeData;
                
                self.numWagerTrial = self.numWagerTrial+1;
                % DEBUG -- how long does this take?
                %topsDataLog.logDataInGroup(topsClock(), 'debugEnterFixation');
                
                % ********** DBGTOB ****************
                %if(self.isFirstFixation)
                %    tetio_startTracking;
                %    self.isFirstFixation = false;
                %end
                
                %this is for 'checkPacedFixation'
                self.fixDuration = 0;
                
                %set the isFixed flag to true; it will be flipped if fixation
                %is broken
                self.isFixed = true;
                
                %hide the text strings
                angleText = self.drawables.getObject(self.angleText);
                angleText.isVisible = false;
                outcomeText = self.drawables.getObject(self.prevOutcomeText);
                outcomeText.isVisible = false;
                
                
                % hide the predictor
                pred = self.drawables.getObject(self.predictor);
                predLine = self.drawables.getObject(self.predLine);
                pred.isVisible = false;
                predLine.isVisible = false;
                
                %hide outcome
                out = self.drawables.getObject(self.outcome);
                outcome = self.logic.getCurrentOutcome();
                rDot = (self.arcRInner + self.arcROuter)/2;
                out.xCenter = rDot*cos(outcome*pi/180);
                out.yCenter = rDot*sin(outcome*pi/180);
                out.isVisible = false;
                
                
                %center the estimator
                est = self.drawables.getObject(self.estimator);
                estLine = self.drawables.getObject(self.estLine);
                est.isVisible = true;
                estLine.isVisible = true;
                est.xCenter = 0;
                est.yCenter = 0; %self.arcRInner/2;
                estLine.xTo = est.xCenter;
                estLine.yTo = est.yCenter;
                
                %                 self.drawables.callByName('drawTask');
                %
                %                 %DEBUG -- give a short gap between predict and fixation
                %                 pause(0.8);
                
                %reveal the fixation window after the short pause;
                fixWindow = self.drawables.getObject(self.fixWindow);
                fixWindow.isVisible = true;
                fixWindow.colors = self.isoColor3;
                
                self.drawables.callByName('drawTask');
                
                topsDataLog.logDataInGroup(tetio_localTimeNow,...
                    'tLocalEnterMainFixation');
                topsDataLog.logDataInGroup(toc(self.blockStartTime),...
                    'tEnterMainFixation');
                
            else
                fixWindow = self.drawables.getObject(self.fixWindow);
                fixWindow.isVisible = false;
            end
            
            if(self.taskType == 3 || self.taskType == 4)
                
                %clear data that might be carried over from paced fixation
                tetio_readGazeData;
                
                %set the isFixed flag to true; it will be flipped if fixation
                %is broken
                self.isFixed = true;
                
                %hide the text strings
                angleText = self.drawables.getObject(self.angleText);
                angleText.isVisible = false;
                outcomeText = self.drawables.getObject(self.prevOutcomeText);
                outcomeText.isVisible = false;
                
                
                % hide the predictor
                pred = self.drawables.getObject(self.predictor);
                predLine = self.drawables.getObject(self.predLine);
                pred.isVisible = false;
                predLine.isVisible = false;
                
                %hide outcome
                out = self.drawables.getObject(self.outcome);
                outcome = self.logic.getCurrentOutcome();
                rDot = (self.arcRInner + self.arcROuter)/2;
                out.xCenter = rDot*cos(outcome*pi/180);
                out.yCenter = rDot*sin(outcome*pi/180);
                out.isVisible = false;
                
                
                %center the estimator
                est = self.drawables.getObject(self.estimator);
                estLine = self.drawables.getObject(self.estLine);
                est.isVisible = true;
                estLine.isVisible = true;
                est.xCenter = 0;
                est.yCenter = 0; %self.arcRInner/2;
                estLine.xTo = est.xCenter;
                estLine.yTo = est.yCenter;
                
                
                fixWindow = self.drawables.getObject(self.fixWindow);
                fixWindow.isVisible = true;
                
                fixWindow.colors = self.isoColor3;
                
                self.drawables.callByName('drawTask');
                
                
            end
            % DEBUG -- how long does this take?
            %topsDataLog.logDataInGroup(topsClock(), 'debugLeaveFixation');
        end
        
        %function that checks self-paced fixation
        function state =  checkPacedFixation(self)
            
            if(self.logic.isWagerTrial)
                
                auxTime = tic;
                
                auxPupilContainerR = []; %array that stores fix samples
                auxPupilContainerL = [];
                auxGazeContainerLx = [];
                auxGazeContainerLy = [];
                auxGazeContainerRx = [];
                auxGazeContainerRy = [];
                auxTimeStamp = [];
                auxTrig = [];
                isFixated = false;
                timeFixed = 0;
                
                while(~isFixated)
                    pause(1.0/60); % so that we have at least 1 sample
                    [leftEye, rightEye, timeStamp, trigSignal] = tetio_readGazeData();
                    
                    auxPupilContainerL = [auxPupilContainerL ;leftEye(:,12)];
                    auxPupilContainerR = [auxPupilContainerR ; rightEye(:,12)];
                    auxGazeContainerLx = [auxGazeContainerLx; leftEye(:,7)];
                    auxGazeContainerRx = [auxGazeContainerRx; rightEye(:,7)];
                    
                    auxGazeContainerLy = [auxGazeContainerLy; leftEye(:,8)];
                    auxGazeContainerRy = [auxGazeContainerRy; rightEye(:,8)];
                    
                    auxTimeStamp = [auxTimeStamp ; timeStamp];
                    auxTrig = [auxTrig; trigSignal];
                    
                    
                    
                    if(length(auxPupilContainerL)>36)
                        % CHECK FOR BROKEN FIXATION
                        maxMissLenL = 0;
                        maxMissLenR = 0;
                        maxOutLenL = 0;
                        maxOutLenR = 0;
                        
                        %coordinates of the eye gaze
                        aux1X = auxGazeContainerLx(end-36:end);
                        aux1Y = auxGazeContainerLy(end-36:end);
                        aux2X = auxGazeContainerRx(end-36:end);
                        aux2Y = auxGazeContainerRy(end-36:end);
                        
                        
                        %find the coordinates of the window - a1 b1, a2 b2
                        %indices for which the gaze is outside the window
                        outIdxL = find(aux1X - 0.6 > 0 | aux1X - 0.4 < 0 | aux1Y - 0.4 < 0 | aux1Y - 0.6 > 0);
                        outIdxR = find(aux2X - 0.6 > 0 | aux2X - 0.4 < 0 | aux2Y - 0.4 < 0 | aux2Y - 0.6 > 0);
                        
                        if(~isempty(outIdxL));
                            %starting points of the sequences of consecutive indices
                            outAuxL = [1 find(diff(outIdxL) > 1)' + 1];
                            outStartL = outIdxL(outAuxL);
                            
                            %ending points of the sequences of consecutive indices
                            outAuxEL = [1 find(diff(outIdxL) > 1)'];
                            outEndL = [outIdxL(outAuxEL(2:end)); outIdxL(end)];
                            
                            maxOutLenL = max(outEndL - outStartL);
                        end
                        
                        
                        if(~isempty(outIdxR))
                            %starting points of the sequences of consecutive indices
                            outAuxR = [1 find(diff(outIdxR) > 1)' + 1];
                            outStartR = outIdxR(outAuxR);
                            %ending points of the sequences of consecutive indices
                            outAuxER = [1 find(diff(outIdxR) > 1)'];
                            outEndR = [outIdxR(outAuxER(2:end)); outIdxR(end)];
                            
                            maxOutLenR = max(outEndR - outStartR);
                        end
                        
                        
                        %pupil diameter
                        pupil_diamL = auxPupilContainerL(end-36:end);
                        pupil_diamR = auxPupilContainerR(end-36:end);
                        
                        %find indices for which pupil diameter is -1
                        missIdxL = find(pupil_diamL < 0);
                        missIdxR = find(pupil_diamR < 0);
                        
                        if(~isempty(missIdxL) && ~isempty(missIdxR))
                            m1L = [1 find(diff(missIdxL) > 1)' + 1];
                            m1R = [1 find(diff(missIdxR) > 1)' + 1];
                            startIdxL = missIdxL(m1L);
                            startIdxR = missIdxR(m1R);
                            
                            m2L = [1 find(diff(missIdxL) > 1)'];
                            m2R = [1 find(diff(missIdxR) > 1)'];
                            endIdxL = [missIdxL(m2L(2:end)); missIdxL(end)];
                            endIdxR = [missIdxR(m2R(2:end)); missIdxR(end)];
                            
                            maxMissLenL = max(endIdxL - startIdxL);
                            maxMissLenR = max(endIdxR - startIdxR);
                            
                        end
                        
                        if (maxMissLenL > 9 || maxMissLenR > 9 || maxOutLenL > 9 || maxOutLenR > 9)
                            isFixated = 0;
                        else
                            isFixated = 1;
                        end
                        
                    end
                end
                
                
                if(isFixated)
                    state = 'outcome';
                    topsDataLog.logDataInGroup(tetio_localTimeNow,...
                        'tLocalLeavePacedFixation');
                    topsDataLog.logDataInGroup(toc(self.blockStartTime),...
                        'tLeavePacedFixation');
                else
                    state = 'else'
                    error('shouldnt be here');
                end
                
                self.pacedEtData = {auxPupilContainerL,auxPupilContainerR,...
                    auxGazeContainerLx, auxGazeContainerLy,...
                    auxGazeContainerRx, auxGazeContainerRy,...
                    auxTimeStamp, auxTrig};
                
                %DEBUG
                fprintf('Time spent in self-paced fixation -----------> %f \n',toc(auxTime));
                
            else
                state = 'outcome';
            end
            
            
            if (self.taskType == 3 || self.taskType == 4)
                pause(0.2);
                auxTime = topsClock;
                [leftEye, rightEye, timeStamp, trigSignal] = tetio_readGazeData();
                
                %DBGTOB
                % CHECK FOR BROKEN FIXATION
                maxMissLenL = 0;
                maxMissLenR = 0;
                maxOutLenL = 0;
                maxOutLenR = 0;
                
                %coordinates of the eye gaze
                aux1X = leftEye(:,7);
                aux1Y = leftEye(:,8);
                aux2X = rightEye(:,7);
                aux2Y = rightEye(:,8);
                
                
                %DEBUG -- in paced fixation check if they are looking
                %inside the box
                inIdxL = find(aux1X - 0.6 < 0 & aux1X - 0.4 > 0 & aux1Y - 0.4 > 0 & aux1Y - 0.6 < 0);
                inIdxR = find(aux2X - 0.6 < 0 & aux2X - 0.4 > 0 & aux2Y - 0.4 > 0 & aux2Y - 0.6 < 0);
                
                
                isFixated = 0;
                
                %keyboard
                
                if(length(inIdxR)>3 && length(inIdxL)>3)
                    isFixated = 1;
                end
                
                %make self-paced fixation time 0.6 seconds
                while(isFixated &&  self.fixDuration < 0.6 )
                    t = topsClock;
                    self.fixDuration = (t - auxTime);
                end
                
                % proceed to outcome if fixated for more than 500 ms(in while
                % loop above)
                if(self.fixDuration > 0.6)
                    state = 'outcome';
                else
                    state = 'else';
                    
                    %leaving this place blank for now
                    
                    %DEBUG -- don't log data from getValue
                    %self.pacedEtData = [];
                end
                
                self.pacedEtData = [self.pacedEtData,...
                    {leftEye, rightEye, timeStamp, trigSignal}];
                %get at least a few samples
                
            end
            
            
        end
        
        %function that checks fixation
        function notFixed =  getEyeData(self)
            
            %DEBUG -- let us start tracking at the beginning of every block
            % and stop at the end. DO NOT DO IT HERE
            %tetio_stopTracking();
            [leftEye, rightEye, timeStamp, trigSignal] = tetio_readGazeData();
            
            %DBGTOB
            % CHECK FOR BROKEN FIXATION
            maxMissLenL = 0;
            maxMissLenR = 0;
            maxOutLenL = 0;
            maxOutLenR = 0;
            
            %coordinates of the eye gaze
            aux1X = leftEye(:,7);
            aux1Y = leftEye(:,8);
            aux2X = rightEye(:,7);
            aux2Y = rightEye(:,8);
            
            %find the coordinates of the window - a1 b1, a2 b2
            %indices for which the gaze is outside the window
            outIdxL = find(aux1X - 0.8 > 0 | aux1X - 0.3 < 0 | aux1Y - 0.2 < 0 | aux1Y - 0.8 > 0);
            outIdxR = find(aux2X - 0.8 > 0 | aux2X - 0.3 < 0 | aux2Y - 0.2 < 0 | aux2Y - 0.8 > 0);
            
            if(~isempty(outIdxL));
                %starting points of the sequences of consecutive indices
                outAuxL = [1 find(diff(outIdxL) > 1)' + 1];
                outStartL = outIdxL(outAuxL);
                
                %ending points of the sequences of consecutive indices
                outAuxEL = [1 find(diff(outIdxL) > 1)'];
                outEndL = [outIdxL(outAuxEL(2:end)); outIdxL(end)];
                
                maxOutLenL = max(outEndL - outStartL);
            end
            
            
            if(~isempty(outIdxR))
                %starting points of the sequences of consecutive indices
                outAuxR = [1 find(diff(outIdxR) > 1)' + 1];
                outStartR = outIdxR(outAuxR);
                %ending points of the sequences of consecutive indices
                outAuxER = [1 find(diff(outIdxR) > 1)'];
                outEndR = [outIdxR(outAuxER(2:end)); outIdxR(end)];
                
                maxOutLenR = max(outEndR - outStartR);
            end
            
            %pupil diameter
            pupil_diamL = leftEye(:,12);
            pupil_diamR = rightEye(:,12);
            
            %find indices for which pupil diameter is -1
            missIdxL = find(pupil_diamL < 0);
            missIdxR = find(pupil_diamR < 0);
            
            if(~isempty(missIdxL) && ~isempty(missIdxR))
                m1L = [1 find(diff(missIdxL) > 1)' + 1];
                m1R = [1 find(diff(missIdxR) > 1)' + 1];
                startIdxL = missIdxL(m1L);
                startIdxR = missIdxR(m1R);
                
                m2L = [1 find(diff(missIdxL) > 1)'];
                m2R = [1 find(diff(missIdxR) > 1)'];
                endIdxL = [missIdxL(m2L(2:end)); missIdxL(end)];
                endIdxR = [missIdxR(m2R(2:end)); missIdxR(end)];
                
                maxMissLenL = max(endIdxL - startIdxL);
                maxMissLenR = max(endIdxR - startIdxR);
                
            end
            
            %label this as a broken fixation if you have >10 bad
            %samples
            brokenFixation = 0;
            if (maxMissLenL > 10 || maxMissLenR > 10 || maxOutLenL > 10 || maxOutLenR > 10)
                brokenFixation = 1;
            end
            
            self.isFixed = ~brokenFixation;
            
            %if fixation was broken, repeat this TAC
            if(brokenFixation)
                self.numBrokenFix = self.numBrokenFix +1;
                wagerTrialList = self.logic.wagerTrialList;
                self.logic.wagerTrialNo = self.logic.wagerTrialNo - 1;
                wagerTrialNo = self.logic.wagerTrialNo;
                self.logic.wagerTAC = wagerTrialList(wagerTrialNo);
                fprintf('------> Broken fixation -----> %d/%d \n',...
                    self.numBrokenFix, self.numWagerTrial);
            end
            
            %record the data for BOTH paced fixation and main fixation
            topsDataLog.logDataInGroup({leftEye, rightEye, timeStamp, trigSignal},'trialETData');
            topsDataLog.logDataInGroup(self.pacedEtData,'pacedETData');
            topsDataLog.logDataInGroup(self.isFixed,'isFixedTrial');
            
            %clear the gaze data containers
            self.pacedEtData = {};
            self.isFirstFixation = true;
        end
        
        
        %function that checks fixation
        function fixed =  checkTrainingFixation(self)
            
            %DBGTOB
            %DEBUG -- let us start tracking at the beginning of every block
            % and stop at the end. DO NOT DO IT HERE
            %tetio_stopTracking();
            [leftEye, rightEye, timeStamp, trigSignal] = tetio_readGazeData();
            
            %DBGTOB
            % CHECK FOR BROKEN FIXATION
            maxMissLenL = 0;
            maxMissLenR = 0;
            maxOutLenL = 0;
            maxOutLenR = 0;
            
            %coordinates of the eye gaze
            aux1X = leftEye(:,7);
            aux1Y = leftEye(:,8);
            aux2X = rightEye(:,7);
            aux2Y = rightEye(:,8);
            
            %find the coordinates of the window - a1 b1, a2 b2
            %indices for which the gaze is outside the window
            outIdxL = find(aux1X - 0.8 > 0 | aux1X - 0.3 < 0 | aux1Y - 0.2 < 0 | aux1Y - 0.8 > 0);
            outIdxR = find(aux2X - 0.8 > 0 | aux2X - 0.3 < 0 | aux2Y - 0.2 < 0 | aux2Y - 0.8 > 0);
            
            if(~isempty(outIdxL));
                %starting points of the sequences of consecutive indices
                outAuxL = [1 find(diff(outIdxL) > 1)' + 1];
                outStartL = outIdxL(outAuxL);
                
                %ending points of the sequences of consecutive indices
                outAuxEL = [1 find(diff(outIdxL) > 1)'];
                outEndL = [outIdxL(outAuxEL(2:end)); outIdxL(end)];
                
                maxOutLenL = max(outEndL - outStartL);
            end
            
            
            if(~isempty(outIdxR))
                %starting points of the sequences of consecutive indices
                outAuxR = [1 find(diff(outIdxR) > 1)' + 1];
                outStartR = outIdxR(outAuxR);
                %ending points of the sequences of consecutive indices
                outAuxER = [1 find(diff(outIdxR) > 1)'];
                outEndR = [outIdxR(outAuxER(2:end)); outIdxR(end)];
                
                maxOutLenR = max(outEndR - outStartR);
            end
            
            %pupil diameter
            pupil_diamL = leftEye(:,12);
            pupil_diamR = rightEye(:,12);
            
            %find indices for which pupil diameter is -1
            missIdxL = find(pupil_diamL < 0);
            missIdxR = find(pupil_diamR < 0);
            
            if(~isempty(missIdxL) && ~isempty(missIdxR))
                m1L = [1 find(diff(missIdxL) > 1)' + 1];
                m1R = [1 find(diff(missIdxR) > 1)' + 1];
                startIdxL = missIdxL(m1L);
                startIdxR = missIdxR(m1R);
                
                m2L = [1 find(diff(missIdxL) > 1)'];
                m2R = [1 find(diff(missIdxR) > 1)'];
                endIdxL = [missIdxL(m2L(2:end)); missIdxL(end)];
                endIdxR = [missIdxR(m2R(2:end)); missIdxR(end)];
                
                maxMissLenL = max(endIdxL - startIdxL);
                maxMissLenR = max(endIdxR - startIdxR);
                
            end
            
            %label this as a broken fixation if you have >10 bad
            %samples
            brokenFixation = 0;
            if (maxMissLenL > 10 || maxMissLenR > 10 || maxOutLenL > 10 || maxOutLenR > 10)
                brokenFixation = 1;
            end
            
            self.isFixed = ~brokenFixation;
            
            
            if(brokenFixation)
                self.isFixed = false;  %record the broken fixation. This will be logged in logETData
            end
            
            fixed = -1;   %DEBUG -- just a dummy output. make sure it doesn't break anything
            %record the data for BOTH paced fixation and main fixation
            topsDataLog.logDataInGroup({leftEye, rightEye, timeStamp, trigSignal},'trialETData');
            topsDataLog.logDataInGroup(self.isFixed,'isFixedTrial');
            
            
        end
        
        
        
        %log all the eye tracker history during second fixation in training
        % (outcome -> estimate). Also record broken fixations during the
        % trial
        function logTrainingETData(self)
            hist = self.eyeTracker.history;
            timedFrames = self.eyeTracker.timedFrames;
            topsDataLog.logDataInGroup({hist, timedFrames},'trainingETData');
            topsDataLog.logDataInGroup(self.isFixed,'isFixedTrial');
        end
        
        function state = isNextStateFixate(self)
            %this function checks which state is entered after 'outcome'
            %for wager trials fixate2 is entered otherwise 'estimate' is
            %enetered
            
            
            if(self.logic.isWagerTrial)
                state = 'fixate2';
            else
                state = 'estimate';
            end
            
        end
        
        
        function state = isNextStateFailure(self)
            %this function checks which state is entered after 'outcome'
            %for wager trials fixate2 is entered otherwise 'estimate' is
            %enetered
            if(self.isFixed)
                state = 'doNothing';
            else
                state = 'failure';
            end
            
        end
        
    end
end