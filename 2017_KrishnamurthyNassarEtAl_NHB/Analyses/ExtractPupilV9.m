%%
clear all
%close all
fileDir = '~/Box Sync/AuditoryTaskV9/RawData/';
writeDir = '~/Dropbox/AuditoryTask/AuditoryTaskV9/Data/ProcessedPupil/New Sampling/';

subjectList = {...
    'ZBD3S3fullTask(16-May-2016)'
    };

visualize = false;

%ENTER FALSE POSITIVE TRIAL NUMBERS HERE
fP = [];

if(length(subjectList)>1)
    error('Extract pupil one subject at a time');
end


%FILTER SPECIFICATION
d=fdesign.lowpass('Fp,Fst,Ap,Ast',4,6,1,40,60);
%myFilt = design(d);
myFilt1 = design(d,'butter','MatchExactly','passband');


extrFileName = [writeDir 'processed' subjectList{1} '.mat'];
numTAC = 6;

aux = load([fileDir subjectList{1} '.mat']);
data = aux.data;



idx = strcmp('trialETData',{data.group});
idx = find(idx);
pupilL = cell(length(idx),1);
pupilR = cell(length(idx),1);
pupTimes = cell(length(idx),1);
for i1 = 1:length(idx)
    auxPup = data(idx(i1)).item;
    pupilL{i1} = auxPup{1}(:,12);
    pupilR{i1} = auxPup{2}(:,12);
    pupTimes{i1} = auxPup{3};
end

idxFixed = strcmp('isFixedTrial',{data.group});
isFixed = cell2mat({data(idxFixed).item});

% get all the outcomes
idx1 = strcmp('outcome',{data.group});
outcome = data(idx1);
outcome = cell2mat({outcome.item});
% outcomeTimes = cell2mat({data(idx1).mnemonic});

idx1 = strcmp('tLocal',{data.group});
outcomeTimes = data(idx1);
outcomeTimes = cell2mat({outcomeTimes.item});
outcomeTimes = int64(outcomeTimes);

% get the indicies of wager outcomes (which of the indices in
% allOutcomes are wager trials?)
idx2 = strcmp('wagerTrialIdx',{data.group});
wagerTrialIdx = data(idx2);
wagerTrialIdx = cell2mat({wagerTrialIdx.item});


%DEBUG -- need to correct for block index since trial no. gets set to
%one at the start of each block
idxTP = strcmp('trialsPerBlock',{data.group});
numTrialsPerBlock = data(idxTP);
numTrialsPerBlock = cell2mat({numTrialsPerBlock.item});
%numTrialsPerBlock = 300;
blockStartIdx = find(diff([500 wagerTrialIdx 1])<0);
for i2 = 1:length(blockStartIdx)-1
    wagerTrialIdx(blockStartIdx(i2):blockStartIdx(i2+1)-1) = ...
        wagerTrialIdx(blockStartIdx(i2):blockStartIdx(i2+1)-1) + ...
        (i2-1)*numTrialsPerBlock;
end

%now record the outcomes corresponding to wager trials
wagerOutcomes = outcome(wagerTrialIdx);



idx3 = strcmp('percept',{data.group});
estimate = data(idx3);
estimate = cell2mat({estimate.item});
estimate = estimate*180/pi;

idx4 = strcmp('prediction',{data.group});
prediction = data(idx4);
prediction = cell2mat({prediction.item});
prediction = prediction*180/pi;

update = [diff(prediction) 0];

idx5=strcmp('wager',{data.group});
wager = data(idx5);
wager = cell2mat({wager.item});

idx6 = strcmp('mean',{data.group});
Mean = data(idx6);
Mean = cell2mat({Mean.item});

wagerMeans = Mean(wagerTrialIdx);
%outcomeTimes = outcomeTimes(wagerTrialIdx);

idx7 = strcmp('std',{data.group});
Std = data(idx7);
Std = cell2mat({Std.item});

idx9 = strcmp('isCpTrial', {data.group});
isCpTrial = data(idx9);
isCpTrial = cell2mat({isCpTrial.item});

%Block starts are considered as change-points because the TAC
%counter resets at the start of each block. This is not a big issue as
%there will only be 4 trials that may not be "true" change-points
step = numTrialsPerBlock;
isCpTrial([1 step+1 2*step+1 3*step+1]) = 1;


cpTrialNo = find(isCpTrial);

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


%Calculate the wagerTAC
wagerTAC = TAC(wagerTrialIdx)+1;

%convert the pupil time stamps into local time
pupTimesMat = cell(length(idx), 1);

% for j = 1:length(idx)
%     for k = 1:length(pupTimes{j})
%     pupTimesMat{j}(k) = tetio_remoteToLocalTime(int64(pupTimes{j}(k)));
%     end
% end

%align the sound time stamps with the pupil time stamps
alignTimes = cell(length(idx), 1);

% for j = 1:length(idx)
%     for k = 1:length(pupTimesMat{j})
%     alignTimes{j}(k) = tetio_remoteToLocalTime(outcomeTimes(j)) - pupTimesMat{j}(k);
%     end
% end



%************************ DO YOU WANT TO VISUALIZE THE DATA ************

%we visualize the interpolated trace and the raw trace on the same plot

numFixed = sum(isFixed==1);
idxFixed = find(isFixed);


pupilLinterp = cell(length(idxFixed),1);
pupilRinterp = cell(length(idxFixed),1);
filtPupL = cell(length(idxFixed),1);
filtPupR = cell(length(idxFixed),1);

for i2 = 1:numFixed
    aux1 = pupilL{idxFixed(i2)};
    aux2 = pupilR{idxFixed(i2)};
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
    
    pupilLinterp{i2} = aux1Interp;
    pupilRinterp{i2} = aux2Interp;
    
    
    %Filter the interpolated traces
    myFilt1.reset;
    prePad = ones(100,1)*aux1Interp(1);
    postPad = ones(100,1)*aux1Interp(end);
    auxFilt = [prePad; aux1Interp; postPad];     % DEBGU CHECK !!!!
    auxFilt = filtfilt(myFilt1.sosMatrix,myFilt1.ScaleValues,auxFilt);
    filtPupL{i2} = auxFilt(101:end-100);
    %filtPupL{i2} = filtfilt(myFilt1.sosMatrix,myFilt1.ScaleValues,aux1Interp);
    
    
    %Right Eye
    %myFilt = design(d);
    myFilt1.reset;
    prePad = ones(100,1)*aux2Interp(1);
    postPad = ones(100,1)*aux2Interp(end);
    auxFilt = [prePad; aux2Interp; postPad];     % DEBGU CHECK !!!!
    auxFilt = filtfilt(myFilt1.sosMatrix,myFilt1.ScaleValues,auxFilt);
    filtPupR{i2} = auxFilt(101:end-100);
    %filtPupR{i2} = filtfilt(myFilt1.sosMatrix,myFilt1.ScaleValues,aux2Interp);
end






if(visualize)
    % Pupil visualization
    idx = strcmp('trialETData',{data.group});
    idx = find(idx);
    aux1L = cell(length(idx),1);
    aux1R = cell(length(idx),1);
    pupTimes = cell(length(idx),1);
    for i1 = 1:length(idx)
        auxPup = data(idx(i1)).item;
        aux1L{i1} = auxPup{1}(:,12);
        aux1R{i1} = auxPup{2}(:,12);
        pupTimes{i1} = auxPup{3};
    end
    
    idxFixed = strcmp('isFixedTrial',{data.group});
    isFixed = cell2mat({data(idxFixed).item});
    idxFixed = find(isFixed);
    
    N = floor(length(isFixed)/16);
    
    for i1 = 1:N+1
        figure;
        for i2 = 1:16
            if((i1-1)*16+i2 <= length(isFixed))
                subplot(4,4,i2);
                plot(aux1L{(i1-1)*16+i2},'r')
                hold on
                plot(aux1R{(i1-1)*16+i2},'b')
                
                if(isFixed((i1-1)*16+i2))
                    auxFixIdx = find(idxFixed == (i1-1)*16+i2);
                    %                     plot(pupilLinterp{auxFixIdx},'r','linewidth',1.5)
                    %                     plot(pupilRinterp{auxFixIdx},'b','linewidth',1.5)
                    plot(filtPupL{auxFixIdx},'r','linewidth',1.5)
                    plot(filtPupR{auxFixIdx},'b','linewidth',1.5)
                end
                
                
                titleStr = ['Trial:' num2str((i1-1)*16+i2)...
                    '- Is fixed : ' ...
                    num2str(isFixed((i1-1)*16+i2)) ' TAC --> ',...
                    num2str(wagerTAC((i1-1)*16+i2)) ];
                title(titleStr);
            end
        end
    end
end




idx1 = strcmp({data.group},'trial states:enter:outcome');
idx2 = strcmp({data.group},'trial states:enter:predict');
idx3 = strcmp({data.group},'trial states:enter:fixate1');

idx4 = strcmp({data.group},'trial states:exit:fixate2');
idx5 = strcmp({data.group},'trial states:enter:wager');



outcomeTimes = {data(idx1).mnemonic};
outcomeTimes = cell2mat(outcomeTimes);

enterPredTimes = {data(idx2).mnemonic};
enterPredTimes = cell2mat(enterPredTimes);
enterPacedFixTimes = {data(idx3).mnemonic};
enterPacedFixTimes = cell2mat(enterPacedFixTimes);
timeToPredict = enterPacedFixTimes - enterPredTimes;

%make sure there's no reset in between blocks
if(any(diff(outcomeTimes)<0))
    error('Time log reset between blocks');
end

%make sure there's no reset in between blocks
if(any(timeToPredict<0))
    error('Time log reset between blocks');
end

% get the indicies of wager outcomes (which of the indices in
% allOutcomes are wager trials?)
idx2 = strcmp('wagerTrialIdx',{data.group});
wagerTrialIdx = data(idx2);
wagerTrialIdx = cell2mat({wagerTrialIdx.item});


%DEBUG -- need to correct for block index since trial no. gets set to
%one at the start of each block
idxTP = strcmp('trialsPerBlock',{data.group});
numTrialsPerBlock = data(idxTP);
numTrialsPerBlock = cell2mat({numTrialsPerBlock.item});
%numTrialsPerBlock = 300;
blockStartIdx = find(diff([500 wagerTrialIdx 1])<0);
for i2 = 1:length(blockStartIdx)-1
    wagerTrialIdx(blockStartIdx(i2):blockStartIdx(i2+1)-1) = ...
        wagerTrialIdx(blockStartIdx(i2):blockStartIdx(i2+1)-1) + ...
        (i2-1)*numTrialsPerBlock;
end

wagerOutcomeTimes = outcomeTimes(wagerTrialIdx);


%time spent in prediction stage
timeToPredict = timeToPredict(wagerTrialIdx);
if(any(timeToPredict < 0.05))
    error('something fishy with prediction duration');
end


% now calculate the time since the last probe. For now, just include 0
% for block boundaries
timeSinceLastProbe = diff([0 wagerOutcomeTimes]);
if(any(timeSinceLastProbe<0))
    error('negative timeSinceLastProbe');
end


% calculate the time spent in the estimation stage
mainFixationExitTimes = {data(idx4).mnemonic}';
mainFixationExitTimes = cell2mat(mainFixationExitTimes);

wagerEnterTimes = {data(idx5).mnemonic}';
wagerEnterTimes = cell2mat(wagerEnterTimes);
wagerEnterTimes = wagerEnterTimes(wagerTrialIdx);

timeInEstimate = wagerEnterTimes - mainFixationExitTimes;
if(any(timeInEstimate < 0.04))
    error('something fishy with estimation duration');
end



% correct the mislabeled trials
if(visualize)
    keyboard
end


if(~visualize)
    
    save(extrFileName,'isCpTrial', 'isFixed','outcome','pupilL', 'pupilR',...
        'numTrialsPerBlock', 'wager','wagerMeans','wagerOutcomes',...
        'wagerTAC','wagerTrialIdx','estimate','prediction',...
        'update','Mean','Std','TAC','cpTrialNo',...
        'outcomeTimes', 'pupTimes', 'pupTimesMat', 'alignTimes','fP',...
        'timeSinceLastProbe', 'wagerOutcomeTimes',...
        'outcomeTimes','timeToPredict','timeInEstimate');
    
    
    fprintf(['Wrote ' extrFileName '\n']);
    
    
end