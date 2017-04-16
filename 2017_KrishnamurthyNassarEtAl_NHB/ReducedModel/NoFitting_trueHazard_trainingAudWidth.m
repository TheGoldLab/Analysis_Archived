%% ************************************************************************************
%                      NO FITTING IN NORMATIVE MODEL WITH TRUE HAZARD AND
%                     AUDITORY NOISE WIDTH FROM TRAINING DATA
% *************************************************************************************



% Updated by MRN on 7/28/16 to:

% 1) verify that its accurately generating data and inferring latent
% variables

% 2) extend the generative process to include an additional "prediction
% noise" term... that is known to the model and accounted for in the prior
% width.




% ---- the hope is that the results using RU and CPP from this will hold up
% with the fit auditory width model (making the case stronger).

% ---- Kamesh March 2016

% ---- Adding additional subjects from new data collection (Kamesh MAY 16)

% load the file with the spatial audtiory noise widths
%Base directory
base_dir = '~/Dropbox/auditoryTaskManuscript/forArchiving/';
trainingFileDir = [base_dir 'Data/trainingData/'];  %directory which has training task data file

traininFN = 'train_likeWidth(17-May-2016).mat'; % single auditory likelhihood?
trainingFile = [trainingFileDir traininFN];
load(trainingFile);


% Where do you want to write the output file?
modelData_dir = [base_dir 'Data/modelData/'];
modelSimFileName = [modelData_dir 'modelSimulation_wPredNoise(' ...
    datestr(now,1) ').mat'];


numReps = 100;

doPredNoise=true; % if this is true, add extra noise to model predictions to account for subject sub-optimality


uniqueStds = [10,20];
uniqueIDs = unique(ID);
uniqueSess = unique(fixSessionID);
probeModelUncSubj = [];
probeModelStabSubj = [];
fitSubjParams = [];
probeModStabSim = [];
probeModUncSim = [];
probeModPredSim = [];
probeModEstSim = [];
subj_sess_norm_negLL =[];
subj_sess_norm_BIC = [];


% Set flags in input structure that will be fixed for all runs:
modelInput = struct;
modelInput.samplePosterior = false;
modelInput.samplePrior = false;
modelInput.use_internalUnc_for_sim = true; %for simulating data use noisy internal rep?
modelInput.subjectiveFlag = false; % Use subject prediction?
modelInput.whichParams=[false, false];  %fit the hazard rate





if doPredNoise
    % Load an old modeling data set (just need frugFun5 predictions):
    modelSimFileName='normative_simulated_behaviour(19-Jul-2016)';
    load(modelSimFileName)
    
    
    % compute prediction errors for subject and model:
    subPredErr=probeOutcomes'-probePrediction';
    subEstErr=probeOutcomes'- probeEstimate';
    modPredErr=probeOutcomes'-modPred(allIsProbe);
    for i0 = 1:length(unique(ID))
        sel_subj = probeSubjID == i0 &   probeDay1Sess' == 0 & probeTAC ~= 1;
        subj_vs_true_all(i0) =nanstd(subPredErr(sel_subj));
        mod_vs_true_all(i0) =nanstd(modPredErr(sel_subj));
    end
    
    % Look to see how much more variance there is in subject prediction errors:
    extraVarianceAll=subj_vs_true_all'.^2 -mod_vs_true_all'.^2;
end

% go back to figure directory:
%cd('~/Dropbox/auditorytask/AuditoryTaskV9/Analysis/Kamesh/Current Analyses/Figures')


for i1 = 1:length(uniqueSess)
    tic
    fprintf('simulating Subj -->  %d  session --> %d\n',ID(i1), i1);
    
    %selection variable for each probe (even non-fixated)
    selProbeSess = probeSessionID == uniqueSess(i1);
    % selection variables for each sound
    selAllBlk =  allSessID == uniqueSess(i1);
    subjID    = unique(probeSubjID(selProbeSess));
    
    
    % reset to data generation mode:
    modelInput.genDataFlag = false;
    
    
    % Put subject data in the modelInput structure:
    modelInput.trialNum = allTrialNum(selAllBlk);
    modelInput.probeTrial = allIsProbe(selAllBlk);
    modelInput.outcome = allOutcome(selAllBlk);
    modelInput.soundBlock = allSoundBlkNum(selAllBlk);
    
    %for each outcome location select the correct auditory noise width and
    %"extra prediction noise"
    modelInput.audWidth = ones(size(modelInput.outcome)).* train_PE_subject(subjID);
    
    
    % MRN: model crappy predictions by adding normally distributed noise
    % (predictionStd specifies width of this noise distriubtion)
    if doPredNoise==true
        modelInput.predictionStd= sqrt(extraVarianceAll(subjID));
    else
        modelInput.predictionStd= 0;
    end
    
    
    
    % select variables for each probe trial
    modelInput.blockNum = probeBlkNum(selProbeSess);
    modelInput.blockStd = probeStd(selProbeSess);
    modelInput.TAC = probeTAC(selProbeSess);
    
    % Feed the subject predictions and estimates
    predAux = probePrediction(selProbeSess)';
    predAux(predAux > 180) = 180;
    predAux(predAux < 0) = 0;
    modelInput.prediction = predAux;
    modelInput.CPP = nan(size(modelInput.blockNum));
    
    estAux = probeEstimate(selProbeSess)';
    estAux(estAux > 180) = 180;
    estAux(estAux < 0) = 0;
    modelInput.perception = estAux;
    
    % first entry is for hazard rate, second for motor std; auditory noise
    % width is NOT a parameter
    modelInput.startPoint = [.15, .1];
    modelInput.params = [.15, .1];    %use the true hazard (no fitting)!
    
    %[negLLik, modelOutput] = auditoryModel_gen_spatialAudWidth(modelInput); %fit_auditoryModel_spatialAudWidth(modelInput);
    [negLLik, modelOutput] = auditoryModel_gen_spatialAudWidth(modelInput);
    
    
    modelOutput.negLogLike = negLLik;
    fitSubjParams(i1) = modelOutput.params(1);
    
    % Now let's gather the stability term and the uncertainty term if the
    % model used the subject's reports
    audNoise = modelOutput.audWidth(modelInput.probeTrial)';
    totUncSubj = modelOutput.totUnc;         %total uncertainty estimated by model given subject reports etc
    modStableTermSubj = 1- modelOutput.CPP;
    uncTermSubj = audNoise.^2 ./ (totUncSubj.^2 + audNoise.^2);
    modePredSubj = modelOutput.prediction;  %this should be identical to subject PE  -- and it is
    modeEstSubj = modelOutput.perception; %this should be identical to subject est error  -- and it is
    
    probeModelStabSubj = [probeModelStabSubj; modStableTermSubj];
    probeModelUncSubj = [probeModelUncSubj; uncTermSubj];
    
    subj_sess_norm_negLL(i1) = modelOutput.negLogLike;
    num_sess_trials = sum(selProbeSess);
    subj_normative_BIC(i1) = 2*modelOutput.negLogLike + log(num_sess_trials) ; % 1 free parameter
    
    
    % Let us also simulate the prediction and perception reports by the
    % model instead of using the subject reports, and then collect the same
    % terms
    %clear modStabSim modUncSim modPredSim modEstSim
    modStabSim = [];
    modUncSim =[];
    modPredSim =[];
    modEstSim = [];
    
    modelInput.params = modelOutput.params;
    for i2 = 1:numReps
        
        %fprintf('Sim rep --> %d\n',i2);
        
        modelInput.genDataFlag = true; %we are simulating preds and est.
        modelInput.prediction = nan(size(modelInput.blockNum));
        modelInput.perception = nan(size(modelInput.blockNum));
        modelInput.CPP = nan(size(modelInput.blockNum));
        
        %clear modelOutput
        modelOutput = [];
        %[negLLik, modelOutput] = auditoryModel_gen_spatialAudWidth(modelInput);
        [negLLik, modelOutput] = auditoryModel_gen_spatialAudWidth_KK(modelInput);
        
        modUncSim(:,i2) = audNoise.^2 ./ ...
            (modelOutput.totUnc.^2 + audNoise.^2);
        modPredSim(:,i2) = modelOutput.prediction;
        modEstSim(:,i2) = modelOutput.perception; % this really is the only useful output here!
        modStabSim(:,i2) = 1 - modelOutput.CPP;
        
    end
    
    probeModStabSim = [probeModStabSim; modStabSim];
    probeModUncSim = [probeModUncSim; modUncSim];
    probeModPredSim = [probeModPredSim; modPredSim]; % simulated behavior!!!
    probeModEstSim = [probeModEstSim; modEstSim];
    
    toc
end

save(modelSimFileName,'probeModelUncSubj', 'probeModelStabSubj', 'fitSubjParams',...
    'probeModStabSim', 'probeModUncSim', 'probeModPredSim', 'probeModEstSim',...
    'subj_sess_norm_negLL', 'subj_sess_norm_BIC');

