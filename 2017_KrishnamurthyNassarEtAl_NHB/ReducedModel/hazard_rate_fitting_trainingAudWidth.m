%% ************************************************************************************
%                      FITTING HAZARD RATE IN NORMATIVE MODEL WITH
%                   SPATIALLY ADJUSTED  AUDITORY NOISE WIDTH FROM TRAINING DATA
% *************************************************************************************

% ---- the hope is that this will give better estimates of RU and CPP

% ---- Kamesh Feb. 2016

% load the file with the spatial audtiory noise widths
load train_likeWidth(22-Feb-2016)

modelSimFileName = ['~/Desktop/kameshModelSims_hazardFit(' ...
    datestr(now,1) ').mat'];

numReps = 100;

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


for i1 = 1:length(uniqueSess)
    
    fprintf('simulating Subj -->  %d  session --> %d\n',ID(i1), i1);
    
    %selection variable for each probe (even non-fixated)
    selProbeSess = probeSessionID == uniqueSess(i1);
    % selection variables for each sound
    selAllBlk =  allSessID == uniqueSess(i1);
    
    clear modelInput modelOutput
    
    %Do you want to use model priors or subject predictions?
    modelInput.samplePosterior = false;
    modelInput.subjectiveFlag = true;
    
    modelInput.genDataFlag = true;  % We first fit models to behaviour and then simulate data

    
    % which parameters are being fit? 
    % first entry is for hazard rate, second for motor std  auditory noise
    % width is not a parameter
    modelInput.whichParams=[true, false];  %fit the hazard rate
    
    modelInput.trialNum = allTrialNum(selAllBlk);
    modelInput.probeTrial = allIsProbe(selAllBlk);
    modelInput.outcome = allOutcome(selAllBlk);
    modelInput.soundBlock = allSoundBlkNum(selAllBlk);
    
    %for each outcome location select the correct auditory noise width
    subjID = unique(probeSubjID(selProbeSess));
    subj_audWidth = train_PE_subject_3_spatial_zones(subjID,:);
    auxOutcome = modelInput.outcome;
    aux_width = nan(length(auxOutcome),1);
%     spatial_idx1 = (auxOutcome >= 0 & auxOutcome <= 30) | ...
%                         (auxOutcome >= 150 & auxOutcome <= 180);
%                     
%     spatial_idx2 = (auxOutcome >= 30 & auxOutcome <= 60) | ...
%                         (auxOutcome >= 120 & auxOutcome <= 150);
%                     
%     spatial_idx3 = (auxOutcome >= 60 & auxOutcome <= 120);
%     
%     aux_width(spatial_idx1) = subj_audWidth(1);
%     aux_width(spatial_idx2) = subj_audWidth(2);
%     aux_width(spatial_idx3) = subj_audWidth(3);



%   
    [~, aux_spatial_idx] = min(abs(bsxfun(@minus,auxOutcome',mean_train_angles')),[],2);
    aux_width = train_PE_subject_spaceBinned_corrected(subjID,aux_spatial_idx);
    modelInput.audWidth = aux_width;
    
    
    % select variables for each probe trial
    modelInput.blockNum = probeBlkNum(selProbeSess);
    modelInput.blockStd = probeStd(selProbeSess);
    modelInput.TAC = probeTAC(selProbeSess);
    
    % Feed the subject predictions and estimates
    predAux = probePrediction(selProbeSess)';
    predAux(predAux > 180) = 180;
    predAux(predAux < 0) = 0;
    modelInput.prediction = predAux;
    
    
    %correct subject estimates based on their biases measured from training
    estAux = probeEstimate(selProbeSess)';
    aux_probe_outcome = probeOutcomes(selProbeSess)';
    [~, aux_spatial_idx] = min(abs(bsxfun(@minus,aux_probe_outcome,mean_train_angles')),[],2);
    estAux = estAux + train_PE_subject_bias_binned(subjID,aux_spatial_idx)';
    
    
    estAux(estAux > 180) = 180;
    estAux(estAux < 0) = 0;
    modelInput.perception = estAux;
        
    % first entry is for hazard rate, second for motor std; auditory noise
    % width is not a parameter
    modelInput.startPoint = [.15, .1];
    
    modelOutput = fit_auditoryModel_spatialAudWidth(modelInput);
    
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
    clear modStabSim modUncSim modPredSim modEstSim
    
    modelInput.params = modelOutput.params;
    for i2 = 1:numReps
        
        fprintf('Sim rep --> %d\n',i2);
        
        modelInput.genDataFlag = true; %we are simulating preds and est.
        modelInput.prediction = nan(size(modelInput.blockNum));
        modelInput.perception = nan(size(modelInput.blockNum));
        
        clear modelOutput
        [negLLik, modelOutput] = auditoryModel_gen_spatialAudWidth(modelInput);
        
        modUncSim(:,i2) = audNoise.^2 ./ ...
            (modelOutput.totUnc.^2 + audNoise.^2);
        modPredSim(:,i2) = modelOutput.prediction;
        modEstSim(:,i2) = modelOutput.perception;
        modStabSim(:,i2) = 1 - modelOutput.CPP;
        
    end
    
    probeModStabSim = [probeModStabSim; modStabSim];
    probeModUncSim = [probeModUncSim; modUncSim];
    probeModPredSim = [probeModPredSim; modPredSim];
    probeModEstSim = [probeModEstSim; modEstSim];
    
    
end

save(modelSimFileName,'probeModelUncSubj', 'probeModelStabSubj', 'fitSubjParams',...
    'probeModStabSim', 'probeModUncSim', 'probeModPredSim', 'probeModEstSim',...
    'subj_sess_norm_negLL', 'subj_sess_norm_BIC');

