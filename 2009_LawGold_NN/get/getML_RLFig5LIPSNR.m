function [SNR, trial] = getML_RLFig5LIPSNR(fname, crtf, nvf, dirf, bins, recompute)

[h] = dirnames;
fn  = ['getML_RLFig5LIPSNR_' fname(1:2) '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];



if recompute
    
    %**** NEW WAY: ESTIMATE SNR INDIRECTLY BY FIRST COMPUTING ROC AREA OF CHOICE ****%
    %****          THEN CONVERT ROC AREA BACK TO D', ASSUMING NOISE IS NORMAL    ****%
    [rocs, rocn, roce] = getML_ROCs(fname, bins, crtf, nvf, dirf, [], 1);
    SNR = norminv(rocs);
    SNR(rocn<=16) = nan;
    
   
    

    %**** OLD WAY: ESTIMATE SNR DIRECTLY FROM MEAN AND SD OF LIP RESPONSES ****%
    %
    % % get SNR
    % warning off
    % [rM, rSD, rN] = getML_mRate(fname, bins, crtf, nvf, dirf, recompute);
    % warning on
    %
    % to get a good estimate of SD, only include cell with 20 or more
    % trials in that condition
    %    
    % rM(rN<=20)  = nan;
    % rSD(rN<=20) = nan;
    %
    %
    % 
    % rM  = rM(1:7,:,:)-rM(8:14,:,:);
    % rSD = sqrt(rSD(1:7,:,:).^2+rSD(8:14,:,:).^2);
    % SNR = rM./rSD;
    
    
    % get trials
    a   = getML_txt(fname);
    ses = a.data{strcmp(a.name,'session')};
    
    load([fname(1:2) 'Combined.mat'])
    crt         = data(:,3);
    Lgd         = crt>=0 & ~isnan(data(:,6));
    [ses_, tr_] = unique(data(Lgd,1));
    
    for i = 1:length(ses)
        trial(i) = tr_(ses_==ses(i));
    end
    
    save(savepath, 'SNR', 'trial')
    
else
    load(savepath)

end
