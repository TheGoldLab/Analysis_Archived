function [lp, lpd, lN] = getML_RLFig2Lapse(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, BINS, SIMNAME, recompute, POW)
% get lapse rate for reinforcement learning simulation data

%%
if nargin<10
    POW = 1;
end

if ~strcmp(SIMNAME, 'NonLinPool')
    [h] = dirnames;
    fn  = ['getML_RLSuppDocLapse_' SIMNAME Monk '.mat'];
    savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
        filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
        filesep, 'Preparation', filesep, 'mat', filesep, fn];
else
    [h] = dirnames;
    fn  = ['getML_RLSuppDocLapse_' SIMNAME '_' num2str(POW) Monk '.mat'];
    savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
        filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
        filesep, 'Preparation', filesep, 'mat', filesep, fn];
end


if recompute
    warning off
    
    lp     = nans(length(BINS),length(ELN));
    lpd    = nans(length(BINS),length(ELN));
    lN     = nans(length(BINS),length(ELN));
    
    
    for i = 1:length(ELN)
        % load data
        coh = [];
        dir = [];
        vt  = [];
        crt = [];
        x   = [];
        for j = 1:length(SIMNUM)
            % get time constant for lapse
            [hdir, ldir, cdir, tdir] = dirnames;
            if ~strcmp(SIMNAME, 'NonLinPool')
                fname = ['/getML_RL' SIMNAME '_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(i)) '_v' num2str(SIMNUM(j)) '.mat'];
            else
                fname = ['/getML_RL' SIMNAME '_' Monk '_' num2str(POW) '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(i)) '_v' num2str(SIMNUM(j)) '.mat'];
            end
            load([tdir fname])
            fprintf([fname '\n']);

            N      = max(FIRA(:,1));
            coh    = [coh, FIRA(:,2)];
            dir    = [dir, FIRA(:,3)];
            vt     = [vt, FIRA(:,4)];
            crt    = [crt, FIRA(:,9)];
            x      = [x, [1:N]'];
        end
           
        % get lapse
        for k = 1:length(BINS)
            % select trials
            LTR = x>BINS(k,1) & x<=BINS(k,2) & coh>0.9 & vt>=1;

            coh_ = coh(LTR);
            vt_  = vt(LTR);
            crt_ = crt(LTR);

            % compute lapse
            lp(k,i)  = binofit(sum(crt_(:)),sum(LTR(:)));
            lN(k,i)  = sum(LTR(:));
            lpd(k,i) = lp(k)*(1-lp(k,i))/sqrt(lN(k,i));
        end

    end
        
    
     
    
%     lp  = nans(length(BINS), length(SIMNUM), length(ELN));
%     lpd = nans(2,length(BINS), length(SIMNUM), length(ELN));
%     lN  = nans(length(BINS), length(SIMNUM), length(ELN));
% 
%     for i = 1:length(SIMNUM)
%         for j = 1:length(ELN)
%             % load FIRA
%             % FIRA = [trial #, coh, dir, vt, pooled response, reward prob, beta, choice, reward]
%             [hdir, ldir, cdir, tdir] = dirnames;
%             fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(ELN(j)) '_v' num2str(SIMNUM(i)) '.mat'];
%             load([tdir fname])
%             fprintf([fname '\n']);
% 
%             % compute running lapse and thresholds
%             N      = max(FIRA(:,1));
%             coh    = FIRA(:,2);
%             dir    = FIRA(:,3);
%             vt     = FIRA(:,4);
%             crt    = FIRA(:,9);
% 
%             
%             for k = 1:length(BINS)
%                 % select trials
%                 LTR  = [1:N]'>BINS(k,1) & [1:N]'<=BINS(k,2) & coh>0.9 & vt>=1;
%                 coh_ = coh(LTR);
%                 vt_  = vt(LTR);
%                 crt_ = crt(LTR);
%                 [lp(k,i,j) lpd(:,k,i,j)] = binofit(sum(crt_(:)),sum(LTR(:)));
%                 lN(k,i,j) = sum(LTR(:));
%             end
%         end
%     end

    save(savepath, 'lp', 'lpd', 'lN')

    warning on
    
else
    load(savepath)

end

