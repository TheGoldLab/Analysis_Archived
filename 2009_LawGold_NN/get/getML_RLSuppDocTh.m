function [fits, fitsd] = getML_RLSuppDocTh(Monk, NSEN, NTUNE, ELN, DIRS, SIMNUM, BINS, SIMNAME, recompute, POW)
% get quickTs fitted parameters for reinforcement learning simulation data
% 
% the last input argument is optional, POW is the optional input for
% nonlinear pooling model


%%
if nargin<10
    POW = 1;
end

if ~strcmp(SIMNAME, 'NonLinPool')
    [h] = dirnames;
    fn  = ['getML_RLSuppDocTh_' SIMNAME '_' Monk '.mat'];
    savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
        filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
        filesep, 'Preparation', filesep, 'mat', filesep, fn];
else
    [h] = dirnames;
    fn  = ['getML_RLSuppDocTh_' SIMNAME '_' num2str(POW) '_' Monk '.mat'];
    savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
        filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
        filesep, 'Preparation', filesep, 'mat', filesep, fn];
end

if recompute
    fits   = nans(3,length(BINS),length(ELN));
    fitsd  = nans(3,length(BINS),length(ELN));

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
           
        % get thresholds
        
        for k = 1:length(BINS)
            % select trials
            LTR = x>BINS(k,1) & x<=BINS(k,2);

            % compute thresholds using time-dependent weibull function (don't
            % fix lapse)
            coh_ = coh(LTR);
            vt_  = vt(LTR);
            dir_ = dir(LTR);
            crt_ = crt(LTR);
            d = [coh_(:) vt_(:) sign(dir_(:)) crt_(:) ones(size(crt_(:)))];

            % compute lapse
            Llp = coh>=0.9 & vt>=1 & crt>=0;
            lapse  = sum(sum(crt(LTR&Llp)))/sum(sum(LTR&Llp));
            
            if isnan(lapse)
                lapse = 1;
            elseif lapse<0.8
                lapse = 0.8;
            end

            citype  = [];
            [fits(:,k,i), sems(:,k,i)] = ctPsych_fit(@quickTsFixLapse, d(:,1:2), d(:,4), [], citype, [], [1,0,0,0], 1-lapse);

        end
        
    end
        
            
    
%     fits  = nans(3,length(BINS), length(SIMNUM), length(ELN));
%     fitsd = nans(3,length(BINS), length(SIMNUM), length(ELN));
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
%                 [i j k]
% 
%                 % select trials
%                 LTR = [1:N]'>BINS(k,1) & [1:N]'<=BINS(k,2);
% 
%                 % compute lapse
%                 Llp = coh>=0.9 & vt>=1 & crt>=0;
%                 lapse  = sum(crt(LTR&Llp))/sum(LTR&Llp);
%                 if isnan(lapse)
%                     lapse = 1;
%                 elseif lapse<0.8
%                     lapse = 0.8;
%                 end
% 
%                 % compute thresholds using time-dependent weibull function
%                 d = [coh(LTR) vt(LTR) sign(dir(LTR)) ...
%                     crt(LTR) ones(size(crt(LTR)))];
%                      
% %                 citype         = [];
% %                 [f s]    = ctPsych_fit(@quickTsFixLapse, d(:,1:3), d(:,4), [], citype, [], [1,0,1,0]);
%                 citype  = [];
%                 [f s] = ctPsych_fit(@quickTsFixLapse, d(:,1:2), d(:,4), [], citype, [], [1,0,0,0], 1-lapse);
% 
%                 
%                 fits(:,k,i,j)  = real(f(:));
%                 fitsd(:,k,i,j) = real(s(:));
%             end
%         end
%     end

    save(savepath, 'fits', 'fitsd')

else
    load(savepath)

end

