function [fits, sems] = getML_neurometricROCcg(fname, bins, crtf, nvf, dirf, recompute)
% fit neural data with cumulative gaussian function for each good session
% INPUTS:
%   fname - 'txt' filename
%   bins  - time bins for which signal strength and fano factor will be
%           calculated ([begin time end time])
%   crtf  - [0 1]=both correct and incorrect, 1=correct only, 0=incorrect only
%   nvf   - [0 1]=both nv and not nv, 1=nv only, 0=not nv only
%   dirf  - 0=sorted by dot/trg direction, 1=sorted by choice

if ~isempty(findstr(fname, 'LIP'))
    dtype = 'LIP';

elseif ~isempty(findstr(fname, 'MT'))
    if ~isempty(findstr(fname, 'PRe'))
        dtype = 'MTPRe';
    else
        dtype = 'MT';
    end    
  
else
    return;
end

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_neurometricROCT_' dtype '_' fname(1:2) '.mat'];

if recompute
    if strcmp(dtype, 'MT') | strcmp(dtype, 'MTPRe') 
        d1name = 'ddir';
        ecname = 'dot_dir';
    elseif strcmp(dtype, 'LIP')
        d1name = 'trg_dir';
        ecname = 'trg1_dir';
    end

    % get average response before dots on
    a      = getML_txt(fname);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    usable = a.data{strcmp(a.name,'usable')};
    uid    = a.data{strcmp(a.name,'uid')};
    d1     = a.data{strcmp(a.name,d1name)};
    tkflg  = a.data{strcmp(a.name,'taskflag')};
    
    coh    = [0 3.2 6.4 12.8 25.6 51.2 99.9]';

    spec   = [1 0 0 0]; 
    fits   = nans(sum(spec),length(usable));
    sems   = nans(sum(spec),length(usable));
   
    global FIRA
    for i = 1:length(fn)
        if usable(i) == 1            
            fprintf('%d: %s\n', i, fn{i})
            openFIRA(fn{i})

            trials      = find(~isnan(FIRA.ecodes.data(:,1)));
            [Ld, Ud]    = selectFIRA_trialsByUniqueID(ecname, trials);
            Ld          = Ld(:,[find(round(Ud)==d1(i)) find(round(Ud)==mod(d1(i)+180,360))]);
            [Lc, Uc]    = selectFIRA_trialsByUniqueID('dot_coh', trials);
            spikei      = getFIRA_spikeByID(uid(i));
            [Lct, Uct]  = selectFIRA_trialsByUniqueID('correct', trials);
            if ismember(0,Uct) & ismember(1,Uct)
                Lct = Lct(:,Uct==0) | Lct(:,Uct==1);
            elseif ismember(0,Uct)
                Lct = Lct(:,Uct==0);
            elseif ismember(1,Uct)
                Lct = Lct(:,Uct==1);
            else
                Lct = logical(zeros(length(Lct),1));
            end
            [Ltk, Utk]  = selectFIRA_trialsByUniqueID('task', trials);
            if tkflg(i)<=2
                Ltk     = Ltk(:,Utk==3);
            else
                Ltk     = Ltk(:,Utk==6);
            end
           
            for j = 1:size(Ld,2)
                Ld(:,j) = Ld(:,j) & Lct & Ltk;       % select for correct/incorrect and task
            end
            bt          = getFIRA_ecodeTimesByName('dot_on',  0, trials);
            et          = getFIRA_ecodeTimesByName('dot_off', 0, trials);
            minn        = 3;
            
             % get grid of ROCs
            [rocs_, ns_, sems_] = getFIRA_neurometricROCT(trials, Ld, Lc(:,ismember(Uc,coh)),...
                spikei, bt, et, bins, minn);


            % get fits
            fit_dat   = grid2mat(rocs_, ns_, bins(:,2)-50, Uc(ismember(Uc,coh)));
            fit_dat   = [fit_dat(:, 1)./100 fit_dat(:, 2)./1000 ones(size(fit_dat, 1), 1) fit_dat(:, [3 4])];
            time_bins = nonanunique(fit_dat(:,2));

            
            %[fits,sems,stats] = quickTs_fit(fit_dat,2,1,[],[],time_bins);
            [fits(:,i),sems(:,i)] = ctPsych_fit(@cGauss, fit_dat(:,1:3), fit_dat(:,4:5), [], spec);
           fits(:,i);
            % plot threshold and slope with time
            if 0
                lc  = [];
                for j = 1:7
                    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
                end
                
                % plot raw data and fitted data
                clf
                hold on
                t  = (bins(:,2)-50)/1000;
                c  = [0 3.2 6.4 12.8 25.6 51.2 99.9];
                
                
                for j = 1:length(Uc)
                    % no lapse, fix sd
                    m     = fits(1,i)*c(j).*t/100;
                    sd    = (10+m).^0.5;
                    val = 0.5+(0.5-normcdf(0,m,sd));

%                     % no lapse, 1 variable sd
%                     m     = fits(1,i)*c(j).*t/100;
%                     sd    = (15+fits(2,i).*m).^0.5;
%                     val = 0.5+(0.5-normcdf(0,m,sd));


%                     %sd    = (fits(2,i)+fits(3,i).*m).^0.5;
%                     
%                     val = 0.5+(0.5-fits(2,i)).*(1-2*normcdf(0,m,sd));
%                     
                    
                    Ic = find(ismember(c ,Uc(j)));
                    plot(t, rocs_(:,j), 'o', 'MarkerFaceColor', lc(Ic,:), 'MarkerEdgeColor', 'none');
                    plot(t, val, 'LineWidth', 2, 'Color', lc(Ic,:));
                end
                hold off
                ylim([0.4 1]);
                % plot fit

                pause
            end

            
        end
    end

    fits = real(fits);
    sems = real(sems);
    save(savepath, 'fits', 'sems')
else

    load(savepath)
end
