function figLCP_pdBeep(num, collectData)
% function figLCP_pdBeep(num, collectData)
%
%  siteData{1}: trialsxcols matrix, cols are:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correcdt
%   4 ... beep on time, wrt to fix start time
%   5 ... trial begin time, wrt to fix start time
%   6 ... trial end time, wrt to fix start time
%   7 ... trial wrt time (cpu clock)
%   8 ... LFP index corresponding to fix start time (coded above)
%
%  siteData{2}: Analog:
%   dim1: trial
%   dim2: sample
%   dim3: 1 = x, 2 = y, 3 = pupil
%   [remember first sample is considered time=0 (i.e., wrt fix start time)]
%     = eyedat(Lgood,:,:);
%
%  siteData{3}: spikes, re-coded wrt fix start time
%       first is multi-unit
%       rest are single-unit
%
%  siteData{4}: LFP, re-coded wrt fix start time
%
% SiteData{5}: pupil events... rows are events, columns are:
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. end time of event (wrt fix start time)
%   4. magnitude at start of event (raw z-score)
%   5. magnitude at end of event (raw z-score)
%   6. magnitude at start of event (corrected z-score)
%   7. magnitude at end of event (corrected z-score)
%   8. time of subsequent max slope
%   9. magnitude of subsequent max slope (corrected z/sample)

%% parse args
if nargin < 1 || isempty(num)
    num = 7;
end

if nargin < 2 || isempty(collectData)
    collectData = false;
end

%% Set up figure
wid  = 17.6; % total width
ht   = 3;
cols = {5,5,1};
axs  = getPLOT_axes(num, wid, ht, cols, 1.5, .5, [], 'Joshi et al', true);
set(axs,'Units','normalized');

%% Stuff we need
% site, monkeys, example_data (monkey file ylims)
sites = {...
    'LC',  {'Oz' 'Cicero'},      [1 13  1 10]; ...
    'IC',  {'Oz' 'Cicero'},      [1 11 10 30]; ...
    'SC',  {'Oz' 'Cicero'},      [1  3  6 15]; ...
    'ACC', {'Sprout' 'Atticus'}, [2  3  3  7]; ...
    'PCC', {'Sprout' 'Cheetah'}, [1 13 13 20]};
num_sites = size(sites, 1);

% time axis for PETH
tmin   = -1000;
tmax   = 1500;
tsize  = 200;
tstep  = 10;
tbins  = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
ntbins = size(tbins,1);
tb     = find(tbins(:,1)==0);
tax    = mean(tbins,2);
pax    = (tmin:tmax)';
npax   = length(pax);

%% Possibly get data
if collectData
    
    udat = cell(num_sites, 4); % spike/sites, pupil/sites
    
    for ss = 1:num_sites
        for mm = 1:size(sites{ss,2},2)
            % load file
            [base_dir, fnames] = getLCP_cleanDataDir(sites{ss,2}{mm}, sites{ss,1});
            if ~isempty(fnames)
                for ff = 10:length(fnames)
                    disp(sprintf('monkey = %s, site = %s, %d/%d files', ...
                        sites{ss,2}{mm}, sites{ss,1}, ff, length(fnames)))
                    load(fullfile(base_dir, fnames{ff}));
                    
                    % beep trials, single-units only
                    Fbeep      = find(isfinite(siteData{1}(:,4)));
                    num_sunits = size(siteData{3}, 2)-1;
                    
                    if ~isempty(Fbeep) && num_sunits > 0
                        ntr   = size(Fbeep,1);
                        lpd   = size(siteData{2},2);
                        sdat  = nans(ntr, ntbins, num_sunits);
                        sbdat = nans(ntr, num_sunits);
                        pdat  = nans(ntr, npax, 2);
                        pbdat = nans(ntr, 1);
                     
                        for bb = 1:length(Fbeep)
                            
                            % spikes wrt beeps
                            btm  = siteData{1}(Fbeep(bb),4);
                            btax = btm+tbins;
                            si   = find(btax(:,1)>=siteData{1}(Fbeep(bb),5), 1);
                            ei   = find(btax(:,2)<=siteData{1}(Fbeep(bb),6), 1, 'last');
                            if ei > si
                                for uu = 1:num_sunits
                                    sp = siteData{3}{Fbeep(bb), uu+1};
                                    for tt = si:ei
                                        sdat(bb,tt,uu) = sum( ...
                                            sp>=btax(tt,1)&sp<=btax(tt,2));
                                    end
                                    % baseline is pre-beep
                                    if si<tb
                                        sbdat(bb,uu) = sum(sp>=btax(si,1)&sp<=btm)...
                                        ./(btm-btax(si,1)).*1000;
                                    end
                                end                                
                            end                            
                            
                            % PD wrt beeps
                            bpax = pax+round(btm);
                            Lp   = bpax>=0&bpax<lpd;
                            if any(Lp)
                                % save PD, slope
                                pdat(bb,Lp,1) = siteData{2}(Fbeep(bb),bpax(Lp)+1,4);
                                pdat(bb,Lp,2) = siteData{2}(Fbeep(bb),bpax(Lp)+1,5);
                                
                                % PD baseline:
                                pbdat(bb) = nanmean(pdat(bb,pax>=0&pax<100,1));
                                % pbdat(bb) = nanmean(pdat(bb,pax>=-1000&pax<0,1));
                            end
                        end
                        % save spike data
                        for uu = 1:num_sunits
                            udat{ss,1} = cat(1, udat{ss,1}, sdat(:,:,uu)./tsize.*1000);
                            udat{ss,2} = cat(1, udat{ss,2}, ...
                                cat(2, repmat([ss mm ff uu], ntr, 1), sbdat(:,uu)));
                        end
                        % save PD data
                        udat{ss,3} = cat(1, udat{ss,3}, pdat);
                        udat{ss,4} = cat(1, udat{ss,4}, ...
                            cat(2, repmat([ss mm ff], ntr, 1), pbdat));
                    end
                end
            end
        end
    end
    % save data to file
    FS_saveProjectFile('2013_LCPupil', 'pdBeep', udat);
else
    % load data from file
    udat = FS_loadProjectFile('2013_LCPupil', 'pdBeep');
end

pv   = 0.05;
y1s  = [1.2 1.2 1.2 2.6 2.6];
%hxs  = -1:0.1:1;
for ss = 1:num_sites

    % get color for this site
    [co,~] = getLCP_colorsAndSymbol(sites{ss,2}(1), sites{ss,1});

    % number of monkeys for this site
    nm = size(sites{ss,2},2);

    %% Mean±sem PD response
    axes(axs(ss)); cla reset; hold on;
    plot(tax([1 end]), [0 0], 'k:');
    plot([0 0], [-2 2], 'k:');
    nms  = nanmean(udat{ss,3}(:,:,1))';
    ses  = nanse(udat{ss,3}(:,:,1))';
    Lg   = isfinite(nms) & isfinite(ses) & sum(isfinite(udat{ss,3}(:,:,1)))'>100;
    h    = patch([pax(Lg); flipud(pax(Lg))], [nms(Lg)+ses(Lg); ...
        flipud(nms(Lg)-ses(Lg))], co{2});
    set(h,'LineStyle','none');
    plot(pax(Lg), nms(Lg), '-', 'Color', co{1});
    title(sprintf('%s', sites{ss,1}))

    % show peaks for the two monkeys
    for mm = 1:nm
        [~,sy] = getLCP_colorsAndSymbol(sites{ss,2}(mm), sites{ss,1});
        Lm  = udat{ss,4}(:,2)==mm;
        mnm = nanmean(udat{ss,3}(Lm,:,1));
        mmx = max(mnm(pax<800));
        mmi = find(mnm==mmx,1);
        h=plot(pax(mmi), mmx, sy, 'Color', co{1});
        if signtest(udat{ss,3}(Lm,mmi,1), udat{ss,4}(Lm,4))<pv
          set(h, 'MarkerFaceColor', co{1});
        end
    end
    axis([pax(1) pax(end) -0.1 y1s(ss)]);%min(nms(Lg))*1.2 max(nms(Lg))*1.2])
        
    %% baseline-subtracted Mean±sem PETH
    axes(axs(ss+num_sites)); cla reset; hold on;
    plot(tax([1 end]), [0 0], 'k:');
    plot([0 0], [-2 25], 'k:');
    spr  = udat{ss,1}-repmat(udat{ss,2}(:,5),1,ntbins);
    nms  = nanmean(spr)';
    ses  = nanse(spr)';
    Lg   = isfinite(nms) & isfinite(ses) & sum(isfinite(udat{ss,1}))'>100;
    h    = patch([tax(Lg); flipud(tax(Lg))], [nms(Lg)+ses(Lg); ...
        flipud(nms(Lg)-ses(Lg))], co{2});
    set(h,'LineStyle','none');
    plot(tax(Lg), nms(Lg), '-', 'Color', co{1});
   
    % check whether first post-beep bin is sig
    for mm = 1:nm
        [~,sy] = getLCP_colorsAndSymbol(sites{ss,2}(mm), sites{ss,1});
        Lm     = udat{ss,2}(:,2)==mm & isfinite(spr(:,tb));
        pval   = mean(spr(Lm,tb));
        h=plot(tax(tb), pval, sy, 'Color', co{1});
        if signrank(spr(Lm,tb))<pv
          set(h, 'MarkerFaceColor', co{1});
        end
        
        % correlation between baseline & baseline-subtracted evoked
        % disp([ss mm corr(spr(Lm,tb), udat{ss,2}(Lm,5), 'type', 'Spearman')])
    end
    axis([tax(1) tax(end) -1 12]);
end
axes(axs(1));
ylabel('PD (z-score)')
axes(axs(6));
xlabel('Time re: beep (ms)');
ylabel('Response (sp/s)')

%% Summary plot of correlation coefficients between trial-by-trial pupil
%   and neural reponses, per monkey and site
axes(axs(end)); cla reset; hold on;
plot([0.5 5.5], [0 0], 'k:');
pv   = 0.05;
Lpax = pax>0&pax<800;
c2dat = nans(num_sites, 2, 2, 2);
for ss = 1:num_sites

    % get matrix of PD baseline, PD response, neural response per trial
    sdat = udat{ss,4}(:,[4 4]);
    for tt = 1:size(udat{ss,3},1)
        sdat(tt,2) = max(nanrunmean(udat{ss,3}(tt,Lpax,2),50));
    end
    
    % get unique units for this site
    uus   = unique(udat{ss,2}(:,2:4),'rows'); % monk/file/unit
    nuus  = size(uus,1);
    uudat = nans(nuus, 3);
    for uu = 1:size(uus,1)
        % get corresponding spike/pupil data
        Ls = udat{ss,2}(:,2)==uus(uu,1) & ... % monk
            udat{ss,2}(:,3)==uus(uu,2) & ...  % file
            udat{ss,2}(:,4)==uus(uu,3);       % unit
        Lp = udat{ss,4}(:,2)==uus(uu,1) & ... % monk
            udat{ss,4}(:,3)==uus(uu,2);       % file
        
        tmpdat      = cat(2, sdat(Lp,:), udat{ss,1}(Ls,tb));
        tmpdat      = tmpdat(isfinite(sum(tmpdat(:,1:3),2)),:);
        [Rho,Pval]  = partialcorr(tmpdat, 'type', 'Spearman');
        uudat(uu,:) = [Rho(2,3) Pval(2,3) uus(uu,1)];
        
%         subplot(2,1,1); cla reset; hold on;
%         plot(tmpdat(:,2), tmpdat(:,3), 'k.');
%         subplot(2,1,2); cla reset; hold on;
%         plot(tmpdat(:,1), tmpdat(:,3), 'k.');
%         subplot(2,1,2); cla reset; hold on;
%         plot(tmpdat(:,1), tmpdat(:,3), 'k.');
    end
    uudat = uudat(isfinite(uudat(:,1)),:);
    
    % plot per pair, monk
    %   1: pupil baseline vs pupil response, controlling for neural response
    %   2: pupil baseline vs neural response, controlling for pupil response
    %   3: pupil response vs neural response, controlling for pupil baseline
    for mm = 1:2
        [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}{mm}, sites{ss,1});
        Lm = uudat(:,3)==mm;
        xb = ss+0.3*(mm-1);
        xs = xb+rand(sum(Lm),1)*0.2-0.1;
        ys = uudat(Lm,1);
        Llo = ys<0;
        Lhi = ys>0;
        plot(xs(Lhi),ys(Lhi),sy,'Color',co{1});
        plot(xs(Llo),ys(Llo),sy,'Color',co{2});
        Lp = uudat(Lm,2)<pv;
        plot(xs(Lhi&Lp), ys(Lhi&Lp), sy, 'Color', co{1}, 'MarkerFaceColor', co{1});
        plot(xs(Llo&Lp), ys(Llo&Lp), sy, 'Color', co{2}, 'MarkerFaceColor', co{2});
        
        h  = plot(xb+[-0.1 0.1], median(ys).*[1 1], 'k-');
        sr = signrank(ys);
        if sr<pv
            set(h, 'LineWidth', 3);
        end
        
        c2dat(ss,mm,:,1) = [sum(ys>0&Lp) sum(isfinite(ys))];
        c2dat(ss,mm,:,2) = [sum(ys<0&Lp) sum(isfinite(ys))];
        text(xb-.15,0.75,sprintf('%d/%d=%d',sum(ys>0&Lp),sum(isfinite(ys)),...
            round(sum(ys>0&Lp)/sum(isfinite(ys))*100)));
        text(xb-.15,-0.75,sprintf('%d/%d=%d',sum(ys<0&Lp),sum(isfinite(ys)),...
            round(sum(ys<0&Lp)/sum(isfinite(ys))*100)));

%         % print stats
%         disp(sprintf('site=%3s, monk=%7s, medp=%.3f, tot=%2d, +=%2d, -=%2d', ...
%             sites{ss,1}, sites{ss,2}{mm}, sr, ...
%             sum(isfinite(ys)), sum(ys>0&Lp), sum(ys<0&Lp)))
    end
end
axis([0.8 5.6 -1 1]);

% panel labels
setPLOT_panelLabel(axs( 1), 1);
setPLOT_panelLabel(axs( 6), 2);
setPLOT_panelLabel(axs(11), 3);

axes(axs(end));
set(gca, 'XTick', []);
ylabel('Correlation coefficient')
xlabel('Monkey')


% chi-2 for LC/IC/SC proportions, Sprout ACC/PCC
inds = [1 2 1; 1 2 2; 1 3 1; 1 3 2; 2 3 1; 2 3 2; 4 5 1];
js   = {'Pos' 'Neg'};
for ii = 1:size(inds,1)
    for jj = 1:2
        n1 = c2dat(inds(ii,1),inds(ii,3),1,jj);
        N1 = c2dat(inds(ii,1),inds(ii,3),2,jj);
        
        n2 = c2dat(inds(ii,2),inds(ii,3),1,jj);
        N2 = c2dat(inds(ii,2),inds(ii,3),2,jj);
        
        p0 = (n1+n2)/(N1+N2);
        n10 = N1*p0;
        n20 = N2*p0;
        observed = [n1 N1-n1 n2 N2-n2];
        expected = [n10 N1-n10 n20 N2-n20];
        chi2stat = sum((observed-expected).^2 ./ expected);
        p = 1 - chi2cdf(chi2stat,1);
        disp(sprintf('%6s, %3s vs %3s, %s, p=%.2f', sites{inds(ii,1),2}{inds(ii,3)}, ...
            sites{inds(ii,1),1}, sites{inds(ii,2),1}, js{jj}, p))
    end
end

% %% Summary box-whisker of per-session correlation coefficients 
% %   of PD response vs neural response magnitudes
% axes(axs(end)); cla reset; hold on;
% plot([0.5 5.5], [0 0], 'k:');
% for ss = 1:num_sites
% 
%     % get per monkey residuals to PD baseline vs peak
%     PC = nans(size(udat{ss,3},1),1);
%     for mm = 1:size(sites{ss,2},2)
%         Lm  = udat{ss,4}(:,2)==mm;
%         % get peak from mean PD response across sessions
%         mnm = nanmean(udat{ss,3}(Lm,:,1));
%         mmx = max(mnm(pax>100&pax<800));
%         mmi = find(pax>100,1)+find(mnm(pax>100&pax<800)==mmx,1)-1;
%         B   = nanmean(udat{ss,3}(Lm,mmi-50:mmi+50,1),2);
%         A   = [udat{ss,4}(Lm,4) ones(sum(Lm),1)];
%         Lg  = isfinite(A(:,1)) & isfinite(B);
%         X   = lscov(A(Lg,:),B(Lg));
%         PC(Lm) = B-A*X;        
%     end
%     
%     % get unique units for this site
%     uus   = unique(udat{ss,2}(:,2:4),'rows'); % monk/file/unit
%     nuus  = size(uus,1);
%     uudat = nans(nuus, 3);
%     for uu = 1:size(uus,1)
%         % get corresponding spike/pupil data
%         Ls = udat{ss,2}(:,2)==uus(uu,1) & ...
%             udat{ss,2}(:,3)==uus(uu,2) & ...
%             udat{ss,2}(:,4)==uus(uu,3);
%         Lp = udat{ss,4}(:,2)==uus(uu,1) & ...
%             udat{ss,4}(:,3)==uus(uu,2);
%         %sp = nanmax(udat{ss,1}(Ls,tb:tb+100),[],2)-udat{ss,2}(Ls,5);
%         sp = udat{ss,1}(Ls,tb);%-udat{ss,2}(Ls,5);
%         pd = PC(Lp);
%         Lg = isfinite(sp) & isfinite(pd);
%         if sum(Lg) > 10
%             % compute correlation between PD resids, spike response
%             % for each unit monk/file/unit
%             [R,P] = corr(sp(Lg), pd(Lg), 'type', 'Spearman');
%             uudat(uu,:) = [R P uus(uu,1)];
%         end
%     end
%     uudat = uudat(isfinite(uudat(:,1)),:);
%     
%     % plot per monk
%     for mm = 1:2
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}{mm}, sites{ss,1});
%         Lm = uudat(:,3)==mm;
%         xb = ss+0.3*(mm-1);
%         xs = xb+rand(sum(Lm),1)*0.2-0.1;
%         ys = uudat(Lm,1);
%         plot(xs,ys,sy,'Color',co{1});
%         Lp = uudat(Lm,2)<pv;
%         plot(xs(Lp), ys(Lp), sy, 'Color', co{1}, 'MarkerFaceColor', co{1});
%         h = plot(xb+[-0.1 0.1], median(ys).*[1 1], 'k-');
%         sr = signrank(ys);
%         if sr<pv
%             set(h, 'LineWidth', 3);
%         end
%         
%         % print stats
%         disp(sprintf('site=%3s, monk=%7s, medp=%.3f, tot=%2d, +=%2d, -=%2d', ...
%             sites{ss,1}, sites{ss,2}{mm}, sr, ...
%             sum(isfinite(ys)), sum(ys>0&Lp), sum(ys<0&Lp)))
%     end
% end

%% Junk below
% 
%     vals = cat(1, vals, uudat(isfinite(uudat(:,1)),[1 3]));
% end
% 
% pv = 0.05;
% xs = unique(vals(:,2));
% axes(axs(end)); cla reset; hold on;
% plot([0 10], [0 0], 'k:');
% h = boxplot(vals(:,1), vals(:,2), 'notch', 'on', ...
%     'positions', xs, 'symbol', '');
% for ii = 1:size(h,2)
%     set(h(:,ii), 'Color', cols{ii,1});
%     Lp = vals(:,2)==xs(ii);
%     if signrank(vals(Lp,1))<pv
%         plot(xs(ii), nanmedian(vals(Lp,1)), '.', ...
%             'Color', cols{ii,1}, 'MarkerSize', 25)
%     end
% end
% xlim([0.5 5.5]);
% 
%     
%     % single-unit values, ps
%     vs  = uudat(:,1);
%     Lps = uudat(:,2)<0.05;
%     
%     Na = hist(vs,hxs);
%     H = bar(hxs,Na,0.4);
%     set(H, 'EdgeColor', 'k', 'FaceColor', 'w');
%     
%     Np = hist(vs(Lps),hxs);
%     H = bar(hxs,Np,0.4);
%     set(H, 'EdgeColor', 'k', 'FaceColor', 'k');
%     xlim([-1 1]);
%     
%     mv = max(Na);
%     for mm = 1:nm
%         vals = vs(uudat(:,3)==mm & isfinite(vs));
%         h=plot(nanmedian(vals), mv+mm, 'o', 'Color', co{mm});
%         sr = signrank(vals);
%         if sr < pv
%             set(h, 'MarkerFaceColor', co{mm});
%         end
%         disp(sprintf('%3s, %7s: %3d tot, %3d+, %3d-, med=%2.2f (p=%.3f)', ...
%             sites{ss,1}, sites{ss,2}{mm}, length(vals), sum(Lps&vs>0), sum(Lps&vs<0), ...
%             nanmedian(vals), sr))
%     end
%     xlim([-1 1])
% end


%% JUNK BELOW

%     %% Try correlating time of PETH peak and time of PD peak
%      % get unique units for this site
%     cuudat = nans(nuus, 3);
%     for uu = 1:size(uus,1)
%         % get corresponding spike/pupil data
%         Fs = find(udat{ss,2}(:,2)==uus(uu,1) & ...
%             udat{ss,2}(:,3)==uus(uu,2) & ...
%             udat{ss,2}(:,4)==uus(uu,3));
%         Fp = find(udat{ss,4}(:,2)==uus(uu,1) & ...
%             udat{ss,4}(:,3)==uus(uu,2));
%         
%         spdat = nans(length(Fs),2);
%         for ii = 1:length(Fs)
%             sp  = udat{ss,1}(Fs(ii),tb:tb+100)-udat{ss,2}(Fs(ii),5);
%             if any(isfinite(sp))
%                 spdat(ii,1) = find(sp==nanmax(sp),1);
%                 pd = udat{ss,3}(Fp(ii),pb:end);
%                 if any(isfinite(pd))
%                     spdat(ii,2) = find(pd==nanmax(pd),1);
%                 end
%             end
%         end
%         Lg = isfinite(spdat(:,1)) & isfinite(spdat(:,2));
%         if sum(Lg) > 10
%             % compute correlation between PD resids, spike response
%             [R,P] = corr(spdat(Lg,1), spdat(Lg,2), 'type', 'Spearman');
%             cuudat(uu,:) = [R P uus(uu,1)];
%         end
%     end
%     vs = cuudat(:,1);
%     for mm = 1:nm
%         vals = vs(cuudat(:,3)==mm & isfinite(vs));
%         disp(sprintf('%3s, %7s: %3d tot, %3d+, %3d-, med=%2.2f (p=%.3f)', ...
%             sites{ss,1}, sites{ss,2}{mm}, length(vals), sum(Lps&vs>0), sum(Lps&vs<0), ...
%             nanmedian(vals), signtest(vals)))
%     end
% end

