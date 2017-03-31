function figLCP_tonic(num, collectData)
% function figLCP_tonic(num, collectData)
%
% Trial-by-trial relationship between tonic neural activity in LC/IC/SC
% (other?) and baseline PD. Columns are brain regions (LC/IC/SC).
%   Top row is mean PD per trial within an example session.
%   Next row is mean spike rate per trial within the same session.
%   Next row is spike rate versus PD for that session.
%   Bottom row is histogram of correlation coefficients from each session.
%       Black/gray are two monkeys.
%       Filled are p<0.01 (I?ll probably redo the stats using shuffled trials).
%
% Uses siteData:
%  siteData{1}: trialsxcols matrix, cols are:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correcdt
%   4 ... beep on time, wrt to fix start time
%   5 ... trial begin time, wrt to fix start time
%   6 ... trial end time, wrt to fix start time
%   7 ... trial wrt time (cpu clock)
%   8 ... LFP index corresponding to fix start time (coded above)
%   9 ... ELESTM on time (when appropriate), wrt fix start time
%
%  siteData{2}: Analog:
%   dim1: trial
%   dim2: sample
%   dim3: 1=x, 2=y, 3=z-pupil, 4=corrected z-pupil, 5=slope
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
%% Parse args
if nargin < 1 || isempty(num)
    num = 2;
end

if nargin < 2 || isempty(collectData)
    collectData = false;
end

%% Set up figure
wid       = 17.6; % total width
hts       = [1.5 1.5 3 3];
cols      = {5,5,5,1};
[axs,~] = getPLOT_axes(num, wid, hts, cols, 2, .5, [], 'Joshi et al', true);
set(axs,'Units','normalized');

% gr = 0.7.*ones(1,3);
% site, monks, example sessions (monk,file,single unit,ylim1,ylim2)
sites = {...
    'LC',  {'Oz' 'Cicero'},      [1 19 1  5  3]; ...
    'IC',  {'Oz' 'Cicero'},      [2 24 2 10 10]; ...
    'SC',  {'Oz' 'Cicero'},      [2  6 2 25 12]; ...
    'ACC', {'Sprout' 'Atticus'}, [2  3 1 30 20]; ... % 1 5 1; 1 6 1; 2 2 1
    'PCC', {'Sprout' 'Cheetah'}, [1  10 1 10 8]};
num_sites = size(sites, 1);

%% Collect data
if collectData
    udat = cell(num_sites, 2);
    for ss = 1:num_sites
        for mm = 1:size(sites{ss,2},2)
            % load file
            [base_dir, fnames] = getLCP_cleanDataDir(sites{ss,2}{mm}, sites{ss,1});
            if ~isempty(fnames)
                for ff = 1:length(fnames)
                    disp(sprintf('monkey = %s, site = %s, %d/%d files', ...
                        sites{ss,2}{mm}, sites{ss,1}, ff, length(fnames)))
                    load(fullfile(base_dir, fnames{ff}));
                    
                    % only single units
                    num_units  = size(siteData{3}, 2)-1;
                
                    if num_units > 0
                        
                        % Use only no-beep, no-stim trials
                        Fnb = find(~isfinite(siteData{1}(:,4)) & ...
                            ~isfinite(siteData{1}(:,9)));
                        num_tr     = length(Fnb);
                        tr_ends    = siteData{1}(Fnb,6)-200;
                        sdat       = nans(num_tr, num_units);
                        pdat       = nans(num_tr, 1);
                        
                        for tt = 1:num_tr
                            
                            % get pupil average(s)
                            pi = max(500, find(isfinite(siteData{2}(Fnb(tt),:,4)),1));
                            pe = find(isfinite(siteData{2}(Fnb(tt),:,4)),1,'last');
                            
                            if pe > pi + 1000
                                
                                % mean PD during trial
                                pdat(tt) = mean(siteData{2}(Fnb(tt),pi:pe,4));
                                
                                for uu = 1:num_units
                                    sp = siteData{3}{Fnb(tt), uu+1};
                                    sdat(tt,uu) = ...
                                        sum(sp>=200&sp<=tr_ends(tt))./...
                                        (tr_ends(tt)-200).*1000;
                                end
                            end
                        end
                        
                        % compute correlations
                        Lpg = isfinite(pdat);
                        for uu = 1:num_units
                            Lu = Lpg & isfinite(sdat(:,uu));
                            if sum(Lu) > 10
                                
                                % correlation of residuals after linear
                                % regressions
                                tax = siteData{1}(Fnb(Lu),7);
                                tax = tax - tax(1);                  
            
                                X = [ones(sum(Lu),1) tax];
                                % X = [ones(sum(Lu),1) (1:sum(Lu))'];
                                % [Bp,BINTp,Rp] = regress(pdat(Lu), X);
                                Bp = robustfit(X(:,2), pdat(Lu));
                                Rp = pdat(Lu)-X*Bp;
                                Bs = robustfit(X(:,2), sdat(Lu,uu));
                                Rs = sdat(Lu,uu)-X*Bs;
                                [R,P] = corr(Rp, Rs, 'type', 'Spearman');
                                udat{ss,1} = cat(1,udat{ss},[R P mm ff uu]);

                                % check for example site
                                if mm==sites{ss,3}(1) && ...
                                        ff==sites{ss,3}(2) && ...                                
                                        uu==sites{ss,3}(3)
                                    
                                    udat{ss,2} = cat(2, ...
                                        tax, ...
                                        pdat(Lu), ...
                                        Rp, ...
                                        sdat(Lu,uu), ...
                                        Rs);
                                end                                
                            end
                        end
                    end
                end
            end
        end
    end
    % save data to file
    FS_saveProjectFile('2013_LCPupil', 'Tonic', udat);
else
    
    % load data from file
    udat = FS_loadProjectFile('2013_LCPupil', 'Tonic');
end

%% Example plotz
for ss = 1:num_sites
    
    % get colors, time axis
    [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(sites{ss,3}(1)), sites{ss,1});
    tax = udat{ss,2}(:,1)./1000./60; % (1:size(udat{ss,2},1))'; % udat{ss,2}(:,1)
    
    %   1: pupil
    axes(axs(ss)); cla reset; hold on;
    plot(tax, udat{ss,2}(:,2), '.', 'Color', co{1}, 'MarkerSize', 3);    
    lsline;
    axis([min(tax) max(tax) -3 3]);
    title(sprintf('%s', sites{ss,1}))

    %   2: spikes
    axes(axs(ss+num_sites)); cla reset; hold on;
    plot(tax, udat{ss,2}(:,4), '.', 'Color', co{1}, 'MarkerSize', 2);
    lsline;
    axis([min(tax) max(tax) 0 sites{ss,3}(4)]);

    %   3: residual pupil vs spike
    axes(axs(ss+2*num_sites)); cla reset; hold on;
    plot(udat{ss,2}(:,3), udat{ss,2}(:,5), '.', 'Color', co{1}, 'MarkerSize', 2);    
    lsline;
    plot([-30 30], [0 0], 'k:'); 
    plot([0 0], [-30 30], 'k:'); 
    axis([-3 3 -sites{ss,3}(5) sites{ss,3}(5)]);
    
    % report partial correlation, to add to fig
    [R,P] = corr(udat{ss,2}(:,3), udat{ss,2}(:,5), 'type', 'Spearman');
    disp(sprintf('site=%s, r=%.2f, p=%.5f', sites{ss,1}, R, P))
end

%% Summary scatter plotz
pv = 0.05;
axes(axs(end)); cla reset; hold on;
plot([0 10], [0 0], 'k:');
c2dat = nans(num_sites, 2, 2, 2);
for ss= 1:num_sites
    for mm = 1:size(sites{ss,2},2)
        [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(mm), sites{ss,1});
        xb = ss+0.3*(mm-1);
        Lm = udat{ss,1}(:,3) == mm;
        ys = udat{ss,1}(Lm,1);        
        xs = xb+rand(sum(Lm),1)*0.2-0.1;
        
        % r is positive/negative
        for pp = 1:2
            if pp == 1
                Lpp = ys>=0;
            else
                Lpp = ys<0;
            end
            
            plot(xs(Lpp), ys(Lpp), sy, 'Color', co{pp});
            Lp = udat{ss,1}(Lm,2)<pv;
            plot(xs(Lpp&Lp), ys(Lpp&Lp), sy, 'Color', co{pp}, 'MarkerFaceColor', co{pp});
        end
        
        % display median, sig test
        h  = plot(xb+[-0.1 0.1], median(ys).*[1 1], 'k-');
        sr = signrank(ys);
        if sr<pv
            set(h, 'LineWidth', 3);
        end
        
        % display, save counts
        c2dat(ss,mm,:,1) = [sum(ys>0&Lp) sum(isfinite(ys))];
        c2dat(ss,mm,:,2) = [sum(ys<0&Lp) sum(isfinite(ys))];
        text(xb-.15,0.55,sprintf('%d/%d=%d',sum(ys>0&Lp),sum(isfinite(ys)),...
            round(sum(ys>0&Lp)/sum(isfinite(ys))*100)));
        text(xb-.15,-0.55,sprintf('%d/%d=%d',sum(ys<0&Lp),sum(isfinite(ys)),...
            round(sum(ys<0&Lp)/sum(isfinite(ys))*100)));
        
        % print stats
        disp(sprintf('site=%3s, monk=%7s, medp=%.3f, tot=%2d, +=%2d, -=%2d', ...
            sites{ss,1}, sites{ss,2}{mm}, sr, ...
            sum(isfinite(ys)), sum(ys>0&Lp), sum(ys<0&Lp)))
        
        % label example session
        if mm == sites{ss,3}(1)
            Le = ys == corr(udat{ss,2}(:,3), udat{ss,2}(:,5), 'type', 'Spearman');
            plot(xs(Le), ys(Le), sy, 'Color', 'k');
        end
        
    end    
end
axis([0.7 5.3 -0.6 0.6]);

axs_pl = 1:ss:length(axs);
for xx = 1:length(axs_pl)
    setPLOT_panelLabel(axs(axs_pl(xx)), xx);
end
axes(axs(1));
ylabel('Pupil diameter (z-score)');
axes(axs(1+num_sites));
ylabel('Spike rate (sp/s');
xlabel('Time from session onset (min)');
axes(axs(1+num_sites*2));
ylabel('Spike rate (sp/s');
xlabel('Pupil diameter (z-score)');
axes(axs(1+num_sites*3));
ylabel('Partial Spearman''s correlation')
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

%% Summary box-n-whisker plot
% first collect data
% vals = [];
% cols = {};
% for ss= 1:num_sites
%     for mm = 1:size(sites{ss,2},2)
%         % get correlation coefficients
%         Lm = udat{ss,1}(:,3) == mm;        
%         vals = cat(1, vals, [udat{ss,1}(Lm,1) ...
%             (ss+0.2*(mm-1)).*ones(sum(Lm),1)]);
%         
%         % get colors
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(mm), sites{ss,1});
%         cols = cat(1, cols, co);
%     end
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
% %% Show stats per site/monkey
% for ss = 1:size(sites,1)
%     vals = udat{ss}(:,1);
%     ps   = udat{ss}(:,2);
%     disp(sprintf('%s/%6s: med=%.2f (p=%.2f), %3d total, %2d-, %2d+', ...
%         sites{ss,1}, sites{ss,2}, ...
%         nanmedian(vals), signtest(vals), ...
%         sum(isfinite(vals)), sum(vals<0&ps<0.05), sum(vals>0&ps<0.05)))
% end
% 
% % summary histograms
% for ss = 1:size(sites,1)
%     figure(fig(sites{ss,3}));
%     axes(axs(sites{ss,3},9+(sites{ss,5}-1)*3+sites{ss,4})); cla reset; hold on;
%         
%     % single-unit values, ps
%     vs  = udat{ss}(:,1);
%     Lps = udat{ss}(:,2)<0.05;
%         
%     N = hist(vs,xs);
%     % plot(nanmedian(vs).*[1 1], [0 max(N)+1], 'k:');
%     H = bar(xs,N,0.3);
%     set(H, 'EdgeColor', 'k', 'FaceColor', 'w');
%         
%     N = hist(vs(Lps),xs);
%     H = bar(xs,N,0.3);
%     set(H, 'EdgeColor', 'k', 'FaceColor', 'k');
%     xlim([-1 1]);
% end

%% BELOW FOR HISTOGRAM
%scales = [-1 1];
% mn = 0;
%         vs = cat(1, udat{1,ss,uu,pp}(:,1), udat{2,ss,uu,pp}(:,1));
%         ps = cat(1, udat{1,ss,uu,pp}(:,2), udat{2,ss,uu,pp}(:,2));
%         Lg = isfinite(vs);
%         vs = vs(Lg);
%         ps = ps(Lg);
%         N = hist(vs,xs);
%         if max(N)>mn
%             mn = max(N)+2;
%         end
%         H = bar(xs,N.*scales(uu));
%         set(H,'FaceColor','none')
%         N = hist(vs(ps<0.05),xs);
%         H = bar(xs,N.*scales(uu));
%         %        title(sprintf('%s/%s, Median=%.2f, p=%.2f', ...
%         %            monkeys{mm}, sites{ss}, median(vals), signtest(vals)))
%     end
%     axis([-1 1 -mn mn]);
% end


% uu=1; % 1=multi, 2=single units
% num_sites=3;
% xs = -1:0.1:1;
% for mm = 1:num_monkeys
%     for ss = 1:num_sites
%         subplot(num_monkeys,num_sites,(mm-1)*num_sites+ss); cla reset; hold on;
%         vals = udat{mm,ss,uu}(:,1);
%         ps   = udat{mm,ss,uu}(:,2);
%         Lg   = isfinite(vals);
%         vals = vals(Lg);
%         ps   = ps(Lg);
%         N = hist(vals,xs);
%         H = bar(xs,N);
%         set(H,'FaceColor','none')
%         N = hist(vals(ps<0.05),xs);
%         H = bar(xs,N);
%         xlim([-1 1]);
%         title(sprintf('%s/%s, Median=%.2f, p=%.2f', ...
%             monkeys{mm}, sites{ss}, median(vals), signtest(vals)))
%     end
% end

