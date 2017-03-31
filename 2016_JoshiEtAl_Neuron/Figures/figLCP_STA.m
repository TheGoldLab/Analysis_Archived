function figLCP_STA(num, ptype, collectData)
% function figLCP_STA(num, ptype, collectData)
%
% STA = spike-triggered average (of pupil-diameter derivative)
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
%   9 ... ELESTM on time (when appropriate), wrt fix start time

%  siteData{2}: Analog:
%   dim1: trial
%   dim2: sample
%   dim3: 1 = x, 2 = y, 3 = pupil
%   [remember first sample is considered time=0 (i.e., wrt fix start time)]
%     = eyedat(Lgood,:,:);
%
%  siteData{3}: spikes, re-coded wrt fix start time
%
%  siteData{4}: LFP

%% Parse args
if nargin < 1 || isempty(num)
    num = 4;
end

if nargin < 2 || isempty(ptype)
    ptype = 'slope'; % 'pupil' or 'slope'
end

if nargin < 3 || isempty(collectData)
    collectData = false;
end

%% Set up figure
wid     = 17.6; % total width
cols    = {5,5,5,5};
hts     = repmat(2.5, 1, length(cols));
[axs,~] = getPLOT_axes(num, wid, hts, cols, 1.1, .5, [], 'Joshi et al', true);
set(axs,'Units','normalized');

%% data to collect/include
% gr = 0.7.*ones(1,3);
% site    monk                  example
sites        = {...
    'LC'    {'Oz' 'Cicero'}         [1 10]; ...
    'IC'    {'Oz' 'Cicero'}         [1 11]; ...
    'SC'    {'Oz' 'Cicero'}         [1 6]; ...
    'ACC'   {'Sprout' 'Atticus'}    [2 4]; ...
    'PCC'   {'Sprout' 'Cheetah'}    [2 6]; ...
    'MT'    {'Cheetah'}             []; ...
    };
num_sites    = size(sites, 1);
num_shuffled = 100;

%% Possibly collect data
if collectData
    udat = cell(num_sites, 2, 2);
    for ss = 1:num_sites
        for mm = 1:size(sites{ss,2},2)
            % load file
            [base_dir, fnames] = getLCP_cleanDataDir(sites{ss,2}{mm}, sites{ss,1});
            if ~isempty(fnames)
                for ff = 1:length(fnames)
                    disp(sprintf('monkey = %s, site = %s, %d/%d files', ...
                        sites{ss,2}{mm}, sites{ss,1}, ff, length(fnames)))
                    load(fullfile(base_dir, fnames{ff}));
                    %  sdat is:pupil measure
                    %   1: per unit (multi first, if appropriate)
                    %   2: NUM_BINS
                    %   3: 1-mean, 2-sem, 3-shuffled mean, 4=ps
                    [sdat,sudat] = getLCP_stpds(siteData, ptype, 'single', num_shuffled);
                    if ~isempty(sdat)
                        % means, shuffled means, ps
                        udat{ss,mm,1} = cat(1, udat{ss,mm,1}, sdat);
                        % list of file index, unit indices
                        udat{ss,mm,2} = cat(1, udat{ss,mm,2}, ...
                            cat(2, ff.*ones(size(sudat)),sudat));
                    end
                end
            end
        end
    end
    % save data to file
    FS_saveProjectFile('2013_LCPupil', sprintf('STA_%s', ptype), udat);
else
    % load data from file
    udat = FS_loadProjectFile('2013_LCPupil', sprintf('STA_%s', ptype));
end

%% per site summary
ns   = 5; % ignore MT
xax  = (-1000:2000)';
emx  = [-0.0004 0.001];
gr   = 0.8*ones(1,3);
maxt = 0;
pv   = 0.05;
for ss = 1:ns
    for mm = 1:2
        if size(udat{ss,mm,1},1) > maxt
            maxt = size(udat{ss,mm,1},1);
        end
    end
end

for ss = 1:ns
    
    % get color, symbols
    [co, ~] = getLCP_colorsAndSymbol(sites{ss,2}{sites{ss,3}(1)}, sites{ss,1});
   
    %% Examples
    axes(axs(ss)); cla reset; hold on;
    plot(xax([1 end]), [0 0], 'k:');
    plot([0 0], emx, 'k:');
    plot(xax, udat{ss,sites{ss,3}(1)}(sites{ss,3}(2),:,3), '-', ...
        'Color', gr);
    plot(xax, udat{ss,sites{ss,3}(1)}(sites{ss,3}(2),:,1), '-', ...
        'Color', co{1});
    axis([xax(1) xax(end) -0.0004 0.001]);
    title(sprintf('%s', sites{ss,1}));
    if ss > 1
        set(gca, 'YTickLabel', []);
    end
    set(gca,'XTickLabel', '');
    
    %% Summary mean±sem realsies-shuffled combined data from both monkeys
    axes(axs(ns+ss)); cla reset; hold on;
    plot(xax([1 end]), [0 0], 'k:');
    plot([0 0], emx, 'k:');
    vals  = cat(1, ...
        udat{ss,1,1}(:,:,1)-udat{ss,1,1}(:,:,3), ...
        udat{ss,2,1}(:,:,1)-udat{ss,2,1}(:,:,3));
    mn    = nanmean(vals)';
    se    = nanse(vals)';
    h=patch([xax; flipud(xax)], [mn-se; flipud(mn+se)], co{2});
    set(h, 'LineStyle', 'none');
    plot(xax, nanmean(vals), '-', 'Color', co{1});
    axis([xax(1) xax(end) -0.00015 0.00025]);
    if ss > 1
        set(gca, 'YTickLabel', []);
    end
    set(gca,'XTickLabel', '');
    
    % label peak
    peaki = find(mn==max(mn),1);
    h=text(450, 0.0002, sprintf('%d ms', xax(peaki)));
    set(h,'FontSize',9);
    % bootstrap!
    N   = 10000;
    bsd = normrnd(repmat(vals(:,peaki),1,N), ...
        repmat(cat(1,udat{ss,1,1}(:,peaki,2),udat{ss,2,1}(:,peaki,2)),1,N));
    if prctile(median(bsd),2.5)>0
        set(h, 'FontWeight', 'bold');
    end
        
    %% Per monk color plotz
    for mm = 1:2
        
        %% Color plotz -- for realsies
        axes(axs(ns*(1+mm)+ss)); cla reset; hold on;
        vals  = udat{ss,mm,1}(:,:,1)-udat{ss,mm,1}(:,:,3);
        [~,I] = sort(max(vals,[],2)-min(vals,[],2));
        vals  = vals(I,:);
        nv    = size(vals,1);
        imagesc(xax,-nv:-1,vals,[-0.0004 0.0004]);
        hold on;
        plot([0 0], [-maxt 1], 'k:', 'LineWidth', 1);
        axis([xax(1) xax(end) -maxt -1]);
        if ss > 1
            set(gca, 'YTickLabel', []);
        end
        if mm == 1
            set(gca, 'XTickLabel', []);
        end
    end
end

% axis labels
axes(axs(1));
ylabel('Change in PD (z/ms)')
axes(axs(16));
ylabel('Unit number')
xlabel('Time re:spike (ms)');

%% Add labels, stats 
tdat = cell(ns,2);
for ss = 1:ns
    for mm = 1:2
        tdat{ss,mm} = nans(size(udat{ss,mm,1},1),2);
        for uu = 1:size(udat{ss,mm,1},1)
            % find longest contiguous + peak
            vals = udat{ss,mm,1}(uu,:,1)-udat{ss,mm,1}(uu,:,3);
            Lp = (vals>0 & udat{ss,mm,1}(uu,:,4)<pv)';
            d = diff([false;Lp;false]);
            p = find(d==1);
            q = min(length(xax),find(d==-1));
            Lgood = xax(q)>-100&xax(p)<700&(q-p)>=75;
            
%             cla reset; hold on;
%             Lp=udat{ss,mm,1}(uu,:,4)<pv;
%             plot(xax, vals, 'k-');
%             plot(xax(Lp), vals(Lp), 'r.');

            if any(Lgood)
                p = p(Lgood);
                q = q(Lgood);
                for pp = 1:sum(Lgood)
                     pvals = smooth(vals(p(pp):q(pp)-1),100,'sgolay');
                     maxv  = max(pvals);
                     maxt  = xax(p(pp)-1+find(pvals==maxv,1));
                     if maxt > 0 
                         break;
                     end
                end
                 
                %plot(xax(p(ix):q(ix)-1), pvals, 'y-');
                %plot(xax(maxi), max(pvals), 'go');
                
                tdat{ss,mm}(uu,:) = [maxt maxv];
            end
            %axis([-1000 2000 -0.0004 0.001]);
            %r = input('next')
        end
    end
end

% For anova
vals   = [];
g_site = {};
g_monk = {};
c2dat  = nans(ns,2,2);
for ss = 1:ns
    for mm = 1:2
        axes(axs(ns*(1+mm)+ss));
        Lg   = isfinite(tdat{ss,mm}(:,1));
        vals = cat(1,vals,tdat{ss,mm}(Lg,1));
        tmp  = cell(sum(Lg),1);
        [tmp{:}] = deal(sites{ss,1});
        g_site = cat(1,g_site,tmp);
        [tmp{:}] = deal(sites{ss,2}{mm});
        g_monk = cat(1,g_monk,tmp);
        
        h=title(sprintf('%s:%d/%d=%dp,t=%d ms', ...
            sites{ss,2}{mm}(1:2), sum(Lg), size(Lg,1), ...
            round(sum(Lg)./size(Lg,1).*100), ...
            round(median(tdat{ss,mm}(Lg,1)))));
        set(h,'FontWeight','normal', 'FontSize', 9);
        c2dat(ss,mm,:) = [sum(Lg) size(Lg,1)];
    end
end

Lg = ismember(g_monk, {'Oz', 'Cicero'});
p=anovan(vals(Lg), {g_site(Lg),g_monk(Lg)},'display','off');
disp(p)

%% chi-2 for LC/IC/SC proportions, Sprout ACC/PCC
inds = [1 2 1; 1 2 2; 1 3 1; 1 3 2; 2 3 1; 2 3 2; 4 5 1];
for ii = 1:size(inds,1)
    
    n1 = c2dat(inds(ii,1),inds(ii,3),1);
    N1 = c2dat(inds(ii,1),inds(ii,3),2);

    n2 = c2dat(inds(ii,2),inds(ii,3),1);
    N2 = c2dat(inds(ii,2),inds(ii,3),2);

    p0 = (n1+n2)/(N1+N2);
    n10 = N1*p0;
    n20 = N2*p0;
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    chi2stat = sum((observed-expected).^2 ./ expected);
    p = 1 - chi2cdf(chi2stat,1);
    disp(sprintf('%6s, %3s vs %3s, p=%.2f', sites{inds(ii,1),2}{inds(ii,3)}, ...
        sites{inds(ii,1),1}, sites{inds(ii,2),1}, p))
end

%
% %% Combined per monk
% ns   = 5; % ignore MT
% xax  = (-1000:2000)';
% %mmx = [-0.00027 0.00033];
% emx  = [-0.0004 0.001];
% pv   = 0.05;
% maxt = 0;
% for ss = 1:ns
%     for mm = 1:2
%         if size(udat{ss,mm,1},1) > maxt
%             maxt = size(udat{ss,mm,1},1);
%         end
%     end
% end
% 
% % for bootstrap
% N    = 1000;
% Lx   = xax<1500;
% maxs = nans(N,ns,2);
% mins = nans(N,ns,2);
% 
% for ss = 1:ns
%         
%     %% Examples
%     axes(axs(ss)); cla reset; hold on;
%     [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}{sites{ss,3}(1)}, sites{ss,1});
%     plot(xax([1 end]), [0 0], 'k:');
%     plot([0 0], emx, 'k:');
%     plot(xax, udat{ss,sites{ss,3}(1)}(sites{ss,3}(2),:,1), '-', ...
%         'Color', co{1});
%     plot(xax, udat{ss,sites{ss,3}(1)}(sites{ss,3}(2),:,3), '-', ...
%         'Color', 0.8*ones(1,3));
%     Lp = udat{ss,sites{ss,3}(1),1}(sites{ss,3}(2),:,4)<pv;
%     plot(xax(Lp), emx(2), 'k.');
%     axis([xax(1) xax(end) -0.0004 0.001]);
%     title(sprintf('%s', sites{ss,1}));
%     if ss > 1
%         set(gca, 'YTickLabel', []);
%     end
%     set(gca,'XTickLabel', '');
%     
%     %% Fraction significant
%     axes(axs(ns+ss)); cla reset; hold on;
%     plot([0 0], [0 1], 'k:');    
%     ps = cat(1, udat{ss,1,1}(:,:,4), udat{ss,2,1}(:,:,4));
%     vs = cat(1, udat{ss,1,1}(:,:,1), udat{ss,2,1}(:,:,1));
%     plot(xax, sum(ps<pv&vs<0)./sum(isfinite(ps)), '-', 'Color', co{2});
%     plot(xax, sum(ps<pv&vs>0)./sum(isfinite(ps)), '-', 'Color', co{1});
%     axis([xax(1) xax(end) 0 0.5]);
%     if ss > 1
%         set(gca, 'YTickLabel', []);
%     end
%     set(gca,'XTickLabel', '');
%     
%     %% Per monk
%     for mm = 1:2
%         
%         %% Color plotz -- real
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}{mm}, sites{ss,1});
%         axes(axs(ns*(1+mm)+ss)); cla reset; hold on;
%         vals  = udat{ss,mm,1}(:,:,1);
%         [Y,I] = sort(max(vals,[],2)-min(vals,[],2));
%         vals  = vals(I,:);
%         nv    = size(vals,1);
%         imagesc(xax,-nv:-1,vals,[-0.0004 0.0004]);
%         hold on;
%         plot([0 0], [-nv 1], 'k--', 'LineWidth', 1);
%         axis([xax(1) xax(end) -maxt -1]);
%         if ss > 1
%             set(gca, 'YTickLabel', []);
%         end
%         if mm == 1
%             set(gca, 'XTickLabel', []);
%         end
%         
%         
%         
%         
%         
%         
% %         if ss==ns && mm==2
% %             colorbar
% %         end
%         
%         %% BOOTSTRAP PEAKS
%         Lg = isfinite(udat{ss,mm,1}(:,1,1));
%         for ii = 1:N
%             
%             vals = mean(normrnd(udat{ss,mm,1}(Lg,:,1), ...
%                 udat{ss,mm,1}(Lg,:,2)));
%             
%             maxs(ii,ss,mm) = xax(find(vals==max(vals(Lx)),1));
%             mins(ii,ss,mm) = xax(find(vals==min(vals(Lx)),1));
%         end
%         
%         disp(sprintf('ss=%d, monk=%d, MAX med=%.2f [%.2f %.2f]', ...
%             ss, mm, median(maxs(:,ss,mm)), ...
%             prctile(maxs(:,ss,mm),25), ...
%             prctile(maxs(:,ss,mm),75)))
%         
%         disp(sprintf('ss=%d, monk=%d, MIN med=%.2f [%.2f %.2f]', ...
%             ss, mm, median(mins(:,ss,mm)), ...
%             prctile(mins(:,ss,mm),25), ...
%             prctile(mins(:,ss,mm),75)))
%         
%         disp(' ')
%         
%         plot(median(maxs(:,ss,mm)), -1, sy, ...
%             'MarkerEdgeColor', 'k', 'MarkerFaceColor', co{1});
%         plot(median(mins(:,ss,mm)), -1, sy, ...
%             'MarkerEdgeColor', 0.7.*ones(1,3), 'MarkerFaceColor', co{2});
%         
%     end
% end
% 
% % panel labels
% axs_pl = 1:ss:length(axs)-ns;
% for xx = 1:length(axs_pl)
%     setPLOT_panelLabel(axs(axs_pl(xx)), xx);
% end

%
% %% Stats on peaks -- compute 95% CI from Monte-Carlo
% N    = 1000;
% Lx   = xax<1500;
% maxs = nans(N,ns,2);
% mins = nans(N,ns,2);
% for ss = 1:ns
%
%     %% Per monk
%     monks = find(strcmp(bas{ss}, sites(:,1)));
%     for mm = 1:2
%         Lg = isfinite(udat{monks(mm)}(:,1,1));
%         for ii = 1:N
%
%             vals = mean(normrnd(udat{monks(mm)}(Lg,:,1), ...
%                 udat{monks(mm)}(Lg,:,2)));
%
%             maxs(ii,ss,mm) = xax(find(vals==max(vals(Lx)),1));
%             mins(ii,ss,mm) = xax(find(vals==min(vals(Lx)),1));
%             %subplot(ns,2,(ss-1)*2+mm); cla reset; hold on;
%             %hist(maxs(:,ss,mm),50);
%         end
%     end
% end
%
% %% report bootstrapped median & 95% CI -- MAX
% for ss = 1:ns
%     monks = find(strcmp(bas{ss}, sites(:,1)));
%     for mm = 1:2
%         disp(sprintf('ss=%d, monk=%d, MAX med=%.2f [%.2f %.2f]', ...
%             ss, mm, median(maxs(:,ss,mm)), ...
%             prctile(maxs(:,ss,mm),25), ...
%             prctile(maxs(:,ss,mm),75)))
%     end
% end
% disp(' ')
%
% %% report bootstrapped median & 95% CI -- MAX
% for ss = 1:ns
%     monks = find(strcmp(bas{ss}, sites(:,1)));
%     for mm = 1:2
%         disp(sprintf('ss=%d, monk=%d, MIN med=%.2f [%.2f %.2f]', ...
%             ss, mm, median(mins(:,ss,mm)), ...
%             prctile(mins(:,ss,mm),25), ...
%             prctile(mins(:,ss,mm),75)))
%     end
% end


% %% Summary across sites
% Lxax  = xax>-500&xax<1000;
% axes(axs(end)); cla reset; hold on;
% plot(xax([1 end]), [0 0], 'k:');
% plot([0 0], [-4 4], 'k:');
% for ss = 1:size(sites,1)
%
%     if ~strcmp(sites{ss,1}, 'MT')
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}, sites{ss,1});
%         vals = udat{ss}(:,Lxax,1)-udat{ss}(:,Lxax,3);
%         Lp   = sum(udat{ss}(:,Lxax,4)<0.05,2)>50;
%         xs   = xax(Lxax)';
%         ys   = nanmean(vals(Lp,:));
%         plot(xs, ys, '-', 'Color', co{2});
%         maxi = find(ys==max(ys),1);
%         plot(xs(maxi), ys(maxi), sy, 'Color', co{1}, ...
%             'MarkerFaceColor', co{1}, 'MarkerSize', 10);
%     end
% end
% axis([-500 1000 -0.00025 0.00045]);
% xlabel('Time re: spike (ms)')
% ylabel('Change in PD (z/ms)')

%% Junk Below
%
%         % monkey-by-monkey peaks
%         for ff = 1:length(Fb)
%             % color and symbol0
%             vals = nanmedian(udat{Fb(ff)}(:,:,1)-udat{Fb(ff)}(:,:,3));
%             vals = nanrunmean(vals(xax>0&xax<1200),5);
%             mnv  = min(vals);
%
%
%     % per monk peaks
%     for ff = 1:size(sites{ss,2},2)
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(ff), sites{ss,1});
%         Lm   = udat{ss,3} == ff;
%         vals = udat{ss,1}(Lm,Ltax,1) - udat{ss,1}(Lm,Ltax,3);
%         Lp   = sum(udat{ss,2}(Lm,Ltax,1)<pv,2)>25;
%
% %% Summary (combined): peak values/times
% axes(axs(end)); cla reset; hold on;
% plot(xax([1 end]), [0 0], 'k:');
% plot([0 0], mmx, 'k:')
% for ss = 1:length(bas)
%
%     if ~strcmp(sites{ss,1}, 'MT')
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,1}, sites{ss,2});
%
%     % monkey-by-monkey peaks
%     for ff = 1:length(Fb)
%         % color and symbol0
%         vals = nanmedian(udat{Fb(ff)}(:,:,1)-udat{Fb(ff)}(:,:,3));
%         vals = nanrunmean(vals(xax>0&xax<1200),5);
%         mnv  = min(vals);
%         mnvi = find(vals==mnv,1);
%         mxv  = max(vals(1:mnvi));
%         mxvi = find(vals==mxv,1);
%         plot([mnvi mxvi], [mnv mxv], [sy '-'], 'MarkerSize', 10, ...
%             'Color', co{1}, 'MarkerFaceColor', co{1});
% %         if signtest(udat{Fb(ff)}(:,mxvi,1), udat{Fb(ff)}(:,mxvi,3))<pv
% %             plot(mxvi, mxv, sy, 'MarkerSize', 10, ...
% %                 'Color', co{1}, 'MarkerFaceColor', co{1});
% %         end
%         plot(mnvi, mnv, sy, 'MarkerSize', 10, ...
%             'Color', co{2}, 'MarkerFaceColor', co{2});
% %         if signtest(udat{Fb(ff)}(:,mnvi,1), udat{Fb(ff)}(:,mnvi,3))<pv
% %             plot(mnvi, mnv, sy, 'MarkerSize', 10, ...
% %                 'Color', co{2}, 'MarkerFaceColor', co{2});
% %         end
%     end
% end
% axis([0 1200 -0.00035 0.00045]);
%
%


%% JUNK BELOW
%
%         plot(xax([1 end]), [0 0], 'k:');
%     plot([0 0], mmx, 'k:')
%
%     Fb = find(strcmp(bas{ss}, sites(:,1)));
%     if length(Fb)==1
%         sdat = udat{Fb};
%     else
%         sdat = cat(1, udat{Fb(1)}, udat{Fb(2)});
%     end
%     mn = nanmedian(sdat(:,:,1)-sdat(:,:,3))';
%     se = nanse(sdat(:,:,2))';
%     h  = patch([xax; flipdim(xax,1)], [mn-se; flipdim(mn+se,1)], ...
%         0.6.*ones(1,3));
%     set(h, 'LineStyle', 'none');
%     plot(xax, mn, 'k-');%, 'Color', sites{ss,4});

% session-by-session peaks
%     for ff = 1:length(Fb)
%         vals = udat{Fb(ff)}(:,:,1)-udat{Fb(ff)}(:,:,3);
%         for ee = 1:size(vals,1)
%             mx = max(vals(ee,:));
%             plot(xax(find(vals(ee,:)==mx,1)), mmx(2)*.9, '.', ...
%                 'Color', co{ff}, 'MarkerSize', 3);
%             mn = min(vals(ee,:));
%             plot(xax(find(vals(ee,:)==mn,1)), mmx(1)*.9, 'x', ...
%                 'Color', co{ff}, 'MarkerSize', 3);
%         end
%     end

% pop sig
%     for tt = 1:size(udat{ss},2)
%         Lg = isfinite(udat{ss}(:,tt,1)) & isfinite(udat{ss}(:,tt,3));
%         if signtest(udat{ss}(Lg,tt,1), udat{ss}(Lg,tt,3)) < 0.05
%             plot(xax(tt), mn(tt), 'r.')
%         end
%     end

%end


% ss = 1;
% for tt = 1:size(udat{ss},1)
%
%     %% Examples:[1 21; 3 11; 5 6; 8 5; 10 13];
%     axes(axs(1)); cla reset; hold on;
%     plot(xax([1 end]), [0 0], 'k:');
%     plot([0 0], emx, 'k:')
%     plot(xax, udat{ss}(tt,:,1), 'k-');
%     plot(xax, udat{ss}(tt,:,3), '-', ...
%         'Color', 0.8*ones(1,3));
%     Lp = udat{ss}(tt,:,4)<pv;
%     if any(Lp)
%         plot(xax(Lp), emx(2), 'k.');
%     end
%     axis([xax(1) xax(end) emx]);
% r = input('next')
% end

%% Seperate per monk
% xax = (-1000:2000)';
% mmx = [-0.0003 0.0004];
% for ss = 1:num_sites
%
%     % mean±sem
%     subplot(2, num_sites, ss); cla reset; hold on;
%     plot(xax([1 end]), [0 0], 'k:');
%     plot([0 0], mmx, 'k:')
%
%     mn = nanmedian(udat{ss}(:,:,1)-udat{ss}(:,:,3))';
%     se = nanse(udat{ss}(:,:,2))';
%     h  = patch([xax; flipdim(xax,1)], [mn-se; flipdim(mn+se,1)], ...
%         0.6.*ones(1,3));
%     set(h, 'LineStyle', 'none');
%     plot(xax, mn, 'k-');%, 'Color', sites{ss,4});
% %     for tt = 1:size(udat{ss},2)
% %         Lg = isfinite(udat{ss}(:,tt,1)) & isfinite(udat{ss}(:,tt,3));
% %         if signtest(udat{ss}(Lg,tt,1), udat{ss}(Lg,tt,3)) < 0.05
% %             plot(xax(tt), mn(tt), 'r.')
% %         end
% %     end
%     axis([xax(1) xax(end) mmx]);
%
%     % fraction sig
%     subplot(2, num_sites, num_sites+ss); cla reset; hold on;
%     plot([0 0], [0 1], 'k:')
%     plot(xax, sum(udat{ss}(:,:,4)<0.05&udat{ss}(:,:,1)<0)./ ...
%         sum(isfinite(udat{ss}(:,:,4))), 'k--');
%     plot(xax, sum(udat{ss}(:,:,4)<0.05&udat{ss}(:,:,1)>0)./ ...
%         sum(isfinite(udat{ss}(:,:,4))), 'k-');
%
%     axis([xax(1) xax(end) 0 1]);
% end


%
% %% Scatters of magnitude/timing of max vals
% xax  = (-500:1000)';
% for aa = 1:num_sites
%     subplot(1, num_sites,aa); cla reset; hold on;
%     plot(xax([1 end]), [0 0], 'k:');
%     plot([0 0], [-1 1], 'k:');
% end
% %ainds = -500:500;
% pthresh = 0.05;
% sy = {'ks' 'ko'};
% for ss = 1:num_sites
%
%     % plot mean±sem Spike-triggered PD
%     subplot(1,num_sites,ss);
%     ntr   = size(udat{ss,1},1);
%     spdat = nans(ntr, 4, 2); % 2:time/mean/sem/p, 3:min/max
%     for tt = 1:size(udat{ss,1},1)
%
%         for mm = 1:2
%             pt = udat{ss,2}(tt,1,mm);
%             if isfinite(pt)
%                 spdat(tt,:,mm) = [ ...
%                     pt, ...
%                     udat{ss,1}(tt,pt-xax(1)+1,1)-udat{ss,1}(tt,pt-xax(1)+1,3), ... % mean at peak
%                     udat{ss,1}(tt,pt-xax(1)+1,2), ... % sem at peak
%                     udat{ss,2}(tt,3,mm)];
%             end
%         end
%     end
%
%     plot(squeeze(spdat(:,1,1:2))', squeeze(spdat(:,2,1:2))', '-', 'Color', 0.5*ones(1,3));
%     for mm = 1:2
%         plot(repmat(spdat(:,1,mm)', 2, 1), ...
%             [spdat(:,2,mm)-spdat(:,3,mm) spdat(:,2,mm)+spdat(:,3,mm)]', '-', 'Color', 0.8.*ones(1,3));
%     end
%     for mm = 1:2
%         plot(spdat(:,1,mm), spdat(:,2,mm), sy{mm});
%         Lp = spdat(:,4,mm)<pthresh;
%         if any(Lp)
%             plot(spdat(Lp,1,mm), spdat(Lp,2,mm), sy{mm}, 'MarkerFaceColor', 'k');
%         end
%         axis([xax(1) xax(end) -0.002 0.0035])
%     end
% end

%     % remove shuffled means
%     vals = udat{ss,1}(:,:,1) - udat{ss,1}(:,:,3);
%     nm   = nanmean(vals)';
%     ses  = nanse(vals)';
%     h=patch([xax; flipdim(xax,1)], [nm-ses; flipdim(nm+ses,1)], ...
%         0.8.*ones(1,3));
%     set(h, 'LineStyle', 'none');
%     plot(xax, nm, 'k-');%, 'Color', sites{ss,4});
%
%     disp(sprintf('%s %s (%d): peak = %d', sites{ss,2}, sites{ss,1}, ...
%         sum(isfinite(vals(:,700))), xax(nm==nanmax(nm),1)))
%     for xx = 1:length(xax)
%         Lg = isfinite(vals(:,xx));
%         if signtest(vals(Lg,xx))<pthresh
%             plot(xax(xx), nm(xx), 'k.');%, 'Color', sites{ss,4});
%         end
%     end
%
%     % plot peaks
%     yv = 0.0005;
%     %     if sites{ss,4}=='k'
%     %         yv = 0.0005;
%     %     else
%     %         yv = 0.00055;
%     %     end
%     plot(udat{ss,2}(:,1), yv, 'ko', 'MarkerSize', 4);
%     % , 'MarkerEdgeColor', sites{ss,4}, ...
%
%     Lp = udat{ss,2}(:,2)<(pthresh/length(xax));
%     if any(Lp)
%         plot(udat{ss,2}(Lp,1), yv, 'ko', ... %'MarkerEdgeColor', sites{ss,4}, ...
%             'MarkerFaceColor', 'k', 'MarkerSize', 4);
%         pcts = prctile(udat{ss,2}(Lp,1), [25 50 90]);
%         plot(pcts([1 3]), [yv yv], 'k-');%, 'Color', sites{ss,4});
%         plot(pcts(2), yv, 'sk', ... %'MarkerEdgeColor', sites{ss,4}, ...
%             'MarkerFaceColor', 'k', 'MarkerSize', 6);
%         disp(sprintf('median (%d/%d=%.2f pct)= %.2f, iqr = [%.2f %.2f]', ...
%             sum(Lp), length(Lp), sum(Lp)./length(Lp).*100, ...
%             pcts(2), pcts(1), pcts(3)))
%     end
%     axis([xax(1) xax(end) -0.0003 0.0006]);
%
%     % plot aligned peaks
%     subplot(6,4,(sites{ss,3}-1)*4+(sites{ss,4}-1)*2+2);
%     vdat = nans(size(vals,1), size(ainds,2));
%     for vv = 1:size(vals,1)
%         ai  = ainds + udat{ss,2}(vv,1);
%         Lai = ai>=xax(1) & ai<=xax(end);
%         vdat(vv,Lai) = vals(vv,abs(xax(1))+1+ai(Lai));
%     end
%     nm   = nanmean(vdat(~Lp,:))';
%     ses  = nanse(vdat(~Lp,:))';
%     h=patch([ainds'; flipdim(ainds',1)], [nm-ses; flipdim(nm+ses,1)], ...
%         0.8.*ones(1,3));
%     set(h, 'LineStyle', 'none');
%     plot(ainds', nm, 'k--');%, 'Color', sites{ss,4});
%
%     if any(Lp)
%         nm   = nanmean(vdat(Lp,:))';
%         ses  = nanse(vdat(Lp,:))';
%         h=patch([ainds'; flipdim(ainds',1)], [nm-ses; flipdim(nm+ses,1)], ...
%             0.8.*ones(1,3));
%         set(h, 'LineStyle', 'none');
%         plot(ainds', nm, 'k-');%, 'Color', sites{ss,4});
%     end
%
% %     plot(ainds, vdat(~Lp,:), '--', 'Color', sites{ss,4})
% %     if any(Lp)
% %         plot(ainds, vdat( Lp,:), '-', 'Color', sites{ss,4})
% %     end
%     axis([-500 500 -.001 .002]);
% end


%
% %% raw/aligned to peak
% xax  = (-500:1000)';
% for aa = 1:24
%     subplot(6,4,aa); cla reset; hold on;
%     plot(xax([1 end]), [0 0], 'k:');
% end
% ainds = -500:500;
% pthresh = 0.05;
% for ss = 1:num_sites
%
%     % plot mean±sem Spike-triggered PD
%     subplot(6,4,(sites{ss,3}-1)*4+(sites{ss,4}-1)*2+1);
%
%     % remove shuffled means
%     vals = udat{ss,1}(:,:,1) - udat{ss,1}(:,:,3);
%     nm   = nanmean(vals)';
%     ses  = nanse(vals)';
%     h=patch([xax; flipdim(xax,1)], [nm-ses; flipdim(nm+ses,1)], ...
%         0.8.*ones(1,3));
%     set(h, 'LineStyle', 'none');
%     plot(xax, nm, 'k-');%, 'Color', sites{ss,4});
%
%     disp(sprintf('%s %s (%d): peak = %d', sites{ss,2}, sites{ss,1}, ...
%         sum(isfinite(vals(:,700))), xax(nm==nanmax(nm),1)))
%     for xx = 1:length(xax)
%         Lg = isfinite(vals(:,xx));
%         if signtest(vals(Lg,xx))<pthresh
%             plot(xax(xx), nm(xx), 'k.');%, 'Color', sites{ss,4});
%         end
%     end
%
%     % plot peaks
%     yv = 0.0005;
%     %     if sites{ss,4}=='k'
%     %         yv = 0.0005;
%     %     else
%     %         yv = 0.00055;
%     %     end
%     plot(udat{ss,2}(:,1), yv, 'ko', 'MarkerSize', 4);
%     % , 'MarkerEdgeColor', sites{ss,4}, ...
%
%     Lp = udat{ss,2}(:,2)<(pthresh/length(xax));
%     if any(Lp)
%         plot(udat{ss,2}(Lp,1), yv, 'ko', ... %'MarkerEdgeColor', sites{ss,4}, ...
%             'MarkerFaceColor', 'k', 'MarkerSize', 4);
%         pcts = prctile(udat{ss,2}(Lp,1), [25 50 90]);
%         plot(pcts([1 3]), [yv yv], 'k-');%, 'Color', sites{ss,4});
%         plot(pcts(2), yv, 'sk', ... %'MarkerEdgeColor', sites{ss,4}, ...
%             'MarkerFaceColor', 'k', 'MarkerSize', 6);
%         disp(sprintf('median (%d/%d=%.2f pct)= %.2f, iqr = [%.2f %.2f]', ...
%             sum(Lp), length(Lp), sum(Lp)./length(Lp).*100, ...
%             pcts(2), pcts(1), pcts(3)))
%     end
%     axis([xax(1) xax(end) -0.0003 0.0006]);
%
%     % plot aligned peaks
%     subplot(6,4,(sites{ss,3}-1)*4+(sites{ss,4}-1)*2+2);
%     vdat = nans(size(vals,1), size(ainds,2));
%     for vv = 1:size(vals,1)
%         ai  = ainds + udat{ss,2}(vv,1);
%         Lai = ai>=xax(1) & ai<=xax(end);
%         vdat(vv,Lai) = vals(vv,abs(xax(1))+1+ai(Lai));
%     end
%     nm   = nanmean(vdat(~Lp,:))';
%     ses  = nanse(vdat(~Lp,:))';
%     h=patch([ainds'; flipdim(ainds',1)], [nm-ses; flipdim(nm+ses,1)], ...
%         0.8.*ones(1,3));
%     set(h, 'LineStyle', 'none');
%     plot(ainds', nm, 'k--');%, 'Color', sites{ss,4});
%
%     if any(Lp)
%         nm   = nanmean(vdat(Lp,:))';
%         ses  = nanse(vdat(Lp,:))';
%         h=patch([ainds'; flipdim(ainds',1)], [nm-ses; flipdim(nm+ses,1)], ...
%             0.8.*ones(1,3));
%         set(h, 'LineStyle', 'none');
%         plot(ainds', nm, 'k-');%, 'Color', sites{ss,4});
%     end
%
% %     plot(ainds, vdat(~Lp,:), '--', 'Color', sites{ss,4})
% %     if any(Lp)
% %         plot(ainds, vdat( Lp,:), '-', 'Color', sites{ss,4})
% %     end
%     axis([-500 500 -.001 .002]);
% end

%
%
% % udat is, for each monkey and site:
% %   sessions x time bins
% WIN_START = -2000;
% WIN_END   = 2000;
% xax       = (WIN_START:WIN_END)';
% %co        = {'k' 0.8.*ones(1,3)}; % colors for Oz, Cicero
% if strcmp(ptype, 'slope')
%     lms = [-0.0005 0.0005];
% else
%     lms = [-0.05 0.15];
% end
% for xx = 1:5
%     subplot(5,1,xx); cla reset; hold on;
% end
% for ss = 1:num_sites
%     subplot(5,1,sites{ss,3});
%     plot(xax([1 end]), [0 0], 'k:');
%     vals = udat{ss}(:,:,1) - udat{ss}(:,:,3);
%     vals = vals(isfinite(sum(vals,2)),:);
%     vals = zscore(vals')';
%     mns  = nanmean(vals)';
% %    ses  = nanse(vals)';
% %    h=patch([xax; flipdim(xax,1)], [mns-ses; flipdim(mns+ses,1)], ...
% %        0.8.*ones(1,3));
% %    set(h, 'LineStyle', 'none');
%     plot(xax, mns, '-', 'Color', sites{ss,4})
%     for xx = 1:length(xax)
%         Lg = isfinite(vals(:,xx));
%         if signtest(vals(Lg,xx))<0.05
%             plot(xax(xx), mns(xx), '.', 'Color', sites{ss,4});
%         end
%     end
%     disp(sprintf('%s, %s, peak = %d', sites{ss,2}, sites{ss,1}, ...
%         xax(find(xax>0&mns==max(mns),1))))
%     %axis([WIN_START WIN_END lms])
% end
%
% % color plots
% for ss = 1:num_sites
%     for mm = 1:num_monkeys
%         subplot(2,3,(mm-1)*3+ss); cla reset; hold on;
%         len=size(udat{mm,ss},1);
%         imagesc(xax, 1:len, udat{mm,ss}(:,:,1)-udat{mm,ss}(:,:,5),[-.001 .001]);
%         axis([xax(1) xax(end) 1 len]);
%     end
% end


% %LC 11, 37
% %IC 2
% % look for single sessions
% ii = [1 5];
% co = {'k' 'r'};
% for uu = 1:66
%     for ss = 1:num_sites
%         for mm = 1:num_monkeys
%             if uu<=size(udat{mm,ss},1)
%                 subplot(2,3,(mm-1)*3+ss); cla reset; hold on;
%                 plot(xax([1 end]), [0 0], 'k:');
%                 for jj = 1:2
%                     mns = udat{mm,ss}(uu,:,ii(jj))';
%                     ses = udat{mm,ss}(uu,:,ii(jj)+1)';
%                     h=patch([xax; flipdim(xax,1)], [mns-ses; flipdim(mns+ses,1)], 0.8.*ones(1,3));%co{mm});
%                     set(h, 'LineStyle', 'none');
%                     plot(xax, mns, 'Color', co{jj});
%                 end
%                 axis([-2000 2000 -.2 .2])
%             end
%         end
%     end
%     r = input(sprintf('%d', uu))
% end


%% JUNK BELOW
%     end
%     %axis([-1000 1200 -.0004 .0008])
%
%     axes(axs(ss+3)); cla reset; hold on;
%     for mm = 1:num_monkeys
%         ps = udat{mm,ss}(:,:,4);
%         Lp = (udat{mm,ss}(:,:,1)-udat{mm,ss}(:,:,5))>0;
%         % positive
%         pnum = sum(ps<0.005&Lp);
%         psum = sum(isfinite(ps)&Lp);
%         plot(xax, (pnum./psum)', '-', 'Color', co{mm});
%         % negative
%         pnum = sum(ps<0.005&~Lp);
%         psum = sum(isfinite(ps)&~Lp);
%         plot(xax, -(pnum./psum)', '--', 'Color', co{mm});
%     end
%     axis([xax(1) xax(end) -1 1])
% end
%
%
% for ss = 1:num_sites
%     axes(axs(ss)); cla reset; hold on;
%     plot(xax([1 end]), [0 0], 'k:');
%     for mm = 1:num_monkeys
%         mns  = nanmean(udat{mm,ss}(:,:,1)-udat{mm,ss}(:,:,5))';% - nanmean(udat{mm,ss}(:,:,5))';
%         ses = nanse(udat{mm,ss}(:,:,2))';
%         h=patch([xax; flipdim(xax,1)], [mns-ses; flipdim(mns+ses,1)], co{mm});
%         set(h, 'LineStyle', 'none');
%         plot(xax, mns, '-', 'Color', co{mm})
%     end
%     axis([-1500 1500 -.05 .15])
%     %axis([-1000 1200 -.0004 .0008])
%
%     axes(axs(ss+3)); cla reset; hold on;
%     for mm = 1:num_monkeys
%         ps = udat{mm,ss}(:,:,4);
%         Lp = (udat{mm,ss}(:,:,1)-udat{mm,ss}(:,:,5))>0;
%         % positive
%         pnum = sum(ps<0.005&Lp);
%         psum = sum(isfinite(ps)&Lp);
%         plot(xax, (pnum./psum)', '-', 'Color', co{mm});
%         % negative
%         pnum = sum(ps<0.005&~Lp);
%         psum = sum(isfinite(ps)&~Lp);
%         plot(xax, -(pnum./psum)', '--', 'Color', co{mm});
%     end
%     axis([xax(1) xax(end) -1 1])
% end
%
%     % color plots
%     for mm = 1:num_monkeys
%         axes(axs(ss+3+(mm-1)*3)); cla reset; hold on;
%         len=size(udat{mm,ss},1);
%         imagesc(xax, 1:len, udat{mm,ss}(:,:,1)-udat{mm,ss}(:,:,5),[-.4 .4]);
%         axis([xax(1) xax(end) 1 len]);
%     end
% end
%
%
%     % min/max
%     axes(axs(ss+3)); cla reset; hold on;
%     plot(xax([1 end]), [0 0], 'k:');
%     for mm = 1:num_monkeys
%         for ee = 1:size(udat{mm,ss},1)
%             vals = nanrunmean(udat{mm,ss}(ee,:,1),20);
%             minv = find(vals==min(vals),1);
%             maxv = find(vals==max(vals),1);
%             plot(xax(minv), max(-0.5, vals(minv)), 's', 'MarkerSize', 5, 'Color', co{mm});
%             plot(xax(maxv), min(0.5, vals(maxv)), 'o', 'MarkerSize', 5, 'Color', co{mm});
%         end
%     end
% end
%
%         mns  = nanmean(udat{mm,ss}(:,:,1))';
%
%         ps = udat{mm,ss}(:,:,4);
%         Lp = udat{mm,ss}(:,:,1)>0;
%         % positive
%         pnum = sum(ps<0.01&Lp);
%         psum = sum(isfinite(ps)&Lp);
%         plot(xax, (pnum./psum)', '-', 'Color', co{mm});
%         % negative
%         pnum = sum(ps<0.01&~Lp);
%         psum = sum(isfinite(ps)&~Lp);
%         plot(xax, -(pnum./psum)', '--', 'Color', co{mm});
%     end
%     axis([xax(1) xax(end) -1 1])
% end
%
%
