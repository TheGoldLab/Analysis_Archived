function figLCP_PETH(num, collectData)
% function figLCP_PETH(num, collectData)
%
% uses just constrict/dilate, nothin' else
%
% collectData is flag, indicating whether to re-create data matrix
%   (otherwise load from file)
%
%  siteData{1}: trialsxcols matrix, cols are:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correct
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
    num = 4;
end

if nargin < 2 || isempty(collectData)
    collectData = false;
end

%% Set up figure
wid     = 17.6; % total width
cols    = {5,5,5,5,5};
hts     = [2 2.5 2.5 2.5 2.5];
[axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Joshi et al', true);
set(axs,'Units','normalized');

%% Stuff we need
% site, monkeys, example_data (monkey file ylims)
sites = {...
    'LC',  {'Oz' 'Cicero'},      [1 13  1 10]; ...
    'IC',  {'Oz' 'Cicero'},      [1 11 10 30]; ...
    'SC',  {'Oz' 'Cicero'},      [1  3  6 15]; ...
    'ACC', {'Sprout' 'Atticus'}, [2  3  3  7]; ...
    'PCC', {'Sprout' 'Cheetah'}, [1  12 4 7]};
num_sites = size(sites, 1);

% time axis for PETH
tmin   = -2000;
tmax   = 1000;
tsize  = 250;
tstep  = 10;
tbins  = [(tmin:tstep:tmax-tsize)' (tmin+tsize:tstep:tmax)'];
ntbins = size(tbins,1);
tax    = mean(tbins,2);
etype  = 8; % 2=peak, 8=max slope
tdat   = nans(1,ntbins,4); % PETH data, two conditions by pupil (see below)

%% Possibly get data
if collectData
    
    % 1:PETH, 2:corr, 3:[monk file unit], 4:pd stats, 5:example raster
    % data, 6: example raster selection array
    udat = cell(num_sites, 6); 
    
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
                    num_sunits = size(siteData{3},2)-1;
                    if num_sunits > 0
                        
                        % get spike data using only no-beep trials
                        Levents = ismember(siteData{5}(:,1), ...
                            find(~isfinite(siteData{1}(:,4)) & ....
                            ~isfinite(siteData{1}(:,9))));
                        Fevents = find(Levents);
                        num_events = length(Fevents);
                        sdat = nans(num_events, ntbins, num_sunits);
                        
                        % get PD stats
                        Fpds = find(Levents(1:end-1) & ...
                            siteData{5}(1:end-1,3) == siteData{5}(2:end,2) & ...
                            siteData{5}(1:end-1,9) > 0);
                        pd_stats = prctile(siteData{5}(Fpds+1,3)-siteData{5}(Fpds,2), ...
                            [0 5 25 50 75 95 100]);
                        
                        % loop through the pupil events
                        for ee = 1:num_events
                            
                            % index in event matrix (row of siteData{5})
                            ev = Fevents(ee);
                            
                            % spikes wrt pupil events
                            etax = siteData{5}(ev,etype)+tbins;
                            si   = find(etax(:,1)>=siteData{1}(siteData{5}(ev,1),5), 1);
                            ei   = find(etax(:,2)<=siteData{1}(siteData{5}(ev,1),6), 1, 'last');
                            if ei > si
                                for uu = 1:num_sunits
                                    sp = siteData{3}{siteData{5}(ev,1), uu+1};
                                    % sp = siteData{3}{siteData{5}(ev,1), 1};
                                    for tt = si:ei
                                        sdat(ee,tt,uu) = sum( ...
                                            sp>=etax(tt,1)&sp<=etax(tt,2));
                                    end
                                end
                            end
                        end
                        
                        % get PETHs, stats
                        
                        % PD magnitude measure -- P-P or max slope
                        if etype == 2
                            % difference between peak and trough
                            Zmag = diff(siteData{5}(Levents,6:7),[],2);
                        else
                            % max slope
                            Zmag = siteData{5}(Levents,9);
                        end
                        
                        % PD at start of event
                        Zi = siteData{5}(Levents,4);
                        
                        % get residuals of Zmag, after linear regression with Z
                        Lr    = [Zmag>0 Zmag<0];
                        A     = [Zi ones(size(Zi))];
                        ZmagR = nans(size(Zmag));
                        for rr = 1:2
                            X = lscov(A(Lr(:,rr),:), Zmag(Lr(:,rr)));
                            ZmagR(Lr(:,rr)) = Zmag(Lr(:,rr))-A(Lr(:,rr),:)*X;
                        end
                        
                        % select dilation/constriction events
                        %   1: big dilation
                        %   2: big constriction
                        Lp = [Lr(:,1)&ZmagR>=0 Lr(:,2)&ZmagR<=0];
                        nlp  = size(Lp,2);
                        
                        for uu = 1:num_sunits
                            
                            % get peths (& sem)
                            tdat(:) = nan;
                            for pp = 1:nlp
                                tdat(1,:,pp) = nanmean(sdat(Lp(:,pp),:,uu))./ ...
                                    tsize.*1000;
                                tdat(1,:,pp+2) = nanse(sdat(Lp(:,pp),:,uu))./ ...
                                    tsize.*1000;
                            end
                            
                            % save peth
                            udat{ss,1} = cat(1, udat{ss,1}, tdat);
                            
                            % compare medians
                            tdat(1,:,1) = nan;
                            for cc = 1:ntbins
                                
                                % non-nan trials
                                Lg = isfinite(sdat(:,cc,uu));
                                if sum(Lg) >= 10
                                    % compare spike-rate distributions for the two PD
                                    % groups..
                                    tdat(1,cc,1) = ranksum(sdat(Lg&Lp(:,1),cc,uu), ...
                                        sdat(Lg&Lp(:,2),cc,uu));
                                 end
                            end
                            
                            % save 'em
                            udat{ss,2} = cat(1, udat{ss,2}, tdat(1,:,1));
                            udat{ss,3} = cat(1, udat{ss,3}, [mm ff uu]);
                            udat{ss,4} = cat(1, udat{ss,4}, pd_stats);
                        end
                    end
                end
            end
        end
    end
    
    % save example rasters
    for ss = 1:num_sites
        [base_dir, fnames] = getLCP_cleanDataDir(...
            sites{ss,2}{sites{ss,3}(1)}, sites{ss,1});
        Fm = find(udat{ss,3}==sites{ss,3}(1));
        ff = udat{ss,3}(Fm(sites{ss,3}(2)),2);
        uu = udat{ss,3}(Fm(sites{ss,3}(2)),3)+1;
        load(fullfile(base_dir, fnames{ff}));
                    
        % get spike data using only no-beep trials
        Levents = ismember(siteData{5}(:,1), ...
            find(~isfinite(siteData{1}(:,4)) & ....
            ~isfinite(siteData{1}(:,9))));
        Fevents = find(Levents);
        num_events = length(Fevents);
        udat{ss,5} = cell(num_events,1);
        
        % loop through the pupil events
        for ee = 1:num_events

            % index in event matrix (row of siteData{5})
            ev = Fevents(ee);
            % spikes wrt pupil events, within tax range
            sp = siteData{3}{siteData{5}(ev,1), uu}-siteData{5}(ev,etype);
            udat{ss,5}{ee} = sp(sp>=tmin&sp<=tmax);
        end
        
        % make the selection array
        % PD magnitude measure -- P-P or max slope
        if etype == 2
            % difference between peak and trough
            Zmag = diff(siteData{5}(Levents,6:7),[],2);
        else
            % max slope
            Zmag = siteData{5}(Levents,9);
        end
        % PD at start of event
        Zi = siteData{5}(Levents,4);
        
        % get residuals of Zmag, after linear regression with Z
        Lr    = [Zmag>0 Zmag<0];
        A     = [Zi ones(size(Zi))];
        ZmagR = nans(size(Zmag));
        for rr = 1:2
            X = lscov(A(Lr(:,rr),:), Zmag(Lr(:,rr)));
            ZmagR(Lr(:,rr)) = Zmag(Lr(:,rr))-A(Lr(:,rr),:)*X;
        end
        
        % Save ZmagR, type (below), sorted indices
        %   1: big dilation
        %   2: big constriction
        Fdc    = (Lr(:,1)&ZmagR>=0)+2.*(Lr(:,2)&ZmagR<=0);
        Fd     = find(Lr(:,1));
        [~,dI] = sort(ZmagR(Fd)); % sorted dilations
        Fc     = find(Lr(:,2));
        [~,cI] = sort(ZmagR(Fc)); % sorted constrictions
        inds   = nans(size(ZmagR));
        inds(Fc) = cI;
        inds(Fd) = dI+max(Fc);
        udat{ss,6} = cat(2, ZmagR, Fdc, inds);
    end
    
    % save data to file
    FS_saveProjectFile('2013_LCPupil', 'PETH', udat);
else
    % load data from file
    udat = FS_loadProjectFile('2013_LCPupil', 'PETH');
end

%% PLOTZ
rinds = [2 1];
emax  = 45;
maxt = 0;
pv = 0.05;
for ss = 1:num_sites
    for mm = 1:2
        if sum(udat{ss,3}(:,1) == mm) > maxt
            maxt = sum(udat{ss,3}(:,1) == mm);
        end
    end
end
for ss = 1:num_sites
    
    %% Example sessions
    [co, ~] = getLCP_colorsAndSymbol(sites{ss,2}{sites{ss,3}(1)}, sites{ss,1});
    
    % raster
    axes(axs(ss)); cla reset; hold on;
    for ii = 1:2
        Fi = find(udat{ss,6}(:,2)==rinds(ii));
        %disp(length(Fi))
        Fi = Fi(randperm(length(Fi),emax));
        for jj = 1:length(Fi)
            xs = udat{ss,5}{Fi(jj)}';
            ys = repmat((ii*2-3)*jj, size(xs));
            plot(xs, ys, '.', 'Color', co{3-ii}, 'MarkerSize', 2);
            %plot([xs; xs], [ys; ys+0.9], '-', 'Color', co{3-ii});
        end
    end
    plot([0 0], [-800 800], 'k:');
    axis([tax(1) tax(end) -emax emax]);
    set(gca,'XTickLabel', '');
    title(sprintf('%s', sites{ss,1}));
    if ss > 1
        set(gca, 'YTickLabel', '');
    end
    
    % peth
    axes(axs(num_sites+ss)); cla reset; hold on;
    plot([0 0], [-30 30], 'k:');
    Fm = find(udat{ss,3}(:,1)==sites{ss,3}(1));
    sd = squeeze(udat{ss,1}(Fm(sites{ss,3}(2)),:,[1 2]));
    plot(tax, sd(:,1), '-', 'Color', co{1});
    plot(tax, sd(:,2), '-', 'Color', co{2});
    %Lp = udat{ss,2}(Fm(sites{ss,3}(2)),:,1)'<pv;
    %plot(tax(Lp), (sites{ss,3}(4)).*ones(sum(Lp),1), 'k.');
    axis([tax([1 end])' sites{ss,3}(3:4)]);
    set(gca,'XTickLabel', '');    

    %% summary mean±sem difference between big dilate/constrict, both monks
    axes(axs(num_sites*2+ss)); cla reset; hold on;
    plot(tax([1 end]), [0 0], 'k:');
    plot([0 0], [-30 30], 'k:');
    vals = udat{ss,1}(:,:,1)-udat{ss,1}(:,:,2);
    mn = nanmean(vals)';
    se = nanse(vals)';
    h=patch([tax; flipud(tax)], [mn-se; flipud(mn+se)], co{2}); %0.8.*ones(1,3));
    set(h, 'LineStyle', 'none');
    plot(tax, nanmean(vals), '-', 'Color', co{1});
    axis([tax(1) tax(end) -1 1.4]);
    set(gca,'XTickLabel', '');
    if ss > 1
        set(gca, 'YTickLabel', '');
    end
    
    % label peak
    peaki = mn==max(mn);
    h=text(-1600, 1, sprintf('%d ms', tax(peaki)));
    set(h,'FontSize',9);
    % bootstrap!
    N   = 1000;
    bsd = normrnd(repmat(vals(:,peaki),1,N), ...
        repmat(cat(1,sqrt(udat{ss,1}(:,peaki,3).^2+udat{ss,1}(:,peaki,4).^2)),1,N));
    if prctile(median(bsd),2.5)>0
        set(h, 'FontWeight', 'bold');
    end

    for mm = 1:size(sites{ss,2},2)
        
        %% per session color plotz of differences
        Lm    = udat{ss,3}(:,1) == mm;
        vals  = udat{ss,1}(Lm,:,1)-udat{ss,1}(Lm,:,2);
        nv    = size(vals,1);
        [~,I] = sort(max(vals,[],2)-min(vals,[],2));
        vals  = vals(I,:);        
        vals  = vals - repmat(nanmean(vals,2),1,size(vals,2));
        
        axes(axs(num_sites*(mm+2)+ss)); cla reset; hold on;
        imagesc(tax,-nv:-1,vals,[-1.5 1.5]);
        hold on;
        plot([0 0], [-maxt -1], 'k:');
        axis([tax(1) tax(end) -maxt -1]);
        if mm==1
            set(gca,'XTickLabel', '');
        end
        if ss > 1
            set(gca, 'YTickLabel', '');
        end

    end
end

% axis labelsaxes(axs(1));
axes(axs(1));
ylabel('Event number')
axes(axs(11));
ylabel('Response (sp/s)')
axes(axs(21));
ylabel('Unit number')
xlabel('Time re:PD event (ms)');

tdat = cell(num_sites,2);
for ss = 1:num_sites
    for mm = 1:2
        Fm    = find(udat{ss,3}(:,1) == mm);
        tdat{ss,mm} = nans(length(Fm),2);
        for uu = 1:length(Fm)
            vals = udat{ss,1}(Fm(uu),:,1)-udat{ss,1}(Fm(uu),:,2);
            Lp = (vals>0 & udat{ss,2}(Fm(uu),:,1)<pv)';
            d = diff([false;Lp;false]);
            p = find(d==1);
            q = min(length(tax),find(d==-1));
            Lgood = tax(q)>-1000&tax(p)<100&(q-p)>=7;
            
            %             cla reset; hold on;
            %             Lp=udat{ss,2}(Fm(uu),:,1)<pv;
            %             plot(tax, vals, 'k-');
            %             plot(tax(Lp), vals(Lp), 'r.');
            
            if any(Lgood)
                p = p(Lgood);
                q = q(Lgood);
                for pp = sum(Lgood):-1:1
                    pvals = smooth(vals(p(pp):q(pp)-1),10,'sgolay');
                    maxv  = max(pvals);
                    maxt  = tax(p(pp)-1+find(pvals==maxv,1));
                    % ix    = pp;
                    if maxt < 0
                        break;
                    end
                end
                
                %                 plot(tax(p(ix):q(ix)-1), pvals, 'b-', 'LineWidth', 2);
                %                 plot(maxt, maxv, 'go', 'MarkerSize', 12);
                
                tdat{ss,mm}(uu,:) = [maxt maxv];
            end
            %             axis([tax(1) tax(end) -6 6]);
            %             r = input('next')
        end
    end
end

% Anova
vals   = [];
g_site = {};
g_monk = {};
c2dat  = nans(num_sites,2,2);
for ss = 1:num_sites
    for mm = 1:2
        axes(axs(num_sites*(2+mm)+ss));
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
       
%         % indicate max
%         mxi = max0 + find(mn(Lmax)==max(mn(Lmax)),1);
%         h=plot(tax(mxi), mn(mxi), 'ko');
%         Lg = isfinite(vals(:,mxi));
%         if signrank(vals(Lg,mxi))<pv
%             set(h,'MarkerFaceColor','k')
%         end
%         text(100, 1, num2str(tax(mxi)))
%         axis([tax(1) tax(end) -1 1.3])
%         
% 
% 
% 
%         if ss==num_sites && mm==2
%             colorbar;
%         end
%         
%         % indicate time of max of mean across sessions
%         meanv = nanmean(vals);
%         semv  = nanse(vals);
%         mcs   = normrnd(repmat(meanv,N,1), repmat(semv,N,1));
%         
%         for ii = 1:N
%             maxs(ii) = tax(find(mcs(ii,:)==max(mcs(ii,:)),1));
%         end
%         
%         %mxv   = max(meanv);
%         %tv    = tax(find(meanv==mxv,1)); 
%         
%         plot(median(maxs(ii)), -1, sy, 'MarkerEdgeColor', 'k', ...
%             'MarkerFaceColor', co{1});
%         disp(sprintf('ss=%s, monk=%s, max=%.2f [%.2f %.2f]', ...
%             sites{ss,1}, sites{ss,2}{mm}, ...
%             median(maxs), prctile(maxs,25), prctile(maxs,75)))
%         if mm==1
%             set(gca,'XTickLabel', '');
%         end
%         if ss > 1
%             set(gca, 'YTickLabel', '');
%         end
% 
%     end
% end
%     
%     
% for ss = 1:num_sites
% 
%     
%     %% Summary, per monk:
%     %   1. counts of sig sites
%     axes(axs(num_sites*2+ss)); cla reset; hold on;
%     plot([0 0], [-2 2], 'k:');
%     % differences between big dilation/constriction
%     vals = udat{ss,1}(:,:,1) - udat{ss,1}(:,:,5); %all: 2-4; big:1-5
%     nses = sum(isfinite(udat{ss,2}(:,:,1)));
%     pyp  = nansum(udat{ss,2}(:,:,2)<pv&vals>0)./nses; %all=1; big=2
%     pyn  = nansum(udat{ss,2}(:,:,2)<pv&vals<0)./nses;
%     Lg   = nses>5;
%     plot(tax(Lg), pyp(Lg), '-', 'Color', co{1}); % nanrunmean(pyp(Lg),1)
%     plot(tax(Lg), pyn(Lg), '-', 'Color', co{2}); % nanrunmean(pyn(Lg),1)
%     axis([tax(1) tax(end) 0 0.5])
%     set(gca,'XTickLabel', '');
%     if ss > 1
%         set(gca, 'YTickLabel', '');
%     end
% end
% 
%     %   2. color map of values
%     for mm = 1:size(sites{ss,2},2)
%         
%         %% per session color plotz of differences
%         Lm    = udat{ss,3}(:,1) == mm;
%         vals  = udat{ss,1}(Lm,:,1)-udat{ss,1}(Lm,:,5); %all: 2-4; big:1-5
%         nv    = size(vals,1);
%         [~,I] = sort(max(vals,[],2)-min(vals,[],2));
%         vals  = vals(I,:);
%         vals  = vals - repmat(nanmean(vals,2),1,size(vals,2));
%         
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}{mm}, sites{ss,1});
%         axes(axs(num_sites*(mm+2)+ss)); cla reset; hold on;
%         imagesc(tax,-nv:-1,vals,[-1.5 1.5]);
%         hold on;
%         plot([0 0], [-nv -1], 'k-', 'LineWidth', 1.5);
%         axis([tax(1) tax(end) -maxt -1]);
%         
%         if ss==num_sites && mm==2
%             colorbar;
%         end
%         
%         % indicate time of max of mean across sessions
%         meanv = nanmean(vals);
%         semv  = nanse(vals);
%         mcs   = normrnd(repmat(meanv,N,1), repmat(semv,N,1));
%         
%         for ii = 1:N
%             maxs(ii) = tax(find(mcs(ii,:)==max(mcs(ii,:)),1));
%         end
%         
%         %mxv   = max(meanv);
%         %tv    = tax(find(meanv==mxv,1)); 
%         
%         plot(median(maxs(ii)), -1, sy, 'MarkerEdgeColor', 'k', ...
%             'MarkerFaceColor', co{1});
%         disp(sprintf('ss=%s, monk=%s, max=%.2f [%.2f %.2f]', ...
%             sites{ss,1}, sites{ss,2}{mm}, ...
%             median(maxs), prctile(maxs,25), prctile(maxs,75)))
%         if mm==1
%             set(gca,'XTickLabel', '');
%         end
%         if ss > 1
%             set(gca, 'YTickLabel', '');
%         end
% 
%     end
% end
% 
% % panel labels
% axs_pl = 1:ss:length(axs)-num_sites;
% for xx = 1:length(axs_pl)
%     setPLOT_panelLabel(axs(axs_pl(xx)), xx);
% end

%% Bootstrap peaks


%% Summary across sites
% Ltax = tax>-1000&tax<500;
% Ftax = find(Ltax);
% axes(axs(end)); cla reset; hold on;
% plot(tax([1 end]), [0 0], 'k:');
% plot([0 0], [-4 4], 'k:');
% for ss = 1:num_sites
%
%     % per monk peaks
%     for ff = 1:size(sites{ss,2},2)
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(ff), sites{ss,1});
%         Lm   = udat{ss,3} == ff;
%         % differences between big dilation/constriction
%         vals = udat{ss,1}(Lm,Ltax,1) - udat{ss,1}(Lm,Ltax,3);
%         Lp   = sum(udat{ss,2}(Lm,Ltax,1)<pv,2)>25;
%         % correlation coefficients
%                 % vals = udat{ss,2}(Lm,Ltax,6);
%                 % Lp   = sum(udat{ss,2}(Lm,Ltax,7)<pv,2)>25;
%         xs   = tax(Ltax)';
%         ys   = nanmean(vals(Lp,:));
%         plot(xs, ys, '-', 'Color', co{1});
%         maxi = find(ys==max(ys),1);
%         plot(xs(maxi), ys(maxi), sy, 'Color', co{1}, ...
%             'MarkerFaceColor', co{1}, 'MarkerSize', 10);
% %         % Stats on partial corrs
% %         pcs = udat{ss,2}(Lm,Ftax(maxi),6);
% %         pcs = pcs(isfinite(pcs));
% %         pcv = udat{ss,2}(Lm,xs<-700,6);
% %         pdat = nans(size(pcv,2),1);
% %         for pp = 1:size(pcv,2)
% %             Lg = isfinite(pcv(:,pp));
% %             pdat(pp) = signrank(pcv(Lg,pp));
% %         end
% %         disp(sprintf('%s, %s: %.2f [%.2f %.2f], p=%.2f, %d/%d', ...
% %             sites{ss,1}, sites{ss,2}{ff}, median(pcs), ...
% %             prctile(pcs,25), prctile(pcs,75), signrank(pcs), ...
% %             sum(pdat<0.05), size(pdat,1)))
%
%     end
% end
% axis([-1000 500 -3 3.5]);
% %axis([-1000 500 -0.1 0.2]);
% xlabel('Time re: spike (ms)')
% ylabel('Difference in response (sp/s)')

%
% %% plot re-centered time-dependent correlation coefficients
% xs = tax(Ltax)./1000;
% for ss = 1:num_sites
%     for mm = 1:size(sites{ss,2},2)
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(mm), sites{ss,1});
%         subplot(2,num_sites,(mm-1)*num_sites+ss); cla reset; hold on;
%         plot(xs([1 end]), [0 0], 'k--');
%         plot([0 0], [-1 1], 'k--');
%         Lm   = udat{ss,3} == mm;
%         vals = udat{ss,2}(Lm,Ltax,6);
%         vals = vals - repmat(nanmean(vals,2),1,size(vals,2));
%         plot(xs, nanmean(vals), '-', 'Color', co{1});
%         axis([-1.0 0.5 -0.07 0.11]);
%     end
% end
%
% %% Stats from smoothed time-dependent correlation coefficients
% fdat = cell(num_sites, 2); % per monk
% xs   = tax(Ltax)./1000;
% for ss = 1:num_sites
%     for mm = 1:size(sites{ss,2},2)
%         Lm   = udat{ss,3} == mm;% & Lp;
%         vals = udat{ss,2}(Lm,Ltax,6);
%         for vv = 1:size(vals,1)
%             svals = smooth(vals(vv,:),'rlowess');
%             cla reset; hold on;
%             plot(xs, vals(vv,:), 'r-');
%             plot(xs, svals, 'k-');
%             ylim([-.2 .2])
%             r = input('next')
%         end
%     end
% end
%
%
% %% Gabor fits to time-dependent correlation coefficients
% fdat = cell(num_sites, 2); % per monk
% xs   = tax(Ltax)./1000;
% for ss = 1:num_sites
%     % per monk
%     for mm = 1:size(sites{ss,2},2)
%         Lm       = udat{ss,3} == mm;% & Lp;
%         vals = udat{ss,2}(Lm,Ltax,6); % correlation coefficients
%         lms  = 1000./udat{ss,4}(Lm,[4 2 6]);
%         fdat{ss,ff} = nans(size(vals,1), 6);
%         for vv = 1:size(vals,1)
%             if ~any(~isfinite(vals(vv,:)))
%                 [fits, stats] = fitLCP_gabor(vals(vv,:)', xs, lms(vv,1), true, true);
%                 title(sprintf('p=%.2f, f=%.2f [%.2f ? %.2f]', ...
%                     stats(2), fits(5), lms(vv,3), lms(vv,2)))
%                 r = input('next')
%             end
%         end
%     end
% end
%

%
% % get stats from fits
% gdat = cell(num_sites, 2); % per site, monk
% xs   = tax(Ltax)./1000;
% Lxs  = [xs<-0.5 xs>=-0.5&xs<0 xs>=0];
% for ss = 1:num_sites
%
%     %Lp  = sum(udat{ss,2}(:,Ltax,7)<pv,2)>25;
%
%     % per monk peaks
%     for ff = 1:size(sites{ss,2},2)
%         Lm       = udat{ss,3} == ff;% & Lp;
%         vals     = udat{ss,2}(Lm,Ltax,6); % correlation coefficients
%         nv       = size(vals,1);
%         gdat{ss,ff} = nans(nv,3,3);
%         for vv = 1:nv
%             ys = fitLCP_gaborVal(fdat{ss,ff}(vv,:), xs);
%             for xx = 1:3
%                 if any(isfinite(ys(Lxs(:,xx))))
%                     ys_x = ys(Lxs(:,xx));
%                     xs_x = xs(Lxs(:,xx));
%                     mxi  = find(ys_x==max(ys_x),1);
%                     mni  = find(ys_x==min(ys_x),1);
%                     gdat{ss,ff}(vv,:,xx) = [mxi mni ys_x(mxi)-ys(mni)];
%                 end
%             end
%         end
%     end
% end
%
% for ss = 1:num_sites
%     for mm = 1:size(gdat,2)
%         subplot(2,num_sites,(mm-1)*num_sites+ss); cla reset; hold on;
%         plot([-0.1 0.9], [-0.1 0.9], 'k--');
%         plot(gdat{ss,mm}(:,3,1), gdat{ss,mm}(:,3,2), 'k.');
%     end
% end
%
%
%
%
%
%                 if abs(xs_x(mxi)-xs_x(mni)) < 0.75/fdat{ss,ff}(vv,5)
%                     gdat{ss,ff}(vv,:) = [mxi mni ys(mxi)-ys(mni)];
%                 end
%             end
%     end
% end
%
% cla reset; hold on;
% for ss = 1:num_sites
%     for mm = 1:2
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(mm), sites{ss,1});
%         Lg = isfinite(gdat{ss,mm}(:,1));
%         plot(mean(tax(gdat{ss,mm}(Lg,1))), mean(gdat{ss,mm}(Lg,3)), ...
%             sy, 'Color', co{1});
%     end
% end
%
%
% % plot em
% cla reset; hold on;
% for ss = 1:num_sites
%     for mm = 1:2
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(mm), sites{ss,1});
%         Lg = isfinite(fdat{ss,mm}(:,1));
%         xs = prctile(fdat{ss,mm}(Lg,3), [25 50 75]);
%         ys = prctile(fdat{ss,mm}(Lg,2), [25 50 75]);
%         plot(xs([1 3]), ys([2 2]), '-', 'Color', co{1});
%         plot(xs([2 2]), ys([1 3]), '-', 'Color', co{1});
%         plot(xs(2), ys(2), sy, 'Color', co{1});
%     end
% end
%
%
% cla reset; hold on;
% for ss = 1:num_sites
%     for mm = 1:2
%         subplot(num_sites,2,(ss-1)*2+mm); cla reset; hold on;
%         plot(fdat{ss,mm}(:,3), fdat{ss,mm}(:,2), 'k.')
%         %        plot(fdat{ss,mm}(:,4), fdat{ss,mm}(:,5), 'k.')
%         axis([-1 .5 0 0.6])
%     end
% end
%
%
% xs = tax(Ltax)./1000;
% for ss = 1:num_sites
%
%     Lp  = sum(udat{ss,2}(:,Ltax,7)<pv,2)>25;
%
%     % per monk peaks
%     for ff = 1:size(sites{ss,2},2)
%         [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}(ff), sites{ss,1});
%         Lm       = udat{ss,3} == ff;% & Lp;
%         vals     = udat{ss,2}(Lm,Ltax,6);
%         Lv       = udat{ss,2}(Lm,Ltax,7)<pv;
%         for vv = 1:size(vals,1)
%             cla reset; hold on;
%             plot([-1000 500], [0 0], 'k:');
%             plot(xs, vals(vv,:), 'k-');
%             plot(xs(Lv(vv,:)), vals(vv,Lv(vv,:)), 'r.')
%             axis([-1000 500 -0.2 0.4])
%
%             if isfinite(fdat{ss,ff}(vv,1))
%                 plot(xs, fitLCP_gaborVal(fdat{ss,ff}(vv,:), xs), 'r-')
%             end
%             axis([-1 .5 -.2 .2])
%             title(sprintf('%.2f  ', fdat{ss,ff}(vv,:)))
%             r = input('next')
%         end
%     end
% end
%
% ss=1
% for vv = 1:size(vals,1)
%     cla reset; hold on;
%     plot([-1000 500], [0 0], 'k:');
%     xs   = tax(Ltax)./1000;
%     plot(xs, vals(vv,:), 'k-');
%     plot(xs(Lv(vv,:)), vals(vv,Lv(vv,:)), 'r.')
%     axis([-1000 500 -0.2 0.4])
%
%     fits = fitLCP_gabor(vals(vv,:)', xs);
%     plot(xs, fitLCP_gaborVal(fits, xs), 'r-')
%     axis([-1 .5 -.2 .2])
%     title(sprintf(
%     r = input('next')
% end
%
%         xs   = tax(Ltax)';
%         ys   = nanmean(vals(Lp,:));
%         plot(xs, ys, '-', 'Color', co{1});
%         maxi = find(ys==max(ys),1);
%         plot(xs(maxi), ys(maxi), sy, 'Color', co{1}, ...
%             'MarkerFaceColor', co{1}, 'MarkerSize', 10);
%     end
% end



%% Junk below
%
%
%         subplot(num_sites, 2, (ss-1)*2+ff); cla reset; hold on;
%         plot([0 0], [-6 6], 'k:');
%         plot([-2000 1000], [0 0], 'k:');
%         plot(tax(Ltax)', vals(Lp,:), 'y--');
%     end
% end
%
%
%         Ln   = udat{ss,2}(Lm,:,1)<pv&vals<0;
%         for ii = 1:size(vals,1)
%             if any(Lp(ii,:))
%                 plot(tax(Lp(ii,:)), vals(ii,Lp(ii,:)), 'y-');
%             end
%             if any(Ln(ii,:))
%                 plot(tax(Ln(ii,:)), vals(ii,Ln(ii,:)), 'y--');
%             end
%         end
%         tmpd = nans(size(vals));
%         tmpd(Lp) = vals(Lp);
%         plot(tax, nanmean(tmpd), 'k-');
%         tmpd = nans(size(vals));
%         tmpd(Ln) = vals(Ln);
%         plot(tax, nanmean(tmpd), '-');
%     end
% end
%
%
%
%
%
%         % vals = nanmean(udat{ss,1}(Lm,:,1) - udat{ss,1}(Lm,:,3))';
%         vals = nanmedian(udat{ss,1}(Lm,:,1) - udat{ss,1}(Lm,:,3))';
%         nm   = nanrunmean(vals,10)';
%
%         plot(tax, udat{ss,1}(Lm,:,1) - udat{ss,1}(Lm,:,3))
%         plot(tax, nm, 'r-', 'LineWidth', 2)
%     end
% end
%
%         nmn  = min(nm(Ltmin));
%         Fmn  = find(nm(Ltmin)==nmn,1)+F0min;
%         if ss==3 && ff==2
%             nmx = max(nm(tax>-300&tax<0));
%             Fmx = find(nm(tax>-300)==nmx,1)+find(tax>=-300,1);
%         else
%             nmx = max(nm(Ltmax));
%             Fmx = find(nm(Ltmax)==nmx,1)+F0max;
%         end
%         plot(tax([Fmx Fmn]), [nmx nmn], '-', 'Color', co{1});
%         plot(tax(Fmx), nmx, sy, 'Color', co{1}, ...
%             'MarkerFaceColor', co{1}, 'MarkerSize', 10);
% %         if signrank(udat{ss,1}(Lm,Fmx,1), udat{ss,1}(Lm,Fmx,3))<pv
% %             plot(tax(Fmx), vals(Fmx), sy, 'Color', co{1}, ...
% %                 'MarkerFaceColor', co{1}, 'MarkerSize', 10);
% %         end
%         plot(tax(Fmn), nmn, sy, 'Color', co{2}, 'MarkerFaceColor', co{2}, ...
%             'MarkerSize', 10);
% %         if signrank(udat{ss,1}(Lm,Fmn,1), udat{ss,1}(Lm,Fmn,3))<pv
% %             plot(tax(Fmn), vals(Fmn), sy, 'Color', co{2}, ...
% %                 'MarkerFaceColor', co{2}, 'MarkerSize', 10);
% %         end
%     end
% end
% axis([-1000 0 -0.7 1.0]);

%% correlation coefficients
%     axes(axs(ss+3*num_sites)); cla reset; hold on;
%     plot(tax([1 end]), [0 0], 'k:');
%     for pp = 1:2
%         meanv = nanmean(udat{ss,2}(:,:,pp*2))';
%         Lg    = isfinite(meanv) & sum(isfinite(udat{ss,2}(:,:,pp*2)))'>5;
%         plot(tax(Lg), meanv(Lg), '-', 'Color', co2{pp});
%         for bb = 1:ntbins
%             Lg = isfinite(udat{ss,2}(:,bb,pp*2));
%             if signtest(udat{ss,2}(Lg,bb,pp*2))<pv
%                 plot(tax(bb), meanv(bb), '.', 'Color', co2{pp});%, ...
% %                    'MarkerSize', 3);
%             end
%         end
%     end
%     axis([xlm -0.05 0.07])
% end

%% JUNK BELOW
% ss=2;
% for vv = 1:127
%     cla reset; hold on;
%     plot([0 0], [-30 30], 'k:');
%     %Fm = find(udat{ss,3}==sites{ss,3}(1));
%     sd = squeeze(udat{ss,1}(vv,:,:));
%     Lp = udat{ss,2}(vv,:,1)'<pv;
%     for ii = 1:3
%         plot(tax, sd(:,ii), '-', 'Color', co{ii});
%     end
%     r = input('next');
% end
%
% %plot(tax(Lp), (sites{ss,3}(4)).*ones(sum(Lp),1), 'k.');
% %axis([tax([1 end])' sites{ss,3}(3:4)]);
%
% Example sessions
% xlm = [-1500 200];
% co  = {'r' [1 .8 .8]; 'g' [.8 1 .8]; 'b' [.8 .8 1]};
% for ss = 1:size(exampleSessions,1)
%     mm = find(strcmp(exampleSessions{ss,1}, monkeys));
%     ss = find(strcmp(exampleSessions{ss,2}, sites));
%     ff = exampleSessions{ss,3};
%     uu = exampleSessions{ss,4};
%
%     [base_dir, fnames] = getLCP_cleanDataDir('Oz', 'LC');
%     load(fullfile(base_dir, fnames{ff}));
%
%     % get spike data using only no-beep trials
%     Levents = ismember(siteData{5}(:,1), ...
%         find(~isfinite(siteData{1}(:,4))));
%     Fevents = find(Levents);
%     num_events = length(Fevents);
%     sdat = nans(num_events, ntbins);
%
%     % loop through the pupil events
%     for ee = 1:num_events
%
%         % index in event matrix (row of siteData{5})
%         ev = Fevents(ee);
%
%         % spikes wrt pupil events
%         etax = siteData{5}(ev,etype)+tbins;
%         si   = find(etax(:,1)>=siteData{1}(siteData{5}(ev,1),5), 1);
%         ei   = find(etax(:,2)<=siteData{1}(siteData{5}(ev,1),6), 1, 'last');
%         if ei > si
%             sp = siteData{3}{siteData{5}(ev,1), uu+1};
%             % sp = siteData{3}{siteData{5}(ev,1), 1};
%             for tt = si:ei
%                 sdat(ee,tt) = sum( ...
%                     sp>=etax(tt,1)&sp<=etax(tt,2));
%             end
%         end
%     end
%
%
%     % get PETHs, stats
%
%     % PD magnitude measure -- P-P or max slope
%     if etype == 2
%         % difference between peak and trough
%         Zmag = diff(siteData{5}(Levents,6:7),[],2);
%     else
%         % max slope
%         Zmag = siteData{5}(Levents,9);
%     end
%     % PD at start of event
%     Z = siteData{5}(Levents,4);
%
%     % get residuals of Zmag, after linear regression with Z
%     Lr    = [Zmag>0 Zmag<0];
%     A     = [Z ones(size(Z))];
%     ZmagR = nans(size(Zmag));
%     for rr = 1:2
%         X = lscov(A(Lr(:,rr),:), Zmag(Lr(:,rr)));
%         ZmagR(Lr(:,rr)) = Zmag(Lr(:,rr))-A(Lr(:,rr),:)*X;
%     end
%
%     % select dilation/constriction events
%     Lp   = [Lr(:,1)&ZmagR>=0 (Lr(:,1)&ZmagR<0)|(Lr(:,2)&ZmagR>0) ...
%         Lr(:,2)&ZmagR<=0];
%     nlp  = size(Lp,2);
%
%     axes(axs(ss)); cla reset; hold on;
%     title(sprintf('Site=%s', sites{ss}))
%     for pp = 1:3
%         nms = nanmean(sdat(Lp(:,pp),:))';
%         ses = nanse(sdat(Lp(:,pp),:))';
%         Lg = isfinite(nms) & isfinite(ses);
%         h=patch([tax(Lg); flipud(tax(Lg))], [nms(Lg)+ses(Lg); ...
%             flipud(nms(Lg)-ses(Lg))], co{pp,2});
%         set(h,'LineStyle','none');
%         plot(tax, nms, '-', 'Color', co{pp,1});
%     end
%     xlim(xlm);
% end

%% Summary plotz, per site
% co  = {'r' 'g' 'b'};
% xlm = [tax(1) tax(end)];
% pv = 0.05;
% for ss = 1:num_sites
%
%     %% mean±sem PETHs, per trial, separated by dilation/constriction
%     subplot(num_sites,3,(ss-1)*3+1); cla reset; hold on;
%     for pp = 1:3
%         nms = nanmedian(udat{ss,1}(:,:,pp))';
%         plot(tax, nms, '-', 'Color', co{pp});
%     end
%     nm1 = udat{ss,1}(:,:,1);
%     nm3 = udat{ss,1}(:,:,3);
%     for bb = 1:ntbins
%         Lg = isfinite(nm1(:,bb)) & isfinite(nm3(:,bb));
%         if signtest(nm1(Lg,bb), nm3(Lg,bb))<pv
%             plot(tax(bb), nanmean(nm1(:,bb)), '.', 'Color', co{1});
%             plot(tax(bb), nanmean(nm3(:,bb)), '.', 'Color', co{3});
%         end
%     end
%     title(sprintf('%s, %s', sites{ss,1}, sites{ss,2}))
%     xlim(xlm);
%
%     %% count of sites with sig difference (signed)
%     subplot(num_sites,3,(ss-1)*3+2); cla reset; hold on;
%     diffs = udat{ss,1}(:,:,3)-udat{ss,1}(:,:,1);
%     nses  = sum(isfinite(udat{ss,2}(:,:,1)));
%     pyp   = nansum(udat{ss,2}(:,:,1)<pv&diffs>0)./nses;
%     pyn   = nansum(udat{ss,2}(:,:,1)<pv&diffs<0)./nses;
%     Lg    = nses>5;
%     plot(tax(Lg), pyp(Lg), 'k-');
%     plot(tax(Lg), pyn(Lg), 'k--');
%     axis([xlm 0 .6])
%
%     %% correlation coefficients
%     subplot(num_sites,3,(ss-1)*3+3); cla reset; hold on;
%     plot(tax([1 end]), [0 0], 'k:');
%     for pp = 3%1:2
%         meanv = nanmean(udat{ss,2}(:,:,pp*2))';
%         Lg    = isfinite(meanv) & sum(isfinite(udat{ss,2}(:,:,pp*2)))'>5;
%         plot(tax(Lg), meanv(Lg), 'k-');%, 'Color', co{(pp-1)*2+1});
%         for bb = 1:ntbins
%             Lg = isfinite(udat{ss,2}(:,bb,pp*2));
%             if signtest(udat{ss,2}(Lg,bb,pp*2))<pv
%                 plot(tax(bb), meanv(bb), 'k.');%, 'Color', co{(pp-1)*2+1});
%             end
%         end
%     end
%     axis([xlm -0.1 0.15])
%end
%
% xlm = [tax(1) tax(end)];
% pv = 0.01;
% for ss = 1:num_sites
%
%     %% mean±sem PETHs, per trial, separated by dilation/constriction
%     subplot(num_sites,1,ss); cla reset; hold on;
%     plot(tax([1 end]), [0 0], 'k:');
%     plot([0 0], [-1 1], 'k:');
%     nm1 = udat{ss,1}(:,:,1);
%     nm3 = udat{ss,1}(:,:,3);
%     nmd = nanmedian(nm1-nm3);
%     plot(tax, nmd, 'k-');
%     for bb = 1:ntbins
%         Lg = isfinite(nm1(:,bb)) & isfinite(nm3(:,bb));
%         if signtest(nm1(Lg,bb), nm3(Lg,bb))<pv
%             plot(tax(bb), nmd(bb), 'k.');
%         end
%     end
%     title(sprintf('%s, %s', sites{ss,1}, sites{ss,2}))
%     xlim(xlm);
% end

%
%      diffs = udat{ss,1}(:,:,3)-udat{ss,1}(:,:,1);
%      nd    = size(diffs,1);
%      for dd = 1:nd
%          mind  = min(diffs(dd,:));
%          mindi = find(diffs(dd,:)==mind,1);
%          maxd  = max(diffs(dd,:));
%          maxdi = find(diffs(dd,:)==maxd,1);
%          plot([mindi maxdi], [mind maxd], 'k.-');
%      end
%
%     adiff = abs(nanmean(udat{ss,1}(:,:,3))-nanmean(udat{ss,1}(:,:,1)));
%     Fm = find(adiff==max(adiff),1);
%     title(sprintf('Time=%.1f', tax(Fm)));
%     plot([0 60], [0 60], 'k:');
%     plot(nm1(:,Fm), nm3(:,Fm), 'k.');
%     Lp = udat{ss,2}(:,Fm,1)<0.01;
%     plot(nm1(Lp,Fm), nm3(Lp,Fm), 'r.');
%     mx = ceil(max([nm1(:,Fm); nm3(:,Fm)])/10)*10;
%     axis([0 mx 0 mx]);
%end


% %% Summary plotz, per site
% co = {'r' 'g' 'b'};
% for mm = 1:num_monkeys
%     for ss = 1:num_sites
%
%         axes(axs((mm-1)*9+ss+3)); cla reset; hold on;
%
%         %% mean±sem PETHs, per trial, separated by dilation/constriction
%         for pp = 1:3
%             nms = nanmean(udat{mm,ss,1}(:,:,pp))';
%             %ses = nanse(udat{mm,ss,1}(:,:,pp))';
%             %Lg = isfinite(nms) & isfinite(ses);
%             %h=patch([tax(Lg); flipud(tax(Lg))], [nms(Lg)+ses(Lg); ...
%             %    flipud(nms(Lg)-ses(Lg))], co{pp,2});
%             %set(h,'LineStyle','none');
%             plot(tax, nms, '-', 'Color', co{pp});
%         end
%         nm1 = udat{mm,ss,1}(:,:,1);
%         nm3 = udat{mm,ss,1}(:,:,3);
%         for bb = 1:ntbins
%             Lg = isfinite(nm1(:,bb)) & isfinite(nm3(:,bb));
%             if signtest(nm1(Lg,bb), nm3(Lg,bb))<0.01
%                 plot(tax(bb), nanmean(nm1(:,bb)), '.', 'Color', co{1});
%                 plot(tax(bb), nanmean(nm3(:,bb)), '.', 'Color', co{3});
%             end
%         end
%         xlim(xlm)
%
%         %% count of sites with sig difference
%         axes(axs((mm-1)*9+num_sites+ss+3)); cla reset; hold on;
%         nses = sum(isfinite(udat{mm,ss,2}(:,:,1)));
%         py   = nansum(udat{mm,ss,2}(:,:,1)<0.01)./nses;
%         Lg   = nses>5;
%         plot(tax(Lg), py(Lg), 'k-');
%         axis([xlm 0 .5])
%
%         %% correlation coefficients
%         axes(axs((mm-1)*9+num_sites*2+ss+3)); cla reset; hold on;
%         plot(tax([1 end]), [0 0], 'k:');
%         for pp = 1:2
%             meanv = nanmean(udat{mm,ss,2}(:,:,pp*2))';
%             Lg    = isfinite(meanv) & sum(isfinite(udat{mm,ss,2}(:,:,pp*2)))'>5;
%             plot(tax(Lg), meanv(Lg), '-', 'Color', co{pp});
%             for bb = 1:ntbins
%                 Lg = isfinite(udat{mm,ss,2}(:,bb,pp*2));
%                 if signtest(udat{mm,ss,2}(Lg,bb,pp*2))<0.01
%                     plot(tax(bb), meanv(bb), '.', 'Color', co{pp});
%                 end
%             end
%         end
%         axis([xlm -0.1 0.1])
%     end
% end



%% Junk below
%
%
%
% %% Plot fraction of sessions with sig difference between dilation/constriction
% %   PETHs, based on shuffleds
% mm=1;
% for ss = 1:num_sites
%     subplot(num_monkeys,num_sites,3+ss); cla reset; hold on;
%     % collect stats
%     nses = size(udat{mm,ss,1},1);
%     stdat = nans(nses,ntbins);
%     for ss = 1:nses
%         for bb = 1:ntbins
%             sd1 = udat{mm,ss,1}(ss,bb,1)-udat{mm,ss,1}(ss,bb,2);
%             sdr = udat{mm,ss,3}(ss,bb,1,)-udat{mm,ss,1}(ss,bb,2);
% %%
%
% % plot mean/sem partial correlation, counts of pvalues
% co = {'k' 0.8.*ones(1,3)};
% st = {'-' '--'};
% for ss = 1:num_sites
%
%     % top is count of sig units
%     subplot(num_monkeys,num_sites,ss); cla reset; hold on;
%     title(sprintf('Site=%s', sites{ss}))
%
%     for mm = 1%:num_monkeys
%         nses = sum(isfinite(udat{mm,ss,2}(:,:,1)));
%         py   = nansum(udat{mm,ss,2}(:,:,1)<0.05)./nses;
%         Lg   = nses>5;
%         plot(tax(Lg), py(Lg), '-', 'Color', co{mm});
%     end
%     axis([-1200 200 0 1])
%
%     % correlations -- separate per dilation/constriction
%     subplot(num_monkeys,num_sites,num_sites+ss); cla reset; hold on;
%     plot(tax([1 end]), [0 0], 'k--');
%
%     for mm = 1%:num_monkeys
%         for pp = 1:2
%             meanv = nanmean(udat{mm,ss,2}(:,:,pp*2))';
%             %            semv  = nanse(udat{mm,ss,2}(:,:,pp*2))';
%             Lg = isfinite(meanv) & sum(isfinite(udat{mm,ss,2}(:,:,pp*2)))'>5;
%             %            h=patch([tax(Lg); flipud(tax(Lg))], [meanv(Lg)+semv(Lg); ...
%             %                flipud(meanv(Lg)-semv(Lg))], 0.9*ones(1,3));
%             %            set(h,'LineStyle','none');
%             plot(tax(Lg), meanv(Lg), st{mm}, 'Color', co{pp});
%         end
%     end
%     axis([-1200 200 -0.1 0.14])
% end
%
% % RANDOMIZED DATA
% % plot mean/sem partial correlation, counts of pvalues
% cop = {'r' [1 .8 .8]; 'b' [.8 .8 1]};
% co = {'k' 0.8.*ones(1,3)};
% st = {'-' '--'};
% for rr = 1:500
%     for ss = 1:num_sites
%
%         % top is mean PETHs
%         subplot(3,num_sites,ss); cla reset; hold on;
%         for pp = 1:2
%             nms = nanmean(udat{mm,ss,3}(:,:,pp,rr))';
%             % ses = nanse(udat{mm,ss,3}(:,:,pp,rr))';
%             %Lg = isfinite(nms) & isfinite(ses);
%             %h=patch([tax(Lg); flipud(tax(Lg))], [nms(Lg)+ses(Lg); ...
%             %    flipud(nms(Lg)-ses(Lg))], co{pp,2});
%             %set(h,'LineStyle','none');
%             plot(tax, nms, '-', 'Color', cop{pp,1});
%         end
%         xlim([-1200 200])
%
%         % next is count of sig units
%         subplot(3,num_sites,ss+3); cla reset; hold on;
%         title(sprintf('Site=%s', sites{ss}))
%
%         for mm = 1%:num_monkeys
%             nses = sum(isfinite(udat{mm,ss,3}(:,:,1,rr)));
%             py   = nansum(udat{mm,ss,3}(:,:,1,rr)<0.05)./nses;
%             Lg   = nses>5;
%             plot(tax(Lg), py(Lg), '-', 'Color', co{mm});
%         end
%         axis([-1200 200 0 1])
%
%         % correlations -- separate per dilation/constriction
%         subplot(3,num_sites,6+ss); cla reset; hold on;
%         plot(tax([1 end]), [0 0], 'k--');
%
%         for mm = 1%:num_monkeys
%             for pp = 1:2
%                 meanv = nanmean(udat{mm,ss,3}(:,:,pp*2,rr))';
%                 %            semv  = nanse(udat{mm,ss,2}(:,:,pp*2))';
%                 Lg = isfinite(meanv) & sum(isfinite(udat{mm,ss,3}(:,:,pp*2,rr)))'>5;
%                 %            h=patch([tax(Lg); flipud(tax(Lg))], [meanv(Lg)+semv(Lg); ...
%                 %                flipud(meanv(Lg)-semv(Lg))], 0.9*ones(1,3));
%                 %            set(h,'LineStyle','none');
%                 plot(tax(Lg), meanv(Lg), st{mm}, 'Color', co{pp});
%             end
%         end
%         axis([-1200 200 -0.1 0.14])
%     end
%     r = input('ext')
% end
%
%
%
%
%
%
%
% % correlations
% co = {'k' 'c' 'r' 'g' 'b' 'y'};
% for mm = 1:num_monkeys
%     for ss = 1:num_sites
%         subplot(num_monkeys,num_sites,(mm-1)*num_sites+ss); cla reset; hold on;
%         title(sprintf('Monkey=%s, Site=%s', monkeys{mm}, sites{ss}))
%         if ~isempty(udat{mm,ss,1})
%             for pp = 3:6
%                 plot(tax, nanmean(udat{mm,ss,uu,1,1}(:,:,pp)), '-', 'Color', co{pp});
%             end
%         end
%         xlim([tax(1) tax(end)])
%     end
% end
% % p values
% co = {'k' 'r'}; % median, correlation
% for mm = 1:num_monkeys
%     for ss = 1:num_sites
%         subplot(num_monkeys,num_sites,(mm-1)*num_sites+ss); cla reset; hold on;
%         title(sprintf('Monkey=%s, Site=%s, %d sessions', ...
%             monkeys{mm}, sites{ss}, size(udat{mm,ss,uu,1,1},1)))
%         for pp = 1:2
%             plot(tax, nansum(udat{mm,ss,uu,pp,2}<0.05), '-', 'Color', co{pp});
%         end
%         xlim([tax(1) tax(end)])
%     end
% end
% % corr
% for mm = 1:num_monkeys
%     for ss = 1:num_sites
%         subplot(num_monkeys,num_sites,(mm-1)*num_sites+ss); cla reset; hold on;
%         title(sprintf('Monkey=%s, Site=%s', monkeys{mm}, sites{ss}))
%         plot(tax, nanmean(udat{mm,ss,uu,2,1}), 'k-');
%         xlim([tax(1) tax(end)])
%     end
% end
%
%
% uu = 1;
% for mm = 1:num_monkeys
%     for ss = 1:num_sites
%         subplot(num_monkeys,num_sites,(mm-1)*num_sites+ss); cla reset; hold on;
%         title(sprintf('Monkey=%s, Site=%s', monkeys{mm}, sites{ss}))
%         if ~isempty(udat{mm,ss,1})
%             %             imagesc(tax, [], nanrunmean( ...
%             %                 (udat{mm,ss,uu}(:,:,3)-udat{mm,ss,uu}(:,:,5))', 5)', [-.4 .4]);
%             %             axis([tax(1) tax(end) 1 size(udat{mm,ss,uu},1)]);
%             % end
%             for ff = 1:size(udat{mm,ss,uu},1)
%                 tdat = nanrunmean( ...
%                     udat{mm,ss,uu}(ff,:,3)-udat{mm,ss,uu}(ff,:,5), 5);
%                 tmax = max(tdat);
%                 plot(tax(find(tdat==tmax,1)), tmax, 'k.');
%             end
%         end
%         axis([tax(1) tax(end) 0 0.7])
%     end
% end
%
%
%
% co = {'r' 'g' 'b' 'k'};
% for mm = 1:num_monkeys
%     for ss = 1:num_sites
%         subplot(num_monkeys,num_sites,(mm-1)*num_sites+ss); cla reset; hold on;
%         for uu = 1:2
%             plot(repmat(tax(udat{mm,ss,uu}(:,1)),1,2)', ...
%                 udat{mm,ss,uu}(:,3:4)', '-', 'Color', co{uu});
%             Lg = udat{mm,ss,uu}(:,3)>udat{mm,ss,uu}(:,4);
%             plot(repmat(tax(udat{mm,ss,uu}(Lg,1)),1,2)', ...
%                 udat{mm,ss,uu}(Lg,3:4)', '-', 'Color', co{uu}, 'LineWidth', 2);
%             plot(tax(udat{mm,ss,uu}(:,1)), udat{mm,ss,uu}(:,3), ...
%                 'x', 'Color', co{uu});
%             plot(tax(udat{mm,ss,uu}(:,1)), udat{mm,ss,uu}(:,2), ...
%                 '.', 'Color', co{uu});
%         end
%         axis([tax(1) tax(end) 0 2])
%     end
% end
%
%
%
%
% % find/save max difference, also whether magnituide of
% % dilations (max slope) reflect mean spike rate
% Ldi = pdat(:,7)>0;
% Mdi = nanmedian(pdat(Ldi,7));
% for uu = 1:num_units
%     ui    = min(uu,2);
%     pdiff = nanmean(sdat(Ldi,:,uu))-nanmean(sdat(~Ldi,:,uu));
%     mpd   = max(pdiff);
%     pi    = find(pdiff==mpd,1);
%     pd1   = nanmean(sdat(Ldi&pdat(:,7)>=Mdi,pi,uu))-nanmean(sdat(~Ldi,pi,uu));
%     pd2   = nanmean(sdat(Ldi&pdat(:,7)<Mdi,pi,uu))-nanmean(sdat(~Ldi,pi,uu));
%     udat{mm,ss,ui} = cat(1, udat{mm,ss,ui}, ...
%         [pi mpd pd1 pd2 ff]);
% end
%
% %                 % separate into dilation/constriction, get PETHs
% %                 Ldi = pdat(:,7)>0;
% %                 Mdi = nanmedian(pdat(Ldi,7));
% %                 Mco = nanmedian(pdat(~Ldi,7));
% %                 Lds = [Ldi ~Ldi Ldi&pdat(:,7)>=Mdi Ldi&pdat(:,7)<Mdi ...
% %                     ~Ldi&pdat(:,7)>=Mco ~Ldi&pdat(:,7)<Mco];
% %
% %                 co = {'r-' 'b-' 'r--' 'r:' 'b--' 'b:'};
% %                 for uu = 1:num_units
% %                     subplot(num_units,1,uu); cla reset; hold on;
% %                     for ll = 1:size(Lds,2)
% %                         plot(tax, nanmean(sdat(Lds(:,ll),:,uu)), co{ll});
% %                     end
% %                     plot(tax, nanmean(sdat(Lds(:,1),:,uu))-nanmean(sdat(Lds(:,2),:,uu)), 'k-', 'LineWidth', 2)
% %                 end
% %
% %                 r = input('next');
% end
% end
% end
% end
%
% uu = 2
% co = {'r' 'k'};
% for mm = 1:num_monkeys
%     for ss = 1:num_sites
%         subplot(num_monkeys,num_sites,(mm-1)*num_sites+ss); cla reset; hold on;
%         title(sprintf('Monkey=%s, Site=%s', monkeys{mm}, sites{ss}))
%         if ~isempty(udat{mm,ss,1})
%             for uu = 1:2
%                 plot(repmat(tax(udat{mm,ss,uu}(:,1)),1,2)', ...
%                     udat{mm,ss,uu}(:,3:4)', '-', 'Color', co{uu});
%                 Lg = udat{mm,ss,uu}(:,3)>udat{mm,ss,uu}(:,4);
%                 plot(repmat(tax(udat{mm,ss,uu}(Lg,1)),1,2)', ...
%                     udat{mm,ss,uu}(Lg,3:4)', '-', 'Color', co{uu}, 'LineWidth', 2);
%                 plot(tax(udat{mm,ss,uu}(:,1)), udat{mm,ss,uu}(:,3), ...
%                     'x', 'Color', co{uu});
%                 plot(tax(udat{mm,ss,uu}(:,1)), udat{mm,ss,uu}(:,2), ...
%                     '.', 'Color', co{uu});
%             end
%         end
%         axis([tax(1) tax(end) 0 2])
%     end
% end
