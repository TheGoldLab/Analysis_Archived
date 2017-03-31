 function figLCP_pupilStats(num, session, collectData)
% function figLCP_pupilStats(num, session, collectData)
%
% Basic pupil stats:
%   A. Mean PD (raw and corrected) per trial as a function of trial 
%       number within a session. 
%   B. Close-up of one trial from A, showing PD versus time. 
%   C-D, histogram of counts of PD dilation and constriction events, 
%       by duration, for the two different monkeys. 
%   E-F, Relationship between baseline PD and PD event magnitude for 
%       the two monkeys.
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
if nargin < 1
    num = 2;
end

if nargin < 2 || isempty(session)
    % monkey, site, session, trial
    session = {'Oz' 'LC' 29 61};
end

monkeys = {...
    'Oz'      {'LC' 'IC' 'SC'}; ...
    'Cicero'  {'LC' 'IC' 'SC'}; ...
    'Sprout'  {'ACC' 'PCC'}; ...
    'Atticus' {'ACC'}; ...
    'Cheetah' {'PCC'}; ...
    };
nm = size(monkeys,1);

%% Set up figure
wid        = 17.6; % total width
ht         = 3.5;
cols       = {1,1,3};
[axs,~] = getPLOT_axes(num, wid, ht, cols, 2, 2, [], 'Joshi et al', true);
set(axs,'Units','normalized');

% Gray color
gr = 0.8*ones(1,3);

%% Load example file
[base_dir, fnames] = getLCP_cleanDataDir(session{1}, session{2});
load(fullfile(base_dir, fnames{session{3}}));

%% First panel is example session, average PD per trial
Fnb    = find(~isfinite(siteData{1}(:,4)) & ~isfinite(siteData{1}(:,9)));
num_nb = length(Fnb);
wrt0   = siteData{1}(1,7);
tax    = 0:size(siteData{2},2)-1;

axes(axs(1)); cla reset; hold on;
aa = [];
for tt = 1:num_nb
    aa = cat(1, aa, siteData{2}(Fnb(tt),:,4)');
end

    plot(((siteData{1}(Fnb(tt),7)-wrt0)+tax)./60000, ...
        siteData{2}(Fnb(tt),:,3), '-', 'Color', gr);
    plot(((siteData{1}(Fnb(tt),7)-wrt0)+tax)./60000, ...
        siteData{2}(Fnb(tt),:,4), 'k-');
end
xlabel('Time (min)');
ylabel('PD (z-score)')
ylim([-4 4]);

%% Second panel is single-trial example
axes(axs(2)); cla reset; hold on;
plot(siteData{2}(Fnb(session{4}),:,3), '-', 'Color', gr);
ys=siteData{2}(Fnb(session{4}),:,4);
plot(ys, 'k-');
Fevt = find(siteData{5}(:,1)==Fnb(session{4}));
for ee = 1:length(Fevt)
    plot(siteData{5}(Fevt(ee),2), siteData{5}(Fevt(ee),6), 'ko')
    plot(siteData{5}(Fevt(ee),8), ys(siteData{5}(Fevt(ee),8)), 'kx')
end
xlabel('Time (msec)');
ylabel('PD (z-score)')

% inset with spectral analysis
axes('position', [.68 .535 .22 .07]); hold on;
ffn = 8192;
ys  = siteData{2}(Fnb(session{4}),:,4);
ys  = ys(isfinite(ys));
[pxx,f] = periodogram(ys,[],ffn,1000,'power');

nnb = length(Fnb);
fdat = nans(nnb, ffn/2+1);
for ii = 1:nnb
    ys  = siteData{2}(Fnb(ii),:,4);
    ys  = ys(isfinite(ys));
    [pxxt,ft] = periodogram(ys,[],ffn,1000,'power');
    fdat(ii,:) = pxxt;
end
mn = nanmean(10*log10(fdat))';
se = nanse(10*log10(fdat))';
h=patch([ft; flip(ft)], [mn-se; flip(mn+se)], 0.8.*ones(1,3));
set(h,'LineStyle', 'none');
plot(ft,mn,'k-', 'LineWidth', 2); % mean across trials
plot(f,10*log10(pxx),'k-'); % example trial
xlim([0 10]);
xlabel('Frequency (Hz)')
ylabel('Power (dB)')


hcn = dsp.ColoredNoise('InverseFrequencyPower',1,'SamplesPerFrame',4096);
x = step(hcn);
Fs = 1000;
[p_pxxt, p_ft] = periodogram(x, [], ffn, 1000, 'power');
cla reset; hold on;
plot(log2(p_ft), 10*log10(p_pxxt))
plot(log2(ft), mn)

[~,fF] = pwelch(x,hamming(128),[],[],Fs,'psd');
PSDPink = 1./fF(2:end);
plot(log2(fF(2:end)),10*log10(PSDPink),'r','linewidth',2)

%% Collect some data:
%   event stats: duration, baseline, size, vmax
%   microsaccade stats: histogram of "phase" wrt pupil event
pbins  = linspace(0,1,5)';
npbins = length(pbins);
if nargin > 3 && collectData
    udat   = cell(nm,2);    
    for mm = 1:nm
        ns = size(monkeys{mm,2},2);
        sdat = zeros(npbins, 2);
        for ss = 1:ns
            [base_dir, fnames] = getLCP_cleanDataDir(monkeys{mm,1}, monkeys{mm,2}{ss});
            for ff = 1:length(fnames)
                disp(sprintf('%s: %s, file %d/%d', ...
                    monkeys{mm,1}, monkeys{mm,2}{ss}, ff, length(fnames)))
                load(fullfile(base_dir, fnames{ff}));
                Lnb = ~isfinite(siteData{1}(:,4)) & ~isfinite(siteData{1}(:,9));
                Fnb = find(Lnb);
                
                % 1: event stats
                udat{mm,1} = cat(1, udat{mm,1}, [...
                    siteData{5}(Lnb,3)-siteData{5}(Lnb,2), ...  % duration
                    siteData{5}(Lnb,6), ...                     % baseline
                    siteData{5}(Lnb,7)-siteData{5}(Lnb,6), ...  % magnitude
                    siteData{5}(Lnb,9)]);                       % vmax
                
                % 2: microsaccade stats
                Lgm = isfinite(siteData{6}(:,6)) & ...
                    ismember(siteData{6}(:,1), Fnb);
                if any(Lgm)
                    vals = siteData{6}(Lgm,6);
                    vals(vals==1) = 0.999;
                    if sum(Lgm) == 1
                        sdat(:,1) = sdat(:,1) + histc(vals, pbins)';
                    else
                        sdat(:,1) = sdat(:,1) + histc(vals, pbins);
                    end
                    % count total # of pupil events
                    sdat(:,2) = sdat(:,2) + ...
                        sum(ismember(siteData{5}(:,1), Fnb));
                end
            end
        end
        udat{mm,2} = [sdat(:,1)./sdat(:,2) sdat(:,2)];
    end
    % save data to file
    FS_saveProjectFile('2013_LCPupil', 'pupilStats', udat);
else
    % load data from file
    udat = FS_loadProjectFile('2013_LCPupil', 'pupilStats');
end

%% Third panel is histograms of PD "event" durations
xs    = (0:20:2000)';
mvy   = [1300 -1600];
hists = zeros(length(xs), 2);
axes(axs(3)); cla reset; hold on;
pv = [];
nv = [];
for mm = 1:nm
    % get monkey symbol
    [~, sy] = getLCP_colorsAndSymbol(monkeys{mm,1}, monkeys{mm,2}{1});
    for ii = 1:2
        if ii == 1
            vals = udat{mm,1}(udat{mm,1}(:,3)>0,1);
            pv = cat(1, pv, [vals mm.*ones(size(vals))]);
        else
            vals = udat{mm,1}(udat{mm,1}(:,3)<0,1);
            nv = cat(1, nv, [vals mm.*ones(size(vals))]);
        end
        hists(:,ii) = hists(:,ii) + hist(vals,xs)';
        plot(nanmedian(vals), mvy(ii), sy, 'Color', 'k', 'MarkerSize', 5);
        disp(sprintf('Monkey %s: med%d=%.1f', monkeys{mm,1}, ...
            ii, nanmedian(vals)))
    end
end
% prctile(cat(1, ...
%     udat{1,1}(udat{1,1}(:,3)<0,1), ...
%     udat{2,1}(udat{2,1}(:,3)<0,1), ...
%     udat{3,1}(udat{3,1}(:,3)<0,1), ...
%     udat{4,1}(udat{4,1}(:,3)<0,1), ...
%     udat{5,1}(udat{5,1}(:,3)<0,1)), [50 25 75])

% add hr marker
H(1) = bar(xs,hists(:,1));
H(2) = bar(xs,-hists(:,2));
set(H,'FaceColor','k')
plot([0 2000], [0 0], 'w-');
axis([0 1200 -1800 1400])
xlabel('Duration (ms)')
ylabel('Count')
kruskalwallis(pv(:,1), pv(:,2))

disp(sprintf('Positive KW test, p=%.2f', kruskalwallis(pv(:,1), pv(:,2), 'off')))
disp(sprintf('Negative KW test, p=%.2f', kruskalwallis(nv(:,1), nv(:,2), 'off')))

%% Fourth panel is density maps of event size vs baseline
% just show example session, report slopes 
axes(axs(4)); cla reset; hold on;
bdat = nans(nm, 2, 2);
for mm = 1:nm
    
    vals = udat{mm,1}(:,2:3);

    % dilation regression
    dval = vals(vals(:,2)>0,:);
    bdat(mm, :, 1) = lscov(cat(2, ones(length(dval),1), dval(:,1)), ...
        dval(:,2));
    
    % constriction regression
    cval = vals(vals(:,2)<0,:);
    bdat(mm, :, 2) = lscov(cat(2, ones(length(cval),1), cval(:,1)), ...
        cval(:,2));
    
    if mm == 1
        plot([-10 10], [0 0], 'k:')
        xs = (-5:.01:5)';
        A  = cat(2, ones(size(xs)), xs);
        ys = A*bdat(mm,:,1)';
        plot(dval(:,1), dval(:,2), 'k.');
        plot(xs(ys>0), ys(ys>0), '-', 'Color', gr, 'LineWidth', 4);

        ys = A*bdat(mm,:,2)';
        plot(cval(:,1), cval(:,2), 'k.')
        plot(xs(ys<0), ys(ys<0), '--', 'Color', gr, 'LineWidth', 4);
        axis([-4 4 -4 4])
    end    
end

%% COMMENTED OUT IS SUMMARY FOR ALL MONKS
% tmpv = [];
% bdat = nans(nm, 2, 2);
% for mm = 1:nm
%     
%     vals = udat{mm,1}(:,2:3);
% 
%     % dilation regression
%     dval = vals(vals(:,2)>0,:);
%     bdat(mm, :, 1) = lscov(cat(2, ones(length(dval),1), dval(:,1)), ...
%         dval(:,2));
%     
%     % constriction regression
%     cval = vals(vals(:,2)<0,:);
%     bdat(mm, :, 2) = lscov(cat(2, ones(length(cval),1), cval(:,1)), ...
%         cval(:,2));
%     
%     % save it all
%     tmpv = cat(1, tmpv, udat{mm,1}(:,2:3));
% end
% 
% axes(axs(4)); cla reset; hold on;
% %% PLOT SEPARATELY!!!
% hist3(tmpv,[200 200]);
% set(gca, 'View', [0 90])
% set(gcf,'renderer','opengl');
% set(get(gca,'child'), ...
%     'LineStyle', 'none', 'FaceColor','interp','CDataMode','auto');
% xs = (-2.5:.01:2.5)';
% A  = cat(2, ones(size(xs)), xs);
% for mm = 1:nm
%     [co, sy] = getLCP_colorsAndSymbol(monkeys{mm,1}, monkeys{mm,2}{1});
%     for ii = 1:2
%         ys = A*bdat(mm,:,ii)';
%         if ii == 1
%             plot(xs(ys>0), ys(ys>0), '-', 'Color', co{mm});
%         else
%             plot(xs(ys<0), ys(ys<0), '--', 'Color', co{mm});
%         end
%     end
% end
% axis([-2.5 2.5 -2.5 2.5])
% xlabel('Baseline PD (z-score)');
% ylabel('Phasic PD (z-score)');

%% Fifth panel is microsaccade
xax = pbins(1:end-1)+mean(diff(pbins))./2;
offs = -.08:.04:.08;
axes(axs(5)); cla reset; hold on;
for mm = 1:nm
    %[co, sy] = getLCP_colorsAndSymbol(monkeys{mm,1}, monkeys{mm,2}{1});
    bar(xax+offs(mm), udat{mm,2}(1:end-1,1), 0.15, 'k');
end
set(gca, 'XTick', 0.2:0.2:0.8, 'XTickLabel', ...
    {'0-90' '90-180' '180-270' '270-360'})
axis([0 1 0 1]);
xlabel('Phase (deg)')
ylabel('Fraction of pupil events')

% Rayleigh's test for uniformity
angs = (pbins(1:end-1)+diff(pbins)/2).*2.*pi;
for mm = 1:nm
    [p,z] = circ_rtest(angs, udat{mm,2}(1:end-1,1)./udat{mm,2}(1:end-1,1));
    disp([mm p])
end



%% JUNK BELOW
%  density maps of vmax vs event size
% for mm = 1:2
%     subplot(4,2,6+mm); cla reset; hold on;
%     hist3(vals{mm}(:,2:3),[200 200]);
%     set(gca, 'View', [0 90])
%     set(gcf,'renderer','opengl');
%     set(get(gca,'child'), ...
%         'LineStyle', 'none', 'FaceColor','interp','CDataMode','auto');
%     %axis([-2.5 2.5 -2.5 2.5])
% end
% 
%             % 2: periodogram per trial, session
%             wrt0  = siteData{1}(1,7);
%             adat = nans(nnb, size(siteData{2},2)+1);
%             for tt = 1:nnb
%                 Ltt = isfinite(siteData{2}(Fnb(tt),:,4));
%                 [A,LAGS] = xcorr(siteData{2}(Fnb(tt),Ltt,4), 'coeff');
%                 adat(tt,1:sum(Ltt)) = A(LAGS>=0);
%             end
%             
%                 subplot(2,1,1); cla reset; hold on;
%                 plot(siteData{2}(Fnb(tt),Ltt,4));
%                 subplot(2,1,2); cla reset; hold on;
%                 plot(XX);
%                 r = input('next')
%             end
%                 
%             
%             tax   = (0:size(siteData{2},2)-1)';
%             ptdat = [];
%             for tt = find(Lnb)'
%                 Ltt = isfinite(siteData{2}(tt,:,4));
%                 ptdat = cat(1, ptdat, ...
%                     [((siteData{1}(tt,7)-wrt0)+tax(Ltt))./1000, ...
%                     siteData{2}(tt,Ltt,4)']);
%             end
%             %[Pf,Ff,ALPHAf] = lomb(ptdat(:,2), ptdat(:,1));
%             [P,F,ALPHA] = fastlomb(ptdat(:,2), ptdat(:,1));
%             vals{mm,2} = cat(1, vals{mm,2}, interp1(F,P,xi,'spline'));
%         end
%     end
% end
