% plot reinforcement learning figure 4 - MT choice probability

% Matlab script for plotting Fig 4
% for submission to Nature Neuroscience
%
% created by jcl on 08/07/08

%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 145/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 90/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 7;    % axis font size
lfs = 7;    % label font size
ffs = 12;   % figure font size
tfs = 10;    % text font size
lw  = 1.5;  % line width
elw = 0.75;
ms  = 2;    % marker size

% define line color
lc1  = [0.2 0.2 0.2];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.2 0.2];
lc2e = [1 0.75 0.75];


asize = 0.21;



%% compute choice probability with decaying correlation structure
% get data

recompute = 0;



Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8;
DIRS      = [-90,90];
SIMNUM    = [81:90];

if strcmp(Monk, 'Cy')
    N = 115985;
end

trialnum  = floor([2000 40000 N]/1000)+1;


if recompute
    % get weights for early/mid/late trial, for each simulation
    w   = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    ind = nans(NSEN, length(SIMNUM), length(trialnum));
    for i = 1:length(SIMNUM)
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(eln_) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);

        for j = 1:length(trialnum)
            w(:,:,i,j) = W(:,:,trialnum(j));
            ind(:,i,j) = MTindex;
        end
    end

    [fitsp, semsp] = getML_neurometricROCT([Monk 'PRe_MT.txt'], [], 1, [0 1], 0, 0);
    [fits, sems]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], [], 1, [0 1], 0, 0);
    [pp, np, fp  ] = getML_SigNoise2([Monk 'PRe_MT.txt'], 1, [0 1], 1, 0);
    [p, n, f]      = getML_SigNoise2([Monk 'TRain_MT.txt'], 1, [0 1], 1, 0);

    % combine pre- and during training data
    PREF = [pp p];
    NULL = [np n];
    FANO = [fp f];
    TH   = getML_RLSimulatedThForMT(Monk, 10000, 0);
    if strcmp(Monk, 'Cy')
        Lgd      = ~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
        Lgd(178) = logical(0);  % this cell is weird, with very small gain to both pref and null
    else
        Lgd      =~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
    end    
    PREF = PREF(:,Lgd);
    NULL = NULL(:,Lgd);
    FANO = FANO(:,Lgd);
    TH   = TH(Lgd);
    
    % sort data by threshold coherence
    [TH,xxx] = sort(TH);
    PREF     = PREF(:,xxx);
    NULL     = NULL(:,xxx);
    FANO     = FANO(:,xxx);


    %%
    % get average weights for MT neurons across simulation
    wavg = nans(length(PREF), NTUNE, length(trialnum));
    for i = 1:length(trialnum)
        for j = 1:length(PREF)
            w_ = [];
            for k = 1:length(SIMNUM)
                L = ind(:,k,i)==j;
                w_          = [w_; w(L,:,k,i)];
            end
            wavg(j,:,i) = nanmean(w_,1);
        end
    end
            
    
    %% assign average weights back to simulations
    wm = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    for i = 1:length(SIMNUM)
        for j = 1:length(trialnum)
            wm(:,:,i,j) = wavg(ind(:,i,j),:,j);
        end
    end    
end




% get weights for early/mid/late trial, for each simulation
cp    = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
MTind = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
if recompute
    for i = 1:length(SIMNUM)
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(eln_) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);
        clear FIRA
        pack

        for j = 1:length(trialnum)
            cp(:,:,i,j)  = getML_RLCPForWeight2(Monk, NSEN, NTUNE, wm(:,:,i,j), ind(:,i,j), 2000, sprintf('%d-%d', i, j), 1);
            MTind(:,:,i,j) = repmat(MTindex',1,NTUNE);
        end
    end
    th = TH(MTind);
else
    for i = 1:length(SIMNUM)
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(eln_) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);
        clear FIRA
        pack

        for j = 1:length(trialnum)
            cp(:,:,i,j)  = getML_RLCPForWeight2(Monk, NSEN, NTUNE, [], [], 2000, sprintf('%d-%d', i, j), 0);
            MTind(:,:,i,j) = repmat(MTindex',1,NTUNE);
        end
    end
end

[fitsp, semsp] = getML_neurometricROCT([Monk 'PRe_MT.txt'], [], 1, [0 1], 0, 0);
[fits, sems]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], [], 1, [0 1], 0, 0);
[pp, np, fp  ] = getML_SigNoise2([Monk 'PRe_MT.txt'], 1, [0 1], 1, 0);
[p, n, f]      = getML_SigNoise2([Monk 'TRain_MT.txt'], 1, [0 1], 1, 0);

% combine pre- and during training data
PREF = [pp p];
NULL = [np n];
FANO = [fp f];
TH   = getML_RLSimulatedThForMT(Monk, 10000, 0);
if strcmp(Monk, 'Cy')
    Lgd      = ~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
    Lgd(178) = logical(0);  % this cell is weird, with very small gain to both pref and null
else
    Lgd      =~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
end
PREF = PREF(:,Lgd);
NULL = NULL(:,Lgd);
FANO = FANO(:,Lgd);
TH   = TH(Lgd);

% sort data by threshold coherence
[TH,xxx] = sort(TH);
PREF     = PREF(:,xxx);
NULL     = NULL(:,xxx);
FANO     = FANO(:,xxx);

th = TH(MTind);


% binned cp by thresholds
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
cpm    = nans(length(THBINS), NTUNE, length(trialnum));
cpd    = nans(length(THBINS), NTUNE, length(trialnum));
cpn    = nans(length(THBINS), NTUNE, length(trialnum));
   
% average weight for each threshold bin
for i = 1:length(trialnum)
    for j = 1:NTUNE
        cp_  = cp(:,j,:,i);
        cp_  = cp_(:);
        th_  = 100*th(:,j,:,i);
        th_  = th_(:);
        
        for k = 1:length(THBINS)
            L          = th_>=THBINS(k,1) & th_<THBINS(k,2);
            cpm(k,j,i) = nanmean(cp_(L));
            cpd(k,j,i) = nanse(cp_(L));
            cpn(k,j,i) = sum(L);
            
        end
    end
end



%% Figure 4a: correlation structure of an example cell with ranked
%  sensitivity of 100, and direction tuning of 0 degree
th = annotation('textbox', [0 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


asize = 0.2;

DTUNE      = 360/NTUNE;
dirtune    = repmat(linspace(-170,180,NTUNE),NSEN,1);

rsen       = 0.18;       % max correlation between neurons with different sensitivity
gsen       = rsen*ones(NSEN,NTUNE);
gsen(gsen<0) = 0;
gsen       = gsen/max(max(gsen));

gtune      = exp(-(abs(dirtune))/30);

corrstruct = rsen.*gsen.*gtune;




% plot
clim_ = [0 rsen];
DTUNE = 360/NTUNE;
pdir  = linspace(-170,180,NTUNE);

% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.08 0.58 asize/1.3 asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:200, corrstruct)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'out', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [1 50:50:200], 'clim', clim_)
xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Ranked sensitivity', 'FontSize', lfs)
title('Correlation structure', 'FontSize', tfs)



%% Figure 4b: Choice probablility early mid and late in training, with
%             decaying correlation structure
th = annotation('textbox', [0.26 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


asize = 0.2;

DTUNE  = 360/NTUNE;
pdir   = linspace(-180,180-DTUNE,NTUNE);
thm    = nangeomean(THBINS,2);
ind    = find(pdir==0); 

bdata  = [0.5248 -0.0064; 0.6153 -0.0290; 0.6674 -0.0397];  % slope of MT CP from real data, value from LAW GOLD 2007 Figure 3



set(axes, 'Units', 'Normalized', 'Position', [0.125+1*(0.02+asize) 0.58 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([thm thm]',[cpm(:,ind,1)-cpd(:,ind,1), cpm(:,ind,1)+cpd(:,ind,1)]','-', 'Color', lc1e, 'LineWidth', elw)
plot(thm,cpm(:,ind,1),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd      = ~isnan(cpm(:,ind,1)) & cpd(:,ind,1)~=0;
[b, bi, h, p, F, df] = nestedFW(cpm(Lgd,ind,1), 1./cpn(Lgd,ind,1), [ones(size(log(thm(find(Lgd))))) log(thm(find(Lgd)))]);
plot(2:100, b(1)+b(2)*log(2:100), 'k-', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickDir', 'out', ...
         'xscale', 'log', 'xlim', [4.5 100], 'ylim', [0.49 0.56], ...
         'XTick', [5 10 20 40 80], 'YTick', [0.5 0.53 0.56])
xlabel('Threshold (% coh)', 'FontSize', lfs)
ylabel('Choice Probability', 'FontSize', lfs)
title('Early', 'FontSize', tfs)
  


set(axes, 'Units', 'Normalized', 'Position', [0.125+2*(0.02+asize) 0.58 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([thm thm]',[cpm(:,ind,2)-cpd(:,ind,2), cpm(:,ind,2)+cpd(:,ind,2)]','-', 'Color', lc1e, 'LineWidth', elw)
plot(thm,cpm(:,ind,2),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd      = ~isnan(cpm(:,ind,2)) & cpd(:,ind,2)~=0;
[b, bi, h, p, F, df] = nestedFW(cpm(Lgd,ind,2), 1./cpn(Lgd,ind,2), [ones(size(log(thm(find(Lgd))))) log(thm(find(Lgd)))]);
plot(2:100, b(1)+b(2)*log(2:100), 'k-', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickDir', 'out', ...
         'xscale', 'log', 'xlim', [4.5 100], 'ylim', [0.49 0.56], ...
         'XTick', [5 10 20 40 80], 'YTick', [0.5 0.53 0.56], 'YTickLabel', {})
xlabel('Threshold (% coh)', 'FontSize', lfs)
title('Mid', 'FontSize', tfs)
  
   


set(axes, 'Units', 'Normalized', 'Position', [0.125+3*(0.02+asize) 0.58 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([thm thm]',[cpm(:,ind,3)-cpd(:,ind,3), cpm(:,ind,3)+cpd(:,ind,3)]','-', 'Color', lc1e, 'LineWidth', elw)
plot(thm,cpm(:,ind,3),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd      = ~isnan(cpm(:,ind,3)) & cpd(:,ind,3)~=0;
[b, bi, h, p, F, df] = nestedFW(cpm(Lgd,ind,3), 1./cpn(Lgd,ind,3), [ones(size(log(thm(find(Lgd))))) log(thm(find(Lgd)))]);
plot(2:100, b(1)+b(2)*log(2:100), 'k-', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickDir', 'out', ...
         'xscale', 'log', 'xlim', [4.5 100], 'ylim', [0.49 0.56], ...
         'XTick', [5 10 20 40 80], 'YTick', [0.5 0.53 0.56], 'YTickLabel', {})
xlabel('Threshold (% coh)', 'FontSize', lfs)
title('Late', 'FontSize', tfs)
  
   




%% compute choice probability with decaying correlation structure
% get data

recompute = 0;



Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8;
DIRS      = [-90,90];
SIMNUM    = [81:90];

if strcmp(Monk, 'Cy')
    N = 115985;
end

trialnum  = floor([2000 40000 N]/1000)+1;


if recompute
    % get weights for early/mid/late trial, for each simulation
    w   = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    ind = nans(NSEN, length(SIMNUM), length(trialnum));
    for i = 1:length(SIMNUM)
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(eln_) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);

        for j = 1:length(trialnum)
            w(:,:,i,j) = W(:,:,trialnum(j));
            ind(:,i,j) = MTindex;
        end
    end

    [fitsp, semsp] = getML_neurometricROCT([Monk 'PRe_MT.txt'], [], 1, [0 1], 0, 0);
    [fits, sems]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], [], 1, [0 1], 0, 0);
    [pp, np, fp  ] = getML_SigNoise2([Monk 'PRe_MT.txt'], 1, [0 1], 1, 0);
    [p, n, f]      = getML_SigNoise2([Monk 'TRain_MT.txt'], 1, [0 1], 1, 0);

    % combine pre- and during training data
    PREF = [pp p];
    NULL = [np n];
    FANO = [fp f];
    TH   = getML_RLSimulatedThForMT(Monk, 10000, 0);
    if strcmp(Monk, 'Cy')
        Lgd      = ~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
        Lgd(178) = logical(0);  % this cell is weird, with very small gain to both pref and null
    else
        Lgd      =~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
    end    
    PREF = PREF(:,Lgd);
    NULL = NULL(:,Lgd);
    FANO = FANO(:,Lgd);
    TH   = TH(Lgd);
    
    % sort data by threshold coherence
    [TH,xxx] = sort(TH);
    PREF     = PREF(:,xxx);
    NULL     = NULL(:,xxx);
    FANO     = FANO(:,xxx);


    %%
    % get average weights for MT neurons across simulation
    wavg = nans(length(PREF), NTUNE, length(trialnum));
    for i = 1:length(trialnum)
        for j = 1:length(PREF)
            w_ = [];
            for k = 1:length(SIMNUM)
                L = ind(:,k,i)==j;
                w_          = [w_; w(L,:,k,i)];
            end
            wavg(j,:,i) = nanmean(w_,1);
        end
    end
            
    
    %% assign average weights back to simulations
    wm = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
    for i = 1:length(SIMNUM)
        for j = 1:length(trialnum)
            wm(:,:,i,j) = wavg(ind(:,i,j),:,j);
        end
    end    
end




% get weights for early/mid/late trial, for each simulation
cp    = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
MTind = nans(NSEN, NTUNE, length(SIMNUM), length(trialnum));
if recompute
    for i = 1:length(SIMNUM)
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(eln_) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);
        clear FIRA
        pack

        for j = 1:length(trialnum)
            cp(:,:,i,j)  = getML_RLCPForWeight3(Monk, NSEN, NTUNE, wm(:,:,i,j), ind(:,i,j), 2000, sprintf('%d-%d', i, j), 1);
            MTind(:,:,i,j) = repmat(MTindex',1,NTUNE);
        end
    end
    th = TH(MTind);
else
    for i = 1:length(SIMNUM)
        [hdir, ldir, cdir, tdir] = dirnames;
        fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(eln_) '_v' num2str(SIMNUM(i)) '.mat'];
        load([tdir fname])
        fprintf([fname '\n']);
        clear FIRA
        pack

        for j = 1:length(trialnum)
            cp(:,:,i,j)  = getML_RLCPForWeight3(Monk, NSEN, NTUNE, [], [], 2000, sprintf('%d-%d', i, j), 0);
            MTind(:,:,i,j) = repmat(MTindex',1,NTUNE);
        end
    end
end

[fitsp, semsp] = getML_neurometricROCT([Monk 'PRe_MT.txt'], [], 1, [0 1], 0, 0);
[fits, sems]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], [], 1, [0 1], 0, 0);
[pp, np, fp  ] = getML_SigNoise2([Monk 'PRe_MT.txt'], 1, [0 1], 1, 0);
[p, n, f]      = getML_SigNoise2([Monk 'TRain_MT.txt'], 1, [0 1], 1, 0);

% combine pre- and during training data
PREF = [pp p];
NULL = [np n];
FANO = [fp f];
TH   = getML_RLSimulatedThForMT(Monk, 10000, 0);
if strcmp(Monk, 'Cy')
    Lgd      = ~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
    Lgd(178) = logical(0);  % this cell is weird, with very small gain to both pref and null
else
    Lgd      =~isnan(PREF(1,:)) & ~isnan(NULL(1,:)) & ~isnan(FANO(1,:));
end
PREF = PREF(:,Lgd);
NULL = NULL(:,Lgd);
FANO = FANO(:,Lgd);
TH   = TH(Lgd);

% sort data by threshold coherence
[TH,xxx] = sort(TH);
PREF     = PREF(:,xxx);
NULL     = NULL(:,xxx);
FANO     = FANO(:,xxx);

th = TH(MTind);


% binned cp by thresholds
bb     = 2.32;
be     = 6.53;
bw     = (6.53-2.32)/50;
THBINS = [2.^[bb:bw:be-bw]' 2.^[bb+bw:bw:be]'];
cpm    = nans(length(THBINS), NTUNE, length(trialnum));
cpd    = nans(length(THBINS), NTUNE, length(trialnum));
cpn    = nans(length(THBINS), NTUNE, length(trialnum));
   
% average weight for each threshold bin
for i = 1:length(trialnum)
    for j = 1:NTUNE
        cp_  = cp(:,j,:,i);
        cp_  = cp_(:);
        th_  = 100*th(:,j,:,i);
        th_  = th_(:);
        
        for k = 1:length(THBINS)
            L          = th_>=THBINS(k,1) & th_<THBINS(k,2);
            cpm(k,j,i) = nanmean(cp_(L));
            cpd(k,j,i) = nanse(cp_(L));
            cpn(k,j,i) = sum(L);
            
        end
    end
end



%% Figure 4c: correlation structure of an example cell with ranked
%  sensitivity of 100, and direction tuning of 0 degree
th = annotation('textbox', [0 0.48 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


asize = 0.2;

rankc      = 100;
tunec      = 0;

ranksen    = repmat([1:200]',1,NTUNE);
DTUNE      = 360/NTUNE;
dirtune    = repmat(linspace(-170,180,NTUNE),NSEN,1);

rsen       = 0.5;       % max correlation between neurons with different sensitivity
rslp       = -rsen/20;  % number of neighbouring neurons that it is correlated to 
gsen       = rsen+rslp*abs(ranksen-rankc);
gsen(gsen<0) = 0;
gsen       = gsen/max(max(gsen));

gtune      = exp(-(abs(dirtune))/30);

corrstruct = rsen.*gsen.*gtune;




% plot
clim_ = [0 rsen];
DTUNE = 360/NTUNE;
pdir  = linspace(-170,180,NTUNE);

% plot weights
set(axes, 'Units', 'Normalized', 'Position', [0.08 0.1 asize/1.3 asize*pw/ph], 'FontName', 'Helvetica')
imagesc(pdir, 1:200, corrstruct)
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'XTick', [-170 -135 -90 -45 0 45 90 135 180], 'TickDir', 'out', ...
         'XTickLabel', {'-170', '', '', '', '0', '', '', '', '180'}, ...
         'YTick', [1 50:50:200], 'clim', clim_)
xlabel('Direction tuning (degree)', 'FontSize', lfs)
ylabel('Ranked sensitivity', 'FontSize', lfs)
%title('Correlation structure', 'FontSize', tfs)



%% Figure 4d: Choice probablility early mid and late in training, with
%             decaying correlation structure
th = annotation('textbox', [0.26 0.48 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'd', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 

asize = 0.2;

DTUNE  = 360/NTUNE;
pdir   = linspace(-180,180-DTUNE,NTUNE);
thm    = nangeomean(THBINS,2);
ind    = find(pdir==0); 

bdata  = [0.5248 -0.0064; 0.6153 -0.0290; 0.6674 -0.0397];  % slope of MT CP from real data, value from LAW GOLD 2007 Figure 3


set(axes, 'Units', 'Normalized', 'Position', [0.125+1*(0.02+asize) 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([thm thm]',[cpm(:,ind,1)-cpd(:,ind,1), cpm(:,ind,1)+cpd(:,ind,1)]','-', 'Color', lc1e, 'LineWidth', elw)
plot(thm,cpm(:,ind,1),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd      = ~isnan(cpm(:,ind,1)) & cpd(:,ind,1)~=0;
[b, bi, h, p, F, df] = nestedFW(cpm(Lgd,ind,1), 1./cpn(Lgd,ind,1), [ones(size(log(thm(find(Lgd))))) log(thm(find(Lgd)))]);
plot(2:100, b(1)+b(2)*log(2:100), 'k-', 'LineWidth', lw)
%plot(2:100, bdata(1,1)+bdata(1,2)*log(2:100), 'r-', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickDir', 'out', ...
         'xscale', 'log', 'xlim', [4.5 100], 'ylim', [0.49 0.545], ...
         'XTick', [5 10 20 40 80], 'YTick', [0.5 0.52 0.54])
xlabel('Threshold (% coh)', 'FontSize', lfs)
ylabel('Choice Probability', 'FontSize', lfs)
%title('Early', 'FontSize', tfs)
  


set(axes, 'Units', 'Normalized', 'Position', [0.125+2*(0.02+asize) 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([thm thm]',[cpm(:,ind,2)-cpd(:,ind,2), cpm(:,ind,2)+cpd(:,ind,2)]','-', 'Color', lc1e, 'LineWidth', elw)
plot(thm,cpm(:,ind,2),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd      = ~isnan(cpm(:,ind,2)) & cpd(:,ind,2)~=0;
[b, bi, h, p, F, df] = nestedFW(cpm(Lgd,ind,2), 1./cpn(Lgd,ind,2), [ones(size(log(thm(find(Lgd))))) log(thm(find(Lgd)))]);
plot(2:100, b(1)+b(2)*log(2:100), 'k-', 'LineWidth', lw)
%plot(2:100, bdata(2,1)+bdata(2,2)*log(2:100), 'r-', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickDir', 'out', ...
         'xscale', 'log', 'xlim', [4.5 100], 'ylim', [0.49 0.545], ...
         'XTick', [5 10 20 40 80], 'YTick', [0.5 0.52 0.54], 'YTickLabel', {})
xlabel('Threshold (% coh)', 'FontSize', lfs)
%title('Mid', 'FontSize', tfs)
  
   


set(axes, 'Units', 'Normalized', 'Position', [0.125+3*(0.02+asize) 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([thm thm]',[cpm(:,ind,3)-cpd(:,ind,3), cpm(:,ind,3)+cpd(:,ind,3)]','-', 'Color', lc1e, 'LineWidth', elw)
plot(thm,cpm(:,ind,3),'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', 'none', 'MarkerSize', ms)
Lgd      = ~isnan(cpm(:,ind,3)) & cpd(:,ind,3)~=0;
[b, bi, h, p, F, df] = nestedFW(cpm(Lgd,ind,3), 1./cpn(Lgd,ind,3), [ones(size(log(thm(find(Lgd))))) log(thm(find(Lgd)))]);
plot(2:100, b(1)+b(2)*log(2:100), 'k-', 'LineWidth', lw)
%plot(2:100, bdata(3,1)+bdata(3,2)*log(2:100), 'r-', 'LineWidth', lw)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickDir', 'out', ...
         'xscale', 'log', 'xlim', [4.5 100], 'ylim', [0.49 0.545], ...
         'XTick', [5 10 20 40 80], 'YTick', [0.5 0.52 0.54], 'YTickLabel', {})
xlabel('Threshold (% coh)', 'FontSize', lfs)
%title('Late', 'FontSize', tfs)
  
   

