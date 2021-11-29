% plot supplemental figure 1

% Matlab script for plotting Fig S1  (Model parameters)
% for submission to Nature Neuroscience
%
% created by jcl on 08/07/08

%%
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 170/mm2inch;    % paper width, make it nature size of one column, 85mm x 85mm
ph = 130/mm2inch;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 8;    % axis font size
lfs = 10;    % label font size
ffs = 12;   % figure font size
tfs = 12;    % text font size
lw  = 1;  % line width
elw = 0.5;
ms  = 5;    % marker size

% define line color
lc1  = [0.2 0.2 0.2];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.2 0.2];
lc2e = [1 0.75 0.75];


asize = 0.22;




%% get mean MT firing rate as a function of coh, and fano factor
IND   = 24;   % use cell 24 for illustration
Monk = 'Cy';
a      = getML_txt([Monk 'PRe_MT.txt']);
fn     = a.data{strcmp(a.name,'dat_fn')};
usable = a.data{strcmp(a.name,'usable')};
uid    = a.data{strcmp(a.name,'uid')};
d1     = a.data{strcmp(a.name,'ddir')};

rM     = nans(14, 1);
rN     = nans(14, 1);
rSD    = nans(14, 1);

pref   = nans(2,1);
null   = nans(2,1);
fano   = nans(1,1);
 

% get rate
openFIRA(fn{IND})
fprintf('%s\n',fn{IND})

% get mean rate from dots_on to dots_off
% selection arrays for coherences
[Lcoh, Ucoh] = selectFIRA_trialsByUniqueID('dot_coh');
Lcoh = Lcoh(:, ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9])); % remove strange coh other than the standard cohs
Ucoh = Ucoh(ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9]));    % in one pre-training file there's a 9% coh condition

% selection arrays for direction
[Ldir, Udir] = selectFIRA_trialsByUniqueID('dot_dir');
Idir         = [find(round(Udir)==d1(IND)) find(round(Udir)==mod(d1(IND)+180,360))];

% selection arrays for task
[Ltk, Utk]  = selectFIRA_trialsByUniqueID('task');
Ltk = ismember(getFIRA_ecodesByName('task','id'), [3,6]);

% selection array for correct/incorrect
Lcrt =  ismember(getFIRA_ecodesByName('correct', 'id'), 1);

% selection arry for viewing time
THRESH = 200;   % minimum viewing time
vt     = getFIRA_ecodeTimesByName('dot_off', 0)-getFIRA_ecodeTimesByName('dot_on', 0);
Ltime  = vt>=THRESH;

% get selection array for mean rates
L = zeros(length(Ldir), 14);
c = [0 3.2 6.4 12.8 25.6 51.2 99.9];
for j = 1:7
    if ismember(c(j),Ucoh)
        L(:,j)   = Ldir(:,Idir(1)) & Lcoh(:,Ucoh==c(j)) & Ltk & Lcrt & Ltime;
        L(:,j+7) = Ldir(:,Idir(2)) & Lcoh(:,Ucoh==c(j)) & Ltk & Lcrt & Ltime;
    end
end

% get rate
for j = 1:14
    r = getFIRA_rate(L(:,j), getFIRA_spikeByID(uid(IND)), ...
        getFIRA_ecodeTimesByName('dot_on', 0), ...
        getFIRA_ecodeTimesByName('dot_off', 0));
    rM(j)  = nanmean(r);
    rN(j)  = sum(~isnan(r));
    rSD(j) = nanstd(r);
end
       

         
% get slope to pref/null and response at 0% coh
mp    = rM(1:7);
sep   = rSD(1:7)%./(rN(1:7).^0.5);
coh   = [0 0.032 0.064 0.128 0.256 0.512 0.999]';

Lgd = ~isnan(mp);
if sum(Lgd)>2
    [b be]   = regressW(mp(Lgd), sep(Lgd), [ones(size(coh(Lgd))) coh(Lgd)]);
    pref(:)  = b;
end

mn     = rM(8:14);
sen    = rSD(8:14)%./(rN(8:14).^0.5);
coh    = [0 0.032 0.064 0.128 0.256 0.512 0.999]';

Lgd = ~isnan(mn);
if sum(Lgd)>2
    [b be]   = regressW(mn(Lgd), sen(Lgd), [ones(size(coh(Lgd))) coh(Lgd)]);
    null(:) = b;
end



% get fano
mf        = rM;
vf        = rSD.^2;
Lgd       = ~isnan(rM);
[b be]    = regressW(vf(Lgd), 1./(rN(Lgd).^0.5), mf(Lgd));
fano(:)   = b;
           

%% plot example MT fit
th = annotation('textbox', [0 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'a', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% plot Cy
coh = 100*[0 0.032 0.064 0.128 0.256 0.512 0.999]';
set(axes, 'Units', 'Normalized', 'Position', [0.08 0.61 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([coh coh]',[mp-sep mp+sep]','-', 'Color', lc1e, 'LineWidth', elw)
plot([coh coh]',[mn-sen mn+sen]','-', 'Color', lc1e, 'LineWidth', elw)
plot(coh,mp,'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', lc1, 'MarkerSize', ms)
plot(coh,mn,'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc1, 'MarkerSize', ms)
plot(coh, null(1)+coh/100*null(2), 'LineStyle', '--', 'LineWidth', lw, 'Color', 'k')
plot(coh, pref(1)+coh/100*pref(2), 'LineStyle', '-', 'LineWidth', lw, 'Color', 'k')
plot(coh, nanmean([pref(1) null(1)])*ones(size(coh)), 'LineStyle', ':', 'LineWidth', lw, 'Color', lc1e)
hold off

set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 100], 'ylim', [0 80], ...
         'XTick', [0 50 100], ...  
         'YTick', [0 40 80]) 
xlabel('% coherence', 'FontSize', lfs)
ylabel('Response (spikes per sec)', 'FontSize', lfs)
text(32, 35, 'k_p', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', tfs)
text(32, -1, 'k_n', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', tfs)
text(90, nanmean([pref(1) null(1)]), 'k_0', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', tfs)



%% plot mean to variance ratio and fano factor
th = annotation('textbox', [0.32 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'b', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% plot Cy
coh = 100*[0 0.032 0.064 0.128 0.256 0.512 0.999]';
set(axes, 'Units', 'Normalized', 'Position', [0.42 0.61 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(mf, vf,'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', lc1, 'MarkerSize', ms)
plot(rM, rM*fano, 'LineStyle', '-', 'LineWidth', lw, 'Color', 'k')
hold off

set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 70], 'ylim', [0 150], ...
         'XTick', [0 35 70], ...  
         'YTick', [0 75 150]) 
xlabel('Mean (spikes per sec)', 'FontSize', lfs)
ylabel('Variance (spikes^2 per sec^2)', 'FontSize', lfs)
text(52, 90, '\phi', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', tfs+2)



%% plot mean response as a function of motion direction
th = annotation('textbox', [0.66 0.95 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'c', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% plot Cy
dirs  = -180:10:170;
atune = normpdf(dirs, 90, 40);
atune = atune/max(max(atune));              % tuning modulation function
acoh  = (pref(2)-null(2)).*atune+null(2);   % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)
k0    = nanmean([pref(1) null(1)]);
        
set(axes, 'Units', 'Normalized', 'Position', [0.76 0.61 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(dirs, k0+0.064*acoh, 'LineStyle', '-', 'LineWidth', lw+0.5, 'Color', [0.8 0.8 0.8])
plot(dirs, k0+0.256*acoh, 'LineStyle', '-', 'LineWidth', lw+0.5, 'Color', [0.5 0.5 0.5])
plot(dirs, k0+0.999*acoh, 'LineStyle', '-', 'LineWidth', lw+0.5, 'Color', [0 0 0])
plot(60:120, (k0+0.999*pref(2))*ones(size([60:120])), 'LineStyle', ':', 'LineWidth', lw, 'Color', 'k')
plot(-90:5, (k0+0.999*null(2))*ones(size([-90:5])), 'LineStyle', ':', 'LineWidth', lw, 'Color', 'k')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [-180 170], 'ylim', [0 70], ...
         'XTick', [-180 -90 0 90 170], ...  
         'YTick', [0 35 70]) 
xlabel('Motion direction (degree)', 'FontSize', lfs)
ylabel('Response (spikes per sec)', 'FontSize', lfs)
% plot annotation arrows 
[xb,yb] = axisXY2figXY(10,k0+0.999*pref(2));
[xe,ye] = axisXY2figXY(42,k0+0.999*pref(2));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5, 'HeadLength', 1.5)
% plot annotation arrows 
[xb,yb] = axisXY2figXY(10,k0+0.999*pref(2));
[xe,ye] = axisXY2figXY(38,k0+0.999*pref(2));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5, 'HeadLength', 1.5)
[xb,yb] = axisXY2figXY(48,k0+0.999*null(2));
[xe,ye] = axisXY2figXY(20,k0+0.999*null(2));
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 1.5, 'HeadLength', 1.5)
% text
text(0, k0+0.999*pref(2)-2, 'k_0+k_p', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', tfs-2)
text(60, k0+0.999*null(2)-2, 'k_0+k_n', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', tfs-2)





%% get cp for random initial weights for different additive decision noise levels
Monk   = 'Cy';
NSEN   = 200;
NTUNE  = 36;
N      = 1000;
DNOISE = [0:1:10];
recompute  = 0;

Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8;
DIRS      = [-90,90];
SIMNUM    = [81];


trialnum  = floor([2000]/1000)+1;
[hdir, ldir, cdir, tdir] = dirnames;
fname = ['/getML_RL_' Monk '_' int2str(NSEN) '_' int2str(NTUNE) '_' num2str(eln_) '_v' num2str(SIMNUM) '.mat'];
load([tdir fname])
fprintf([fname '\n']);

[cpm cpstd] = getML_RLFigS1AverageCP(Monk, NSEN, NTUNE, N, W(:,:,trialnum(1)), DNOISE, recompute);


%% plot cp for neurons tuned to 90 degree as a function of decision noise
th = annotation('textbox', [0 0.42 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'd', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% plot Cy
set(axes, 'Units', 'Normalized', 'Position', [0.08 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(DNOISE, 0.5*ones(size(DNOISE)),':', 'Color', 'k', 'LineWidth', elw)
plot([DNOISE; DNOISE], [cpm-cpstd./sqrt(200); cpm+cpstd./sqrt(200)], '-', 'Color', lc1e, 'LineWidth', elw)
plot(DNOISE, cpm,'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', lc1, 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 8], 'ylim', [0.49 0.53], ...
         'XTick', [0 2 4 6 8], ...  
         'YTick', [0 0.5 0.515 0.53]) 
xlabel('SD of decision noise', 'FontSize', lfs)
ylabel('Average choice probability', 'FontSize', lfs)

% plot annotation arrows 
[xb,yb] = axisXY2figXY(5,0.511);
[xe,ye] = axisXY2figXY(5,0.506);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 2, 'HeadLength', 2)

      
      
      
      
%% mean of pooled signal with optimal weight as a function of normalization factor 
NORMSCALE = 0.0001*2.^[1:0.5:8];
NORMSCALE = [NORMSCALE 0.02];
recompute = 0;
[mpool]   = getML_RLFigS1WNormalization(NORMSCALE, recompute);      


th = annotation('textbox', [0.32 0.42 0.05 0.05]);
set(th, 'LineStyle', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'String', 'e', 'FontSize', ffs, 'FontWeight', 'bold', 'FontName', 'Helvetica') 


% plot Cy
set(axes, 'Units', 'Normalized', 'Position', [0.42 0.1 asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot(0:0.005:0.02, mpool(end)*ones(size(0:0.005:0.02)),':', 'Color', 'k', 'LineWidth', elw)
plot(0.02*ones(size(linspace(0,mpool(end),100))), linspace(0,mpool(end),100),':', 'Color', 'k', 'LineWidth', elw)
plot(NORMSCALE, mpool, 'o', 'MarkerFaceColor', lc1, 'MarkerEdgeColor', lc1, 'MarkerSize', ms)
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, ...
         'xlim', [0 0.03], 'ylim', [0 150], ...
         'XTick', [0:0.01:0.03], ...  
         'YTick', [0 50 100 150]) 
xlabel('Weight normalization', 'FontSize', lfs)
ylabel('Response (spikes per sec)', 'FontSize', lfs)
% plot annotation arrows 
[xb,yb] = axisXY2figXY(0.02,130);
[xe,ye] = axisXY2figXY(0.02,110);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', 2, 'HeadLength', 2)
































