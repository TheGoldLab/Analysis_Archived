% test how correlation affects the time coarses of learning, asymptotic
% performance amd choice probability

%% get data
%% test correlations on asymptotic thresholds 
% rmax
Monk = 'Cy';
DIRS = 0;
COHS = [0 3.2 6.4 12.8 25.6 51.2 99.9]'/100;
a_rmax = [0:0.1:0.9];
bdir = 30;
bsen = 20;
N    = 500;
recompute = 0;

p = getML_RLFigS3CorrAsym(Monk, DIRS, COHS, a_rmax, bdir, bsen, N, recompute);

a_rm  = nans(2,length(a_rmax));
ad_rm = nans(2,length(a_rmax));

for i = 1:length(a_rmax)
    [a_rm(:,i), ad_rm(:,i)] = ctPsych_fit(@quick2,COHS, [p(:,i) N*ones(size(p(:,i)))]);
end



% bdir
Monk = 'Cy';
DIRS = 0;
COHS = [0 3.2 6.4 12.8 25.6 51.2 99.9]'/100;
rmax = 0.5;
a_bdir = [10:20:170];
bsen = 20;
N    = 500;
recompute = 0;

p = getML_RLFigS3CorrAsym(Monk, DIRS, COHS, rmax, a_bdir, bsen, N, recompute);
p = squeeze(p);

a_bd  = nans(2,length(a_bdir));
ad_bd = nans(2,length(a_bdir));

for i = 1:length(a_bdir)
    [a_bd(:,i), ad_bd(:,i)] = ctPsych_fit(@quick2,COHS, [p(:,i) N*ones(size(p(:,i)))]);
end



% bsen
Monk = 'Cy';
DIRS = 0;
COHS = [0 3.2 6.4 12.8 25.6 51.2 99.9]'/100;
rmax = 0.5;
bdir = 30;
a_bsen = [10:20:170];
N    = 500;
recompute = 0;

p = getML_RLFigS3CorrAsym(Monk, DIRS, COHS, rmax, bdir, a_bsen, N, recompute);
p = squeeze(p);

a_bs  = nans(2,length(a_bsen));
ad_bs = nans(2,length(a_bsen));

for i = 1:length(a_bsen)
    [a_bs(:,i), ad_bs(:,i)] = ctPsych_fit(@quick2,COHS, [p(:,i) N*ones(size(p(:,i)))]);
end




%% test correlation on choice probabilities
% rmax
Monk = 'Cy';
DIRS = 0;
COHS = [0 3.2 6.4 12.8 25.6 51.2 99.9]'/100;
cp_rmax = [0:0.15:0.9];
bdir = 30;
bsen = 20;
N    = 500;
recompute = 0;

[cp_rm, cpd_rm] = getML_RLFigS3CorrCP(Monk, DIRS, cp_rmax, bdir, bsen, N, recompute);



% bdir
Monk = 'Cy';
DIRS = 0;
COHS = [0 3.2 6.4 12.8 25.6 51.2 99.9]'/100;
rmax = 0.5;
cp_bdir = [10 30 50 80 110 140 170];
bsen = 20;
N    = 500;
recompute = 0;

[cp_bd, cpd_bd] = getML_RLFigS3CorrCP(Monk, DIRS, rmax, cp_bdir, bsen, N, recompute);



% bsen
Monk = 'Cy';
DIRS = 0;
COHS = [0 3.2 6.4 12.8 25.6 51.2 99.9]'/100;
rmax = 0.5;
bdir = 30;
cp_bsen = [10 20 40 60 90 120 160];
N    = 500;
recompute = 0;

[cp_bs, cpd_bs] = getML_RLFigS3CorrCP(Monk, DIRS, rmax, bdir, cp_bsen, N, recompute);
cp_bs  = squeeze(cp_bs);
cpd_bs = squeeze(cpd_bs);



%% test correlations on time coarses of learning
% rmax
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
t_rmax    = [0.1:0.2:0.9];
bdir      = 30;
bsen      = 20;
RECOMPUTE = 0;

% for rr1 = 1:length(t_rmax)
%     for rr2 = 1:length(bdir)
%         for rr3 = 1:length(bsen)
%             time1            = clock;
%             [FIRA, W, MTindex] = getML_RLSimPerceptualLearningTestCorr(Monk, NSEN, NTUNE, eln_, DIRS, SIMNUM, t_rmax(rr1), bdir(rr2), bsen(rr3), RECOMPUTE);
%             time2            = clock;
%             time2-time1
%             %clear FIRA W MTind
%             pack
%         end
%     end
% end

[bt_rm, btd_rm, bl_rm, bld_rm] = getML_RLFigS3Tau(Monk, NSEN, NTUNE, eln_, DIRS, SIMNUM, t_rmax, bdir, bsen, RECOMPUTE);


% bdir
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
rmax      = 0.5;
t_bdir    = [15 30 60 100 150];
bsen      = 20;
RECOMPUTE = 0;

% for rr1 = 1:length(rmax)
%     for rr2 = 1:length(t_bdir)
%         for rr3 = 1:length(bsen)
%             time1            = clock;
%             [FIRA, W, MTindex] = getML_RLSimPerceptualLearningTestCorr(Monk, NSEN, NTUNE, eln_, DIRS, SIMNUM, rmax(rr1), t_bdir(rr2), bsen(rr3), RECOMPUTE);
%             time2            = clock;
%             time2-time1
%             %clear FIRA W MTind
%             pack
%         end
%     end
% end
[bt_bd, btd_bd, bl_bd, bld_bd] = getML_RLFigS3Tau(Monk, NSEN, NTUNE, eln_, DIRS, SIMNUM, rmax, t_bdir, bsen, RECOMPUTE);


% bsen
Monk      = 'Cy';
NSEN      = 200;
NTUNE     = 36;
eln_      = [7]*1e-8; %for correlated response
DIRS      = [-90 90];
SIMNUM    = [81];
rmax      = 0.5;
bdir      = 30;
t_bsen    = [10 20 50 80 120];
RECOMPUTE = 0;

% for rr1 = 1:length(rmax)
%     for rr2 = 1:length(bdir)
%         for rr3 = 1:length(t_bsen)
%             time1            = clock;
%             [FIRA, W, MTindex] = getML_RLSimPerceptualLearningTestCorr(Monk, NSEN, NTUNE, eln_, DIRS, SIMNUM, rmax(rr1), bdir(rr2), t_bsen(rr3), RECOMPUTE);
%             time2            = clock;
%             time2-time1
%             %clear FIRA W MTind
%             pack
%         end
%     end
% end

[bt_bs, btd_bs, bl_bs, bld_bs] = getML_RLFigS3Tau(Monk, NSEN, NTUNE, eln_, DIRS, SIMNUM, rmax, bdir, t_bsen, RECOMPUTE);



%% plot
% create figure
fh = figure; 
mm2inch = 25.4;
pw = 5;  % paper width, make it nature size of one column, 85mm x 85mm
ph = 8;   % paper height
wysifig(fh, pw, ph) % set it to US letter size

% define figure parameters
alw = 0.5;  % axis line width
afs = 8;    % axis font size
lfs = 8;    % label font size
ffs = 12;   % figure font size
tfs = 18;    % text font size
lw  = 1.5;  % line width
elw = 0.5;
ms  = 5;    % marker size

% define line color
lc1  = [0.2 0.2 0.2];
lc1e = [0.75 0.75 0.75];
lc2  = [1 0.2 0.2];
lc2e = [1 0.75 0.75];


asize = 0.22;



%% tau lapse
set(axes, 'Units', 'Normalized', 'Position', [0.1+0*(asize+0.08) 0.1+3*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([t_rmax; t_rmax], shiftdim(bld_rm(3,:,:),1), '-k', 'Color', lc1e)
plot(t_rmax, bl_rm(3,:), 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 1], 'ylim', [0 1000], ...
         'XTick', [0 0.5 1], 'XTickLabel', {'', '', ''}, ...  
         'YTick', [0 500 1000], 'YTickLabel', {'0', '5', '1'}) 
[xb,yb] = axisXY2figXY(0.5,133);
[xe,ye] = axisXY2figXY(0.5,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
ylabel('\tau_{la}, trials (x1000)', 'FontSize', lfs)
%title('rho_{max}', 'FontSize', lfs)      
      


set(axes, 'Units', 'Normalized', 'Position', [0.1+1*(asize+0.08) 0.1+3*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([t_bdir; t_bdir], shiftdim(bld_bd(3,:,:),1), '-k', 'Color', lc1e)
plot(t_bdir, bl_bd(3,:), 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 180], 'ylim', [0 1000], ...
         'XTick', [0 90 180], 'XTickLabel', {'', '', ''}, ...  
         'YTick', [0 500 1000], 'YTickLabel', {'', '', ''}) 
[xb,yb] = axisXY2figXY(30,133);
[xe,ye] = axisXY2figXY(30,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)



set(axes, 'Units', 'Normalized', 'Position', [0.1+2*(asize+0.08) 0.1+3*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([t_bsen; t_bsen], shiftdim(bld_bs(3,:,:),1), '-k', 'Color', lc1e)
plot(t_bsen, bl_bs(3,:), 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 180], 'ylim', [0 1000], ...
         'XTick', [0 90 180], 'XTickLabel', {'', '', ''}, ...  
         'YTick', [0 500 1000], 'YTickLabel', {'', '', ''}) 
[xb,yb] = axisXY2figXY(20,133);
[xe,ye] = axisXY2figXY(20,0);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
%title('b_{sen}', 'FontSize', lfs)


%% tau th
set(axes, 'Units', 'Normalized', 'Position', [0.1+0*(asize+0.08) 0.1+2*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([t_rmax; t_rmax], shiftdim(btd_rm(3,:,:),1), '-k', 'Color', lc1e)
plot(t_rmax, bt_rm(3,:), 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 1], 'ylim', [15000 25000], ...
         'XTick', [0 0.5 1], 'XTickLabel', {'', '', ''}, ... 
         'YTick', [15000 20000 25000], 'YTickLabel', {'15', '20', '25'}) 
[xb,yb] = axisXY2figXY(0.5,16333);
[xe,ye] = axisXY2figXY(0.5,15000);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
ylabel('\tau_{th}, trials (x1000)', 'FontSize', lfs)


set(axes, 'Units', 'Normalized', 'Position', [0.1+1*(asize+0.08) 0.1+2*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([t_bdir; t_bdir], shiftdim(btd_bd(3,:,:),1), '-k', 'Color', lc1e)
plot(t_bdir, bt_bd(3,:), 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 180], 'ylim', [15000 25000], ...
         'XTick', [0 90 180], 'XTickLabel', {'', '', ''}, ...  
         'YTick', [15000 20000 25000], 'YTickLabel', {'', '', ''}) 
[xb,yb] = axisXY2figXY(30,16333);
[xe,ye] = axisXY2figXY(30,15000);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)



set(axes, 'Units', 'Normalized', 'Position', [0.1+2*(asize+0.08) 0.1+2*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([t_bsen; t_bsen], shiftdim(btd_bs(3,:,:),1), '-k', 'Color', lc1e)
plot(t_bsen, bt_bs(3,:), 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 180], 'ylim', [15000 25000], ...
         'XTick', [0 90 180], 'XTickLabel', {'', '', ''}, ...  
         'YTick', [15000 20000 25000], 'YTickLabel', {'', '', ''}) 
[xb,yb] = axisXY2figXY(20,16333);
[xe,ye] = axisXY2figXY(20,15000);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)


%% asym
set(axes, 'Units', 'Normalized', 'Position', [0.1+0*(asize+0.08) 0.1+1*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([a_rmax; a_rmax], [a_rm(1,:)-ad_rm(1,:); a_rm(1,:)+ad_rm(1,:)], '-k', 'Color', lc1e)
plot(a_rmax, a_rm(1,:), 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 1], 'ylim', [0.06 0.12], ...
         'XTick', [0 0.5 1], 'XTickLabel', {'', '', ''}, ... 
         'YTick', [0.06 0.09 0.12], 'YTickLabel', {'6', '9', '12'}) 
[xb,yb] = axisXY2figXY(0.5,0.068);
[xe,ye] = axisXY2figXY(0.5,0.06);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
ylabel('% coh', 'FontSize', lfs)

      
set(axes, 'Units', 'Normalized', 'Position', [0.1+1*(asize+0.08) 0.1+1*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([a_bdir; a_bdir], [a_bd(1,:)-ad_bd(1,:); a_bd(1,:)+ad_bd(1,:)], '-k', 'Color', lc1e)
plot(a_bdir, a_bd(1,:), 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 180], 'ylim', [0.06 0.12], ...
         'XTick', [0 90 180], 'XTickLabel', {'', '', ''}, ... 
         'YTick', [0.06 0.09 0.12], 'YTickLabel', {'', '', ''}) 
[xb,yb] = axisXY2figXY(30,0.068);
[xe,ye] = axisXY2figXY(30,0.06);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)


set(axes, 'Units', 'Normalized', 'Position', [0.1+2*(asize+0.08) 0.1+1*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([a_bsen; a_bsen], [a_bs(1,:)-ad_bs(1,:); a_bs(1,:)+ad_bs(1,:)], '-k', 'Color', lc1e)
plot(a_bsen, a_bs(1,:), 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 180], 'ylim', [0.06 0.12], ...
         'XTick', [0 90 180], 'XTickLabel', {'', '', ''}, ...
         'YTick', [0.06 0.09 0.12], 'YTickLabel', {'', '', ''}) 
[xb,yb] = axisXY2figXY(20,0.068);
[xe,ye] = axisXY2figXY(20,0.06);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)





%% cp
set(axes, 'Units', 'Normalized', 'Position', [0.1+0*(asize+0.08) 0.1+0*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([cp_rmax; cp_rmax], [cp_rm-cpd_rm cp_rm+cpd_rm]', '-k', 'Color', lc1e)
plot(cp_rmax, cp_rm, 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 1], 'ylim', [0.45 0.65], ...
         'XTick', [0 0.5 1], ...  
         'YTick', [0.45 0.55 0.65], 'YTickLabel', {'0.45', '0.55', '0.65'}) 
[xb,yb] = axisXY2figXY(0.5,0.4767);
[xe,ye] = axisXY2figXY(0.5,0.45);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
ylabel('Fraction correct', 'FontSize', lfs)
xlabel('\rho_{max}', 'FontSize', lfs)



set(axes, 'Units', 'Normalized', 'Position', [0.1+1*(asize+0.08) 0.1+0*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([cp_bdir; cp_bdir], [cp_bd-cpd_bd; cp_bd+cpd_bd], '-k', 'Color', lc1e)
plot(cp_bdir, cp_bd, 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 180], 'ylim', [0.45 0.65], ...
         'XTick', [0 90 180], ...
         'YTick', [0.45 0.55 0.65], 'YTickLabel', {'', '', ''}) 
[xb,yb] = axisXY2figXY(30,0.4767);
[xe,ye] = axisXY2figXY(30,0.45);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
xlabel('b_{dir} (deg)', 'FontSize', lfs)



set(axes, 'Units', 'Normalized', 'Position', [0.1+2*(asize+0.08) 0.1+0*(asize+0.08)*pw/ph asize asize*pw/ph], 'FontName', 'Helvetica')
hold on
plot([cp_bsen; cp_bsen], [cp_bs-cpd_bs cp_bs+cpd_bs]', '-k', 'Color', lc1e)
plot(cp_bsen, cp_bs, 'ok', 'MarkerSize', ms, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
hold off
set(gca, 'box', 'off', 'LineWidth', alw, 'FontSize', afs, 'TickLength', [0.03 0.025], ...
         'xlim', [0 180], 'ylim', [0.45 0.65], ...
         'XTick', [0 90 180], ...
         'YTick', [0.45 0.55 0.65], 'YTickLabel', {'', '', ''})
[xb,yb] = axisXY2figXY(20,0.4767);
[xe,ye] = axisXY2figXY(20,0.45);
anh     = annotation('arrow', [xb xe], [yb ye]);
set(anh, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-', 'HeadStyle', 'cback2', ...
          'HeadWidth', ms, 'HeadLength', ms)
xlabel('b_{sen}', 'FontSize', lfs)


