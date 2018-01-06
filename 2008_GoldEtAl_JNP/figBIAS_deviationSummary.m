function fig_ = figBIAS_deviationSummary(num)
% function fig_ = figBIAS_deviationSummary(num)
%
if nargin < 1 || isempty(num)
    num = 11;
end

% get monk names
[monks,monkn,mse] = getBIAS_monks('FEF');

%% Set up Fig
% units should be in inches, from wysifig
wid  = 5.8; % total width
hts  = 1.3;
cols = {3,3};
[axs,fig_] = getBIAS_axes(num, wid, hts, cols);
set(axs,'Units', 'Normalized');

% get the data
%pmfdat = FS_loadProjectFile('2008_Bias', 'figBIAS_pmfParams');
fcs    = getBIAS_fcs;
pmfdat = FS_loadProjectFile('2008_Bias', 'figBIAS_pmfCompare');
devdat = FS_loadProjectFile('2008_Bias', 'figBIAS_deviationSummary');

if isempty(devdat)

    dind   = 20;
    devdat = cell(monkn,2);

    for mm = 1:monkn

        dat          = FS_getDotsTrainingData(monks{mm});
        sessions     = unique(dat(:,1));
        num_sessions = length(sessions);
        
        bci          =  fcs{mm}.*dat(:,9); % 'biased-choice index': wk-filt x choice
        ts           = (max(0.1, dat(:,6))-0.1).^0.5;
        xdat         = [bci dat(:,dind)];
        ydat         = [dat(:,5).*ts, ts];

        Learly       = dat(:,1) < sessions(round(num_sessions/4));
        Lgood        = dat(:,2) < 3 & dat(:,3) >= 0 & ...
            isfinite(bci) & isfinite(dat(:,dind)) & ...
            ((Learly & dat(:,6) <= 0.8) | (~Learly & dat(:,5) < 0.9 & dat(:,6) <= 0.5));
        Lchc = [dat(:,9)==-1 dat(:,9)==1 isfinite(dat(:,9))];

        devdat{mm,1} = nans(num_sessions,3,3);
        devdat{mm,2} = nans(num_sessions,3,3);

        for ss = 1:num_sessions
            Lses = Lgood & dat(:,1) == sessions(ss);
            disp([ss sessions(ss) sum(Lses)]);

            if sum(Lses) > 100
                
                for cc = 1:3
                    Lc = Lchc(:,cc) & Lses;
                    
                    if sum(Lc) > 50
                        
                        % save mean (fit) biases
                        devdat{mm,1}(ss,:,cc) = [ ...
                            nanmean(abs(pmfdat{mm,4}(Lc,3,5))), ...
                            nanmean(dat(Lc,dind)), ...
                            nanmedian(dat(Lc,dind))];
                        
                        % get partial correlations, confidence interval
                        [R1,P1] = partialcorr(xdat(Lc,:), ydat(Lc,:), 'type', 'Spearman');
                        V       = (1+R1(2,1).^2/2)*(1-R1(2,1).^2).^2./(sum(Lc)-3);
                        devdat{mm,2}(ss,:,cc) = [R1(2,1), P1(2,1), V];
                    end
                end
            end
        end
    end
    
    % save it to disk
    FS_saveProjectFile('2008_Bias', 'figBIAS_deviationSummary', devdat);
end

%% PLOTZ

%% Example data
ex_monk  = 'Atticus';
ex_ss    = 60;
ex_chc   = 1;
dat      = FS_getDotsTrainingData(ex_monk);
sessions = nonanunique(dat(:,1));
Lgood    = dat(:,1)==sessions(ex_ss) & dat(:,2) <= 2 & isfinite(dat(:,20));

%% Deviation end points
gr = 0.5*ones(1,3);
axes(axs(4)); cla reset; hold on;
Lrt = Lgood & dat(:,3)==1 & dat(:,9)==1;
Llt = Lgood & dat(:,3)==1 & dat(:,9)==-1;
plot(dat(Lrt,21), dat(Lrt,22), 'k.', 'MarkerSize', 3);
plot(dat(Llt,21), dat(Llt,22), 'kx', 'MarkerSize', 5);
plot(0,0,'o', 'Color', gr, 'MarkerFaceColor', gr);
axis([-20 20 -20 20])
xlabel('Horizontal pos (deg)');
ylabel('Vertical pos (deg)');

%% plot deviations, correct only
% cor coh time choice dev
axes(axs(2)); cla reset; hold on;
Lgood = Lgood & dat(:,9) == ex_chc;
edat  = dat(Lgood,[3 5 6 20]);
cohs  = nonanunique(edat(:,2));
ncohs = length(cohs);
co    = repmat((0.6:-0.1:0)', 1, 3); % {'k' 'y' 'c' 'm' 'r' 'g' 'b'};
tmm   = [150 700];
plot(tmm, [0 0], 'k:');
for cc = 1:ncohs
    dc = edat(edat(:,2)==cohs(cc), :);
    [Y,I] = sort(dc(:,3));
    plot(Y*1000,nanrunmean(dc(I,4), 15), '-', 'Color', co(cc,:));
end
ylabel('Deviation (deg)');
xlabel('Viewing time (msec)');
axis([tmm -0.5 2.0]);

% filtered choices, signed by choice (same/diff)
axes(axs(5)); cla reset; hold on;
mm    = find(strcmp(ex_monk, getBIAS_monks));
Lgood = Lgood & isfinite(fcs{mm});
xs    = fcs{mm}(Lgood).*dat(Lgood,9); % 'biased-choice index': wk-filt x choice
ys    = dat(Lgood,20);
plot(xs, ys, 'k.', 'MarkerSize', 4);
h=lsline;
[R,P] = corr(xs,ys,'type','Spearman');
disp(sprintf('%.3f p=%.3f', R, P));
set(h,'LineStyle', '--','LineWidth', 2);
plot([-5 5], [0 0], 'k:');
plot([0 0], [-5 5], 'k:');
ylabel('Deviation (deg)');
xlabel('Sequential bias');
axis([-.3 .3 -2.3 4.5])

for mm = 1:2
    xax = (1:mse(mm))';
    
    xs = reshape(devdat{mm,1}(xax,1,1:2),[],1);
    ys = reshape(devdat{mm,2}(xax,1,1:2),[],1).*reshape(sign(devdat{mm,1}(xax,2,1:2)),[],1);
    ps = reshape(devdat{mm,2}(xax,2,1:2),[],1);
    es = reshape(sqrt(devdat{mm,2}(xax,3,1:2)),[],1);
    %xs(xs>1.5) = 1.5;
    
    plotBIAS_rho(xs, ys, es, ps, axs(mm*3), 0.7);    
    xlim([0 2]);
    xlabel('Bias |sp/s|')
    ylabel('Correlation')
end

