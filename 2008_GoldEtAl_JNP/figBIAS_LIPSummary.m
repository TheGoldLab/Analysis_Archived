function fig_ = figBIAS_LIPSummary(num, ex_monk, ex_ss)
%function fig_ = figBIAS_LIPSummary(num, ex_monk, ex_ss)
%

if nargin < 1 || isempty(num)
    num = 13;
end
if nargin < 2 || isempty(ex_monk)
    ex_monk = 'Cyrus';
end
if nargin < 3 || isempty(ex_ss)
    ex_ss = 93;
end

monkc = {'Cyrus' 'Cy' 3; 'ZsaZsa' 'ZZ' 4};
monkn = 2;

%% Set up Fig
% units should be in inches, from wysifig
wid  = 6; % total width
hts  = 0.8;
cols = {1,6,6,6};
[axs,fig_] = getBIAS_axes(num, wid, hts, cols);

% get the data
lipdat = FS_loadProjectFile('2008_Bias', 'figBIAS_LIPSummary');

if isempty(lipdat)

    lipdat = cell(size(monkc,1), 3);

    % R = [cell #, session #, coh, pref_dir, correct, slope];
    slopes = FS_loadProjectFile('2008_Bias', 'figBIAS_LIPSlope');

    [a,b] = unix('hostname');
    %ddir = ['/Volumes/XServerData/Physiology/MTLIP/' monk{1} '/TRain/'];
    if strncmp(b,'labs-computer', 10)
        ddir = '/Volumes/George/Data/MTLIP';
    else
        ddir = '/Users/jigold/Desktop/data/MTLIP/';
    end

    for mm = 1:monkn

        a        = getML_txt([ddir monkc{mm,2} 'TRain_LIP.txt']);
        fn       = a.data{strcmp(a.name,'dat_fn')};
        usable   = find(a.data{strcmp(a.name,'usable')}==1);
        uid      = a.data{strcmp(a.name,'uid')};

        lipdat{mm,1} = nans(length(usable), 2, 7);
        lipdat{mm,2} = nans(length(usable), 3, 7, 6); % b/sem/n; 7 selectors; trg/pre/peri/post/sac
        lipdat{mm,3} = {};

        for ss = 1:length(usable) %ss=100 %length(ses):-1:1

            disp(fn{usable(ss)})
            openFIRA([ddir monkc{mm,1} '/' fn{usable(ss)}]);

            % get spike id, slopes
            sid   = getFIRA_spikeByID(uid(usable(ss)));
            sldat = slopes{mm}(slopes{mm}(:,1)==usable(ss),6);

            [lipdat{mm,1}(ss,:,:),lipdat{mm,2}(ss,:,:,:),exd,wk] = ...
                getBIAS_LIPdata(sid, sldat, ex_ss == ss && strcmp(ex_monk, monkc{mm,1}));

            % save wk in slopes... trust me
            slopes{mm}(slopes{mm}(:,1)==usable(ss),7) = wk;

            if ~isempty(exd)
                lipdat{mm,3} = exd;
            end
        end
    end

    FS_saveProjectFile('2008_Bias', 'figBIAS_LIPSummary', lipdat);
end

%% Plotz
%
% get xample data
exdat = lipdat{1,3};

% PSTH
axes(axs(1)); cla reset; hold on;
rmax  = 75;
ucohs = [0 3.2 6.4 2.8 25.6 51.2 99.9];
ps    = [100 390; 410 700; 950 1250; 2800 3100; 3150 3300];
hs    = nans(size(ps,1)+1,1);
for pp = 1:size(ps,1)
    hs(pp)=patch([ps(pp,1) ps(pp,2) ps(pp,2) ps(pp,1)], [1 1 100 100], [0.9 0.9 0.9]);
end
hs(end) = patch([1400 1700 1700 1400], [1 1 100 100], [0.8 0.8 0.8]);
set(hs,'LineStyle', 'none');
lst = {'-' '--'};
co  = repmat(((length(ucohs)-1:-1:0)./length(ucohs))',1,3);
bs = exdat{2};
es = exdat{3};
lo = [0 850 2400];
for bb = 1:length(bs)
    for pp = 1:2
        for cc = 1:length(ucohs)
            plot(lo(bb)+mean(bs{bb},2), es{bb}(:,cc,pp), lst{pp}, 'Color', co(cc,:));
        end
    end
end
set(gca, ...
    'XTick',      [400      1250     1750   2500     3100], ...
    'XTickLabel', {'trg on' 'dot on' '500' 'dot off' 'fp off'});
axis([0 3400 0 rmax]);
xlabel('Time (ms)')
ylabel('Response (sp/s)');

%% Examples
regs = exdat{1};
yl   = [70 70 70 240 110 120];
for xx = 1:length(regs)
    axes(axs(1+xx)); cla reset; hold on;
    plot([-0.3 0.3], [0 0], 'k:');
    plot([0 0], [-5 5], 'k:');
    plot(regs{xx}(:,1), regs{xx}(:,2), 'k.');
    [R,RLO,RUP,P] = corrcoef(regs{xx}(:,1), regs{xx}(:,2));

    disp(sprintf('%.2f,p=%.2f', R(2,1),P(2,1)))
    axis([-0.3 0.3 0 yl(xx)]);
end

%% summaries
% 1=all, 2=tin, 3=tout, 4=tin/cor, 5=tin/err, 6=tout/cor, 7=tout/err
inds = [4 6];
mx  = 1.8;
for mm = 1:2
    xs = reshape(lipdat{mm,1}(:,1,inds),[],1);
    xs(xs>mx)=mx;
    for pp = 1:6
        ys = reshape(lipdat{mm,2}(:,1,inds,pp),[],1);
        ps = reshape(lipdat{mm,2}(:,2,inds,pp),[],1);
        es = reshape(sqrt(lipdat{mm,2}(:,3,inds,pp)),[],1);
        plotBIAS_rho(xs, ys, es, ps, axs((mm-1)*6+pp+7), 0.7);
        xlim([0 mx]);
    end
end
