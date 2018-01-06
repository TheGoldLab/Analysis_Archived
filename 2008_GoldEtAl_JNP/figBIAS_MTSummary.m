function fig_ = figBIAS_MTSummary(num, ex_monk, ex_ss)
% function fig_ = figBIAS_MTSummary(num, ex_monk, ex_ss)
%

if nargin < 1 || isempty(num)
    num = 12;
end
if nargin < 2 || isempty(ex_monk)
    ex_monk = 'Cyrus';
end
if nargin < 3 || isempty(ex_ss)
    ex_ss = 82;
end

% get monks
monkc = {'Cyrus' 'Cy' 3; 'ZsaZsa' 'ZZ' 4};
monkn = 2;

%% Set up Fig
% units should be in inches, from wysifig
wid  = 3;   % total width
hts  = 0.8;
cols = {1,3,3,3};
[axs,fig_] = getBIAS_axes(num, wid, hts, cols);

%% get the data
mtdat = FS_loadProjectFile('2008_Bias', 'figBIAS_MTSummary');

if isempty(mtdat)

    mtdat  = cell(monkn,3);
    ddir   = '/Users/jigold/Desktop/data/MTLIP/';
    % nrand  = 1000;
    % rs     = nans(nrand,3);
    
    for mm = 1:monkn

        a          = getML_txt([ddir monkc{mm,2} 'TRain_MT.txt']);
        fn         = a.data{strcmp(a.name,'dat_fn')};
        uid        = a.data{strcmp(a.name,'uid')};
        usable     = find(a.data{strcmp(a.name,'usable')}+a.data{strcmp(a.name,'train')}==2);
        num_usable = length(usable);
        sessions   = a.data{strcmp(a.name,'session')};

        mtdat{mm,1} = nans(num_usable, 2, 8);
        mtdat{mm,2} = nans(num_usable, 3, 8, 3); % b/sem/n; 8 selectors; pre/peri/post
        mtdat{mm,3} = {};

        for ss = 1:num_usable %ss=45 %length(ses):-1:1

            disp(fn{usable(ss)})
            openFIRA([ddir monkc{mm,1} '/' fn{usable(ss)}]);
            sid = getFIRA_spikeByID(uid(usable(ss)));

            % call getBIAS_MTdata to do the work
            [mtdat{mm,1}(ss,:,:), mtdat{mm,2}(ss,:,:,:), exd] = ...
                getBIAS_MTdata2(sid, ex_ss == sessions(usable(ss)) && strcmp(ex_monk, monkc{mm,1}));
            
            % get example data
            if ~isempty(exd)
                mtdat{mm,3} = exd;
            end 
        end
    end

    % save it to disk
    FS_saveProjectFile('2008_Bias', 'figBIAS_MTSummary', mtdat);
end

exdat = mtdat{strcmp(ex_monk, monkc(:,1)), 3};

%% PSTH
axes(axs(1)); cla reset; hold on;
ucohs = [0 3.2 6.4 12.8 25.6 51.2 99.9];
pc    = [0.9 0.9 0.9];
h(1)=patch([100 400 400 100],     [1 1 100 100], pc);
h(2)=patch([500 1475 1475 500],   [1 1 100 100], pc);
h(3)=patch([1900 2200 2200 1900], [1 1 100 100], pc);
set(h,'LineStyle', 'none');
lst = {'-' '--'};
co  = repmat(((length(ucohs)-1:-1:0)./length(ucohs))',1,3);
for pp = 1:2
    for cc = 1:length(ucohs)
        plot(mean(exdat{4},2), exdat{5}(:,cc,pp), lst{pp}, 'Color', co(cc,:));
        plot(mean(exdat{6},2)+2200, exdat{7}(:,cc,pp), lst{pp}, 'Color', co(cc,:));
    end
end
set(gca,'XTick', [400 900 1400 1650 2200], 'XTickLabel', {'dot on' 500 '' 'dot off' 'fp off'});
axis([0 2400 0 100]);
xlabel('Time (ms)')
ylabel('Response (sp/s)');

%% Examples
for xx = 1:3
    axes(axs(1+xx)); cla reset; hold on;
    plot(exdat{xx}(:,1), exdat{xx}(:,2), 'k.', 'MarkerSize', 4);
    h = lsline;
    set(h, 'LineWidth', 2, 'LineStyle', '--');
    plot([-1 1], [0 0], 'k:');
    plot([0 0], [-10 10], 'k:');
    if xx == 1 || xx == 3
        axis([-.2 .2 0 50]);
    else
        axis([-.2 .2 0 100]);
    end
end

%% Plot summaries -- histograms
%
% ind: 1-all, 2-pref, 3-null, 4-pref/cor, 5-pref/err, 6-null/cor, 7-null/err, 8-0% coh
inds = [4 6];
mx  = 2.5;
for mm = 1:2
    xs = reshape(mtdat{mm,1}(:,1,inds),[],1);
    xs(xs>mx)=mx;
    for pp = 1:3
        ys = reshape(mtdat{mm,2}(:,1,inds,pp),[],1);
        ps = reshape(mtdat{mm,2}(:,2,inds,pp),[],1);
        es = reshape(sqrt(mtdat{mm,2}(:,3,inds,pp)),[],1);
        plotBIAS_rho(xs, ys, es, ps, axs((mm-1)*3+pp+4), 0.7);
        xlim([0 mx]);
    end
end
