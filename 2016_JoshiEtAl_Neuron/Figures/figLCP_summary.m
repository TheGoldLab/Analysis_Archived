function figLCP_summary(num)

if nargin < 1 || isempty(num)
    num = 7;
end

%% Set up figure
wid     = 11.6; % total width
cols    = {1,1,1};
hts     = repmat(2.5, 1, length(cols));
[axs,~] = getPLOT_axes(num, wid, hts, cols, .8, .5, [], 'Joshi et al', true);
set(axs,'Units','normalized');


PV         = 0.05;
NUM_SITES  = 5;
NUM_MONKS  = 2;
TFIGS      = [3 4 5];
NUM_TTESTS = length(TFIGS);

%% data structure for timing measures
tdat = cell(NUM_SITES, NUM_TTESTS, NUM_MONKS);

%% FIG 3: STA
% load data from file
%  udat{ss,mm,1} is:pupil measure
%   1: per unit (multi first, if appropriate)
%   2: NUM_BINS
%   3: 1-mean, 2-sem, 3-shuffled mean, 4=ps
%  udat{ss,mm,2} is list of file index, unit indices
ptype = 'slope';
udat = FS_loadProjectFile('2013_LCPupil', sprintf('STA_%s', ptype));
xax  = (-1000:2000)';
for ss = 1:NUM_SITES
    for mm = 1:NUM_MONKS
        tdat{ss,1,mm} = nans(size(udat{ss,mm,1},1),3);
        for uu = 1:size(udat{ss,mm,1},1)
            % find longest contiguous + peak
            vals = udat{ss,mm,1}(uu,:,1)-udat{ss,mm,1}(uu,:,3);
            Lp = (vals>0 & udat{ss,mm,1}(uu,:,4)<PV)';
            d = diff([false;Lp;false]);
            p = find(d==1);
            q = min(length(xax),find(d==-1));
            Lgood = xax(q)>-100&xax(p)<700;
            if any(Lgood)
                p = p(Lgood);
                q = q(Lgood);
                Lpq = (q-p)>=75;
                if sum(Lpq)>1
                    ixs = xax(p(Lpq)+round((q(Lpq)-p(Lpq))./2));
                    if any(ixs)>0
                        ix = find(ixs==min(ixs(ixs>0)));
                    else
                        ix = find(ixs==max(ixs));
                    end
                    Fpq = find(Lpq);
                    ix  = Fpq(ix);
                    maxlen = q(ix)-p(ix);
                end
                
                [maxlen,ix] = max(q-p);
                pvals = smooth(vals(p(ix):q(ix)-1),100,'sgolay');
                maxi = p(ix)-1+find(pvals==max(pvals),1);
                
                %                 cla reset; hold on;
                %                 Lp=udat{ss,mm,1}(uu,:,4)<PV;
                %                 plot(xax, vals, 'k-');
                %                 plot(xax(Lp), vals(Lp), 'r.');
                %                 plot(xax(p(ix):q(ix)-1), pvals, 'y-');
                %                 plot(xax(maxi), max(pvals), 'go');
                tdat{ss,1,mm}(uu,:) = [maxlen, xax(maxi), max(pvals)];
            end
            % r = input('next')
        end
    end
end

%% FIG 4: PETH
% load data from file
% udat =
% 1:PETH, 2:corr, 3:[monk file unit], 4:pd stats, 5:example raster
% data, 6: example raster selection array
% load data from file
udat = FS_loadProjectFile('2013_LCPupil', 'PETH');
xax  = (-1875:10:875)';
for ss = 1:NUM_SITES
    for mm = 1:NUM_MONKS
        Fm    = find(udat{ss,3}(:,1) == mm);
        tdat{ss,2,mm} = nans(length(Fm),3);
        for uu = 1:length(Fm)            
            vals = udat{ss,1}(Fm(uu),:,1)-udat{ss,1}(Fm(uu),:,2);
            Lp = (vals>0 & udat{ss,2}(Fm(uu),:,1)<PV)';
            d = diff([false;Lp;false]);
            p = find(d==1);
            q = min(find(d==-1), length(xax));
            Lgood = xax(q)>-800&xax(p)<100;

            %cla reset; hold on;
            %Lp=udat{ss,2}(Fm(uu),:,1)<PV;
            %plot(xax, vals, 'k-');
            %plot(xax(Lp), vals(Lp), 'r.');
                
            if any(Lgood)
                p = p(Lgood);
                q = q(Lgood);
                [maxlen,ix] = max(q-p);
                pvals = smooth(vals(p(ix):q(ix)-1),100,'sgolay');
                maxi = p(ix)-1+find(pvals==max(pvals),1);
                tdat{ss,2,mm}(uu,:) = [maxlen, xax(maxi), max(pvals)];

                %plot(xax(p(ix):q(ix)-1), pvals, 'y-');
                %plot(xax(maxi), max(pvals), 'go');
            end
            %r = input('next')
        end
    end
end


%% collect and plot peak times
sites        = {...
    'LC'    {'Oz' 'Cicero'} ; ...
    'IC'    {'Oz' 'Cicero'} ; ...
    'SC'    {'Oz' 'Cicero'}     ; ...
    'ACC'   {'Sprout' 'Atticus'}; ...
    'PCC'   {'Sprout' 'Cheetah'}; ...
    'MT'    {'Cheetah'}         ; ...
    };
cvs = [50 5];
yl  = [-300 900; -800 400];
for ff = 1:2
    axes(axs(ff)); cla reset; hold on;
    for ss = 1:NUM_SITES
        for mm = 1:NUM_MONKS
            [co, sy] = getLCP_colorsAndSymbol(sites{ss,2}{mm}, sites{ss,1});
            
            Lg = tdat{ss,ff,mm}(:,1) > cvs(ff) & tdat{ss,ff,mm}(:,2)>yl(ff,1) & ...
                tdat{ss,ff,mm}(:,2)<yl(ff,2);
            vs = tdat{ss,ff,mm}(Lg,2);
            xb = ss+0.3*(mm-1);
            xs = xb+rand(sum(Lg),1)*0.2-0.1;
            plot(xs,vs,sy,'Color',co{1},'MarkerFaceColor',co{1}, 'MarkerSize', 4);
            plot(xb+[-0.1 0.1], median(vs).*[1 1], 'k-', 'LineWidth', 2);
            axis([0.8 5.5 yl(ff,:)])
            text(xb-0.2, yl(ff,2)-50, ...
                sprintf('%d/%d', sum(Lg), sum(isfinite(tdat{ss,ff,mm}(:,2)))))
        end
    end
end


%
% SFIGS      = [2 3 4 5 6];
% NUM_STESTS = length(SFIGS);
% %% data structure for fraction of sites with sig effects
% sdat = nans(NUM_SITES, 4, NUM_STESTS, NUM_MONKS); % 2nd is +count, -count, total, pop sig
% 
% %% FIG 2: TONIC
% % load data from file
% udat = FS_loadProjectFile('2013_LCPupil', 'Tonic');
% for ss= 1:NUM_SITES
%     for mm = 1:NUM_MONKS
%         Lm = udat{ss,1}(:,3) == mm;
%         ys = udat{ss,1}(Lm,1);
%         Lp = udat{ss,1}(Lm,2)<PV;
%         sdat(ss,:,1,mm) = [ ...
%             sum(Lp&ys>0), ...
%             sum(Lp&ys<0), ...
%             sum(Lm), ...
%             signrank(ys)];
%     end
% end
