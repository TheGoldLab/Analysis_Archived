function fig_ = figBIAS_biasMagnitude(num)
% function fig_ = figBIAS_biasMagnitude(num)
%

if nargin < 1 || isempty(num)
    num = 10;
end

[monks,monkn,mse] = getBIAS_monks;

%% Set up Fig
% units should be in inches, from wysifig
wid  = 6.5; % total width
hts  = 1.0;
cw   = 4;
cols = {cw,cw};
[axs,fig_] = getBIAS_axes(num, wid, hts, cols);

%% get the data...
biasdat  = FS_loadProjectFile('2008_Bias', 'figBIAS_biasMagnitude');

if isempty(biasdat)

    biasdat = cell(monkn, 1);
    pmfdat  = FS_loadProjectFile('2008_Bias', 'figBIAS_pmfCompare');

    for mm = 1:monkn

        dat          = FS_getDotsTrainingData(monks{mm});
        Lgood        = dat(:,2) <= 2 & dat(:,3)>=0 & isfinite(pmfdat{mm,4}(:,1,1));
        sessions     = unique(dat(:,1));
        num_sessions = length(sessions);

        ps           = 0.5+(abs(0.5-pmfdat{mm,4}(:,1,1)));
        L0           = dat(:,5)==0;
        Lp           = [L0, ...
            ~L0&ps>0.5&ps<=0.6, ~L0&ps>0.55&ps<=0.65, ...
            ~L0&ps>0.6&ps<=0.7, ~L0&ps>0.65&ps<=0.75, ...
            ~L0&ps>0.7&ps<=0.8, ~L0&ps>0.75&ps<=0.85, ...
            ~L0&ps>0.8&ps<=0.9, ~L0&ps>0.85&ps<=0.95, ...
            ~L0&ps>0.9];
        nlp           = size(Lp,2);
        biasdat{mm,1} = nans(num_sessions, nlp, 2); % 1 t all, corr, err; dc; wk; wk+dc; fcs

        for ss = 1:num_sessions

            Lses  = Lgood & dat(:,1) == sessions(ss);
            disp([ss sessions(ss) sum(Lses)])

            for pp = 1:nlp
                Fp = find(Lses([end; (1:end-1)']) & Lses & Lp(:,pp));

                if length(Fp) > 20
                    
                    biasdat{mm}(ss,pp,1) = 0.5+abs(0.5-sum(dat(Fp,9)==...
                        sign(pmfdat{mm,4}(Fp,3,6)))./length(Fp));

                    biasdat{mm}(ss,pp,2) = 0.5+abs(0.5-sum(dat(Fp,8)==...
                        sign(pmfdat{mm,4}(Fp,3,6)))./length(Fp));
                end
            end
        end
    end

    % save it to disk
    FS_saveProjectFile('2008_Bias', 'figBIAS_biasMagnitude', biasdat);
end

%% PLOTZ
xax = (0.5:0.05:0.95)';
for mm = 1:4
    ses = (1:mse(mm))';

    % vs coh
    axes(axs(mm)); cla reset; hold on;
    ys = biasdat{mm}(ses,:,1);
    ds = biasdat{mm}(ses,:,2);
    plot(xax, nanmedian(ys), 'k.-', 'LineWidth', 2);
    plot(xax, nanmedian(ds), 'k.--', 'LineWidth', 1);

    % indicate sig diff
    for xx=1:10
        Lg = isfinite(ys(:,xx)) & isfinite(ds(:,xx));
        p1 = signtest(ys(Lg,xx), ds(Lg,xx));
        if p1<0.01
            plot(xax(xx), nanmedian(ys(:,xx)), 'k*');
        end
    end
    
    axis([0.5 0.95 0.5 0.65]);

    % vs session
    axes(axs(mm+4)); cla reset; hold on;
    vs = nanmedian(ys(:,1:3)');
    plot(ses, vs, 'k.', 'MarkerSize', 7);
    %    plot(ses, nanmedian(ds(:,1:3)'), 'b.', 'MarkerSize', 1);
    
    Lg = isfinite(vs);
    A = [ones(sum(Lg),1) ses(Lg)];
    [b,bint] = regress(vs(Lg)', A);
    disp([b(2) bint(2,:)])
    
    lsline
    axis([1 ses(end) 0.5 0.8]);
end
