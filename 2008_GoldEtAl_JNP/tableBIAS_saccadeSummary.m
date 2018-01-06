function tableBIAS_saccadeSummary
% function tableBIAS_saccadeSummary
%

% if nargin < 1 || isempty(num)
%     num = 1;
% end

%% Set up Fig
% units should be in inches, from wysifig
% wid  = 5.0; % total width
% hts  = 1;
% cols = {3,3,3,3};
% [axs,fig_] = getBIAS_axes(num, wid, hts, cols);

% get monk names
[monks,monkn] = getBIAS_monks;

% get the data
pmfdat = FS_loadProjectFile('2008_Bias', 'figBIAS_pmfCompare');
sacdat = FS_loadProjectFile('2008_Bias', 'tableBIAS_saccadeSummary');

if isempty(sacdat)

    sacdat = cell(monkn, 2);
    fcs    = getBIAS_fcs;
    inds   = [10 14 15 16]; % latency; vavg; vmax; accuracy
    ni     = length(inds);

    for mm = 1:monkn

        dat          = FS_getDotsTrainingData(monks{mm});
        sessions     = unique(dat(:,1));
        num_sessions = length(sessions);

        bci          = fcs{mm}.*dat(:,9); % 'biased-choice index': wk-filt x choice
        ts           = (max(0.1, dat(:,6))-0.1).^0.5;
        ydat         = [dat(:,5).*ts, ts];

        Learly       = dat(:,1) < sessions(round(num_sessions/4));
        Lgood        = dat(:,2) < 3 & dat(:,3) >= 0 & isfinite(bci) & ...
            ((Learly & dat(:,6) <= 0.8) | (~Learly & dat(:,5) < 0.9 & dat(:,6) <= 0.5));
        if mm <= 2
            Lgood = Lgood & dat(:,17) == 0;
        end
        Lchc         = [dat(:,9)==-1 dat(:,9)==1 isfinite(dat(:,9))];
        Lsacs        = isfinite(dat(:,inds));

        sacdat{mm,1} = nans(num_sessions,3,3,ni);
        sacdat{mm,2} = nans(num_sessions,3,3,ni);

        for ss = 1:num_sessions
            Lses = Lgood & dat(:,1) == sessions(ss);
            disp([ss sessions(ss) sum(Lses)]);

            if sum(Lses) > 100

                for nn = 1:ni
                    Ln = Lses & Lsacs(:,nn);

                    for cc = 1:3
                        Lc = Ln & Lchc(:,cc);

                        if sum(Lc) > 50

                            % save mean (fit) biases
                            sacdat{mm,1}(ss,:,cc,nn) = [ ...
                                nanmean(abs(pmfdat{mm,4}(Lc,3,5))), ...
                                nanmean(dat(Lc,inds(nn))), ...
                                nanmedian(dat(Lc,inds(nn)))];

                            % get partial correlations, confidence interval
                            [R1,P1] = partialcorr([bci(Lc) dat(Lc,inds(nn))], ydat(Lc,:), 'type', 'Spearman');
                            V       = (1+R1(2,1).^2/2)*(1-R1(2,1).^2).^2./(sum(Lc)-3);
                            sacdat{mm,2}(ss,:,cc,nn) = [R1(2,1), P1(2,1), V];
                        end
                    end
                end
            end
        end
    end

    % save it to disk
    FS_saveProjectFile('2008_Bias', 'tableBIAS_saccadeSummary', sacdat);
end

%% FOR A TABLE
inds = 1:2;
for mm = 1:monkn
    for dd = 1:4
        ys = reshape(sacdat{mm,2}(:,1,inds,dd),[],1);
        ps = reshape(sacdat{mm,2}(:,2,inds,dd),[],1);
        es = reshape(sacdat{mm,2}(:,3,inds,dd),[],1);
        Lf = isfinite(ys);
        p  = signtest(ys(Lf));
        np = sum(ys>0&ps<0.05);
        nm = sum(ys<0&ps<0.05);
        disp(sprintf('mm=%d,pp=%d: %.3f [%.3f], p=%.3f, n=%d, np=%d, nm=%d', ...
            mm, dd, median(ys(Lf)), iqr(ys(Lf)), p, sum(Lf), np, nm))
    end
end

return


% bci  = fcs{mm}.*dat(:,9); % 'biased-choice index': wk-filtxchoice
cohs = nonanunique(dat(:,5));
ncoh = length(cohs);
Lcoh = false(size(dat,1), ncoh);
for cc = 1:ncoh
    Lcoh(:,cc) = dat(:,5) == cohs(cc);
end
Lgood = isfinite(bci) & dat(:,2) < 3 & dat(:,3) >= 0;
if mm <= 2
    Lgood = Lgood & dat(:,17) == 0;
end
Lsacs = isfinite(dat(:,inds));
for ii = 1:length(inds)
    dat(~Lsacs(:,ii), inds(ii)) = nan;
end

sacdat{mm,1} = nans(num_sessions,3,3,3,4); % r/cis/p; cor/err/both; lat/vavg/vmax/accuracy
sacdat{mm,2} = nans(size(dat,1),4);      % de-trended data

for ss = 1:num_sessions

    Lses = Lgood & dat(:,1) == sessions(ss);
    disp([ss sessions(ss) sum(Lses)]);

    % de-trend the saccade parameters by choice/correct/coh/time
    for ii = 1:3 % choice
        for jj = 1:3 % correct/error
            for dd = 1:4
                Ld = Lses & Lchc(:,ii) & Lcor(:,jj) & Lsacs(:,dd);

                if ii < 3 && jj < 3
                    for cc = 1:ncoh % coherence
                        Fdt = find(Ld & Lcoh(:,cc));
                        if length(Fdt) > 10
                            [Y,I] = sort(dat(Fdt,6));
                            dfi   = dat(Fdt(I),inds(dd));
                            sacdat{mm,2}(Fdt(I),dd) = ...
                                zscore(dfi-nanrunmean(dfi,ceil(length(Fdt)/4)));
                        else
                            sacdat{mm,2}(Fdt,dd) = nanzscore(dat(Fdt,inds(dd)));
                        end
                    end
                end

                if sum(Ld) > 5
                    [b,sem] = lscov([ones(sum(Ld),1) bci(Ld)],  sacdat{mm,2}(Ld,dd));
                    sacdat{mm,1}(ss,:,ii,jj,dd) = [b(2) sem(2) sum(Ld)];
                end
            end
        end
    end
end

% save it to disk
FS_saveProjectFile('2008_Bias', 'figBIAS_saccadeSummary', sacdat);


%% FOR A TABLE
inds = 1:2;
for mm = 1:monkn
    for dd = 1:4
        ys = reshape(sacdat{mm,1}(:,1,inds,3,dd),[],1);
        es = reshape(sacdat{mm,1}(:,2,inds,3,dd),[],1);
        ns = reshape(sacdat{mm,1}(:,3,inds,3,dd),[],1);
        Lf = isfinite(ys) & isfinite(es) & ns > 10;
        p  = signtest(ys(Lf));
        np = sum((ys - es*1.96) > 0);
        nm = sum((ys + es*1.96) < 0);
        disp(sprintf('mm=%d,pp=%d: %.3f [%.3f], p=%.3f, n=%d, np=%d, nm=%d', ...
            mm, dd, median(ys(Lf)), iqr(ys(Lf)), p, sum(Lf), np, nm))
    end
end

return


inds = 1:2;
for mm = 1:monkn
    for dd = 1:4
        subplot(4,4,dd+(mm-1)*4); cla reset; hold on;
        ys = reshape(sacdat{mm,1}(:,1,inds,3,dd),[],1);
        hist(ys, -1:.1:1);
        title(sprintf('%.3f,p=%.3f', nanmedian(ys), signtest(ys(isfinite(ys)))))
        xlim([-1 1]);
    end
end



%% REPORT CORRELAIONS between GLOBAL BAIS, latency differences
% for mm =1:monkn
%
%     xs = sesdat{mm,2}(:,5,end);
%     ys = figdat{mm,2}(:,1,1)-figdat{mm,2}(:,3,1);
%     Lgood = isfinite(xs) & isfinite(ys);
%     [R,P,RLO,RUP] = corrcoef(xs(Lgood), ys(Lgood));
%     disp(sprintf('%s: %.2f [%.2f %.2f], p=%.3f', monks{mm}, R(2,1), RLO(2,1), RUP(2,1), P(2,1)))
% end

%% plot data from both monks, separate by correct/error
axlim  = [-.1 .6 -.55 .55];
ps     = loadBIAS_figData('getBIAS_pData');
dt     = [1 2 4]; % data types (latency, vavg, accuracy)
for mm = 1:monkn
    for dd = 1:3 % data types
        plotBIAS_regress(axs((mm-1)*3+dd), sesdat{mm,1}, ...
            [figdat{mm}(:,1:3,3,dt(dd)) ps{mm}(:,2+dt(dd))], axlim);
        %            figdat{mm}(:,1:3,3,dt(dd)), axlim);
    end
end

axes(axs(1));
title('latency');
axes(axs(2));
title('vavg');
axes(axs(3));
title('accuracy');
