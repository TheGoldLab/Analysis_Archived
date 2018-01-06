function fig_ = figBIAS_wkFits(num, ex_monk, ex_ss)
% function fig_ = figBIAS_wkFits(num, ex_monk, ex_ss)
%

if nargin < 1 || isempty(num)
    num = 5;
end
if nargin < 2 || isempty(ex_monk)
    ex_monk = 'Atticus';
end

if nargin < 2 || isempty(ex_ss)
    ex_ss = 14;
end

% units should be in inches, from wysifig
wid  = 6.0; % total width
ht   = 1;
cols = {2, [0.6 0.4]};
[axs,fig_] = getBIAS_axes(num, wid, ht, cols);

% get monkeys, colors
[monks, monkn, ~] = getBIAS_monks;

% get the data
wkdat = FS_loadProjectFile('2008_Bias', 'figBIAS_wkFits');

if isempty(wkdat)

    wkdat  = cell(length(monks), 5);

    fun    = @ddExp3;
    WKSZ   = 600;
    WKSTEP = 10;
    nvins  = 7; % 6 below plus sum
    nfits  = 2; % raw, fit

    for mm = 1:monkn

        dat          = FS_getDotsTrainingData(monks{mm});
        Lgood        = dat(:,2) <= 2 & sum(isfinite(dat(:,[3 5 6 8 9])),2)==5;
        sessions     = unique(dat(:,1));
        num_sessions = length(sessions);

        % wk x vins x fits x session
        wkdat{mm,1} = nans(WKSZ, nvins, nfits, num_sessions); % wk kernels
        wkdat{mm,2} = nans(WKSZ, nvins, nfits, num_sessions); % wk correlation coefficients
        wkdat{mm,3} = nans(size(dat,1), 2);                   % wk-filtered data, scaled wk-filtered data, resids from simple model, 
        wkdat{mm,4} = nans(num_sessions, 4);
        wkdat{mm,5} = nans(num_sessions, 4, 3); % 4 param kExp2, 3 vins

        % get chc array for vins
        % 1 both
        % 2 cor
        % 3 err
        % 4 cor x rew
        % 5 cor x diff
        % 6 err x diff
        %   (note that the 7th is sum [2 3])
%        p     = ddExp3([30 0 0], dat(Lgood,[5 6 8 9]));
        p     = ddExp3([30 0.01 0], dat(Lgood,[5 6]));
        diffs = zeros(size(dat,1),1);
        diffs(Lgood) = p;
        vins = dat(:,9*ones(1,6));
        vins(dat(:,3)==0,[2 4 5]) = 0;
        vins(dat(:,3)==1,[3 6])   = 0;
        vins(:,4) = vins(:,4).*dat(:,4);
        vins(:,5) = vins(:,5).*diffs;
        vins(:,6) = vins(:,6).*diffs;

        % dfp = nans(num_sessions, 2);
        
        % loop through the sessions
        for ss = 1:num_sessions

            Fses = find(Lgood & dat(:,1) == sessions(ss));
            disp([ss sessions(ss) length(Fses)])

            if length(Fses) > 200

                % fit to simple model
                [f,s,t,p,r] = ctPsych_fit(fun, dat(Fses,[5 6 8 9]), dat(Fses, 3), [], []);
                % [f2,s2,t2,p2,r2] = ctPsych_fit(@ddExp4z, dat(Lses,[5 6 8 9]), dat(Lses, 3), [], []);
                % ctPsych_plot(fun, f, dat(Lses,[5 6 8 9]), dat(Lses, 3))

                % save resids
                wkdat{mm,3}(Fses,2) = r(:,4);
                
                % get wk data using different vins, try all fits
                % [a,b,c,d] = ...
                [wkdat{mm,1}(:,1:6,:,ss), wkdat{mm,2}(:,:,:,ss), c, d] = ...
                    getBIAS_wk(vins(Fses(1:end-1),:), r(2:end,4), WKSZ, WKSTEP, [2 3], {'kExp2'});

                % save first three fits (both/cor/err)
                wkdat{mm,5}(ss,:,:) = permute(d{1}(:,1:3), [3 1 2]);
                
                % recompute with correct len (max r) &
                % save filtered choices for both/fit & r (cis)
                lag_m2 = find(wkdat{mm,2}(:,end,end,ss)==max(wkdat{mm,2}(:,end,end,ss)),1)+WKSTEP;
                [wk, rr, ff, d] = getBIAS_wk(vins(Fses(1:end-1),2:3), ...
                    r(2:end,4), lag_m2, 0, [1 2], {'kExp2'});
                wkdat{mm,3}(Fses(1),1)     = 0;
                wkdat{mm,3}(Fses(2:end),1) = ff(:,3,2);
                wkdat{mm,4}(ss,:)          = [rr(:,3,2)' lag_m2];
            end
        end
    end

    % save it to disk
    FS_saveProjectFile('2008_Bias', 'figBIAS_wkFits', wkdat);
end

%                 [a, b, c, d] = ...
%                     getBIAS_wk(vins(Fses(1:end-1),2), r(2:end,4), WKSZ, WKSTEP, [], {'kExp2', 'kExp1', 'kPow'});
% 
%                 % stats on fits
%                 N      = find(b(:,1,2) == max(b(:,1,2)), 1);
%                 Kexp2  = 5; % kExp2 = 4 params + 1
%                 Kpow   = 3; % kPow  = 2 params + 1
%                 ris    = 1:N;
%                 SSexp2 = sum((a(ris,1,2) - a(ris,1,1)).^2);
%                 SSexp1 = sum((a(ris,1,3) - a(ris,1,1)).^2);
%                 SSpow  = sum((a(ris,1,4) - a(ris,1,1)).^2);
% 
%                 % compute corrected AIC for non-nested models
%                 %   kExp2 and kPow (from GraphPad book): 
%                 %   AIC            = N * log(SS/N)+2K+2K(K-1)/(N-K-1)
%                 %   P(AIC)         = 1/(exp(0.5*AIC)+1)
%                 %   Evidence ratio = p(AIC1)/p(AIC2) 
%                 %                  = exp(-0.5(AIC1-AIC2))
%                 Aexp2  = N*log(SSexp2/N)+2*Kexp2+2*Kexp2*(Kexp2-1)./(N-Kexp2-1);
%                 Apow   = N*log(SSpow /N)+2*Kpow +2*Kpow *(Kpow -1)./(N-Kpow -1);
%                 er     = exp(-0.5*(Aexp2-Apow));
%                 
%                 % compute F-test for nested models (kExp1 and kExp2)
%                 DFexp1 = N - 2; % degrees of freedom = # data points - # parameters
%                 DFexp2 = N - 4;
%                 F      = ((SSexp1-SSexp2)/(DFexp1-DFexp2))/(SSexp2/DFexp2);
%                 
%                 dfp(ss,:) = [fpdf(F, DFexp1, DFexp2), er];
%                 
%                 disp(dfp(ss,:))
%             end
%         end
%         
%         size(find(dfp(:,1)<0.05))
%         size(find(isfinite(dfp(:,1))))
%         size(find(dfp(:,2)>20))
%         size(find(isfinite(dfp(:,2))))
       
                % cla reset; hold on;
                % ind = 2;
                % errorbar(a(1:50,ind,1), e(1:50,ind), 'k-')
                % plot(a(1:50,ind,2), 'g-')

%% Report distributions of correlation coefficients per monkey
for mm = 1:monkn
    ns = size(wkdat{mm,2},4);
    mr = nans(ns,1);
    for ss = 1:ns
        mr(ss) = max(wkdat{mm,2}(:,end,2,ss));
    end
    disp(sprintf('%s: %.2f [%.2f %.2f]', monks{mm}, prctile(mr, [50 2.5 97.5])))
end

%%%%
%% ROW 1 -- examples
%
% Panels A&B
%%%
% wk, raw, correct (A), error (B)
tmax  = 30;
for pp = 1:2
    axes(axs(pp)); cla reset; hold on;
    plot([0 tmax+2], [0 0], 'k:');
    plot(wkdat{strmatch(ex_monk, monks),1}(:,1+pp,1,ex_ss), 'k.', 'MarkerSize', 13, 'LineWidth', 2);
    plot(wkdat{strmatch(ex_monk, monks),1}(:,1+pp,2,ex_ss), 'k--', 'LineWidth', 2);
    axis([0 tmax+2 -0.1 0.1]);
    xlabel('Trial lag')
    if pp == 1
        ylabel('Filter weight')
    end
end

%%%%
%% ROW 2 -- R vs WK Length
%%%%

%% Panel C: r vs lag
%
% summary of rs vs wk length
axes(axs(3)); cla reset; hold on;
mst = {'-' '-' '-', '-'};
co  = {[.5 .5 .5], 'k'};
for mm = 1:monkn
    Lwk = isfinite(wkdat{mm,2}(:,2,1,1));
    xax = (1:size(wkdat{mm,2},1))';
    xax = xax(Lwk);
    % rs x vins x fits x session
    rdat  = permute(wkdat{mm,2}(Lwk,2,[1 2],:), [4 1 3 2]);
    for rr = 1:2
        nm = nanmean(rdat(:,:,rr));
        ns = nanse(rdat(:,:,rr));
        plot([xax xax]', [nm-ns; nm+ns], '-', 'Color', co{rr}, 'LineWidth', mm/2);
        plot(xax, nm, '-', 'Color', co{rr}, 'LineWidth', mm/2);% 'MarkerFillColor);e
    end
end
axis([0 400 0 0.8]);
xlabel('Trial lag')
ylabel('Correlation coefficient')

%% Panel D: Boxplot of max lags
%
% boxplots
lags = [];
for mm = 1:monkn
    lags = cat(1, lags, [mm*ones(size(wkdat{mm,4},1),1) wkdat{mm,4}(:,4)]);
    disp(sprintf('%s: %.2f +- %.2f', monks{mm}, ...
        nanmedian(lags(lags(:,1)==mm,2)), iqr(lags(lags(:,1)==mm,2))))
end
lags = lags(isfinite(lags(:,2)), :);
axes(axs(4)); cla reset; hold on;
h=boxplot(lags(:,2),lags(:,1), ...
    'notch',    'on', ...
    'widths',   0.5, ...
    'colors',   'k', ...
    'symbol',   'k+', ...
    'labels',   {'1' '2' '3' '4'});
for mm = 1:4
    set(h(isfinite(h(:,mm)),mm), 'Color', 'k');
end
set(h(1:2,:),'LineStyle','-');
ylim([-10 600]);

return

%% removed different fits, below

%%%%
%% ROW 3: different fits
%%%%

%%%
%% Panels E-H (2nd row): r comparisons of different vins
%%%
% compare different vins
mst = {'.' 'x' '*', 'o'};
vis = [1 7; 2 4; 2 5; 3 6];
xlabs = {'All', 'Correct', 'Correct', 'Error'};
ylabs = {'Correct + Error', 'Correct x Rew', 'Correct x Diff', 'Error x Diff'};
fi   = 2;
for ii = 1:size(vis,1)
    axes(axs(ii+4)); cla reset; hold on;
    plot([-.6 .6], [-.6 0.6], 'k:');
    plot([-.6 .6], [0 0], 'k-');
    plot([0 0], [-.6 .6], 'k-');
    tts = nans(1,4);
    for mm = 1:4
        
        % find max per session, per vin
        % rs x vins x fits x session -> sess x vins
        rdat = nans(size(wkdat{mm,1},4), 2);
        for ss = 1:size(wkdat{mm,1},4)
            rdat(ss, :) = [ ...
                max(wkdat{mm,2}(:, vis(ii,1), fi, ss)), ...
                max(wkdat{mm,2}(:, vis(ii,2), fi, ss))];

        end
        rdat    = rdat(isfinite(rdat(:,1)), :);
%        [h,p]   = ttest2(rdat(:,1), rdat(:,2));
         p      = ranksum(rdat(:,1), rdat(:,2));
        tts(mm) = p;
        plot(rdat(:,1), rdat(:,2), mst{mm}, ...
            'MarkerSize', 2, 'Color', 'k');
    end
    disp(tts)
    axis([-.2 .6 -.2 .6]);
    xlabel(xlabs{ii});
    ylabel(ylabs{ii});
end

% compare fits

