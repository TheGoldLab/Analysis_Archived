function fig_ = figBIAS_pmfParams(num)
% function fig_ = figBIAS_pmfParams(num)
%

if nargin < 1 || isempty(num)
   num = 9;
end

[monks,monkn,~] = getBIAS_monks;

%% Set up Fig
% units should be in inches, from wysifig
wid  = 6.5; % total width
hts  = 1.0;
cw   = 4;
cols = {cw,cw,cw};
[axs,fig_] = getBIAS_axes(num, wid, hts, cols);
set(axs,'Units', 'Normalized');

%% get the data...
pmfdat = FS_loadProjectFile('2008_Bias', 'figBIAS_pmfCompare');

%% Plotit
%inds = [1 1; 6 2; 8 2];
ylm  = [75 1.5 1.2; 50 1 1.2; 90 1.5 1.2; 30 2 1.2];
fi   = 6;
for mm = 1:monkn
   
   % ses = (1:mse(mm))';
   ses = (1:size(pmfdat{mm,1},1))';
   sax = [ses ses]';
   A   = [ses ones(size(ses))];
   Llapse = pmfdat{mm,2}(ses,6,1,1) < 0.1;
   
   for ii = 1:3
      axes(axs(mm+(ii-1)*4)); cla reset; hold on;
      if ii == 1
         ys = min(ylm(mm,ii), pmfdat{mm,2}(ses,1,1,fi));
         es = max(0.01,       pmfdat{mm,2}(ses,1,2,fi));
      else
         ys = min(ylm(mm,ii), pmfdat{mm,3}(ses,ii-1,fi));
         es = max(0.01,       pmfdat{mm,3}(ses,3   ,fi));
      end
      %ys = min(ylm(mm,ii), pmfdat{mm,2}(ses,inds(ii,1),inds(ii,2),1));
      %es = max(0.01, pmfdat{mm,2}(ses,inds(ii,1),inds(ii,2),2));
      ye = [ys-es ys+es]';
      % plot(sax, ye, '-', 'Color', 0.8*ones(1,3));
      % plot(ses, ys, '.', 'Color', 0.6*ones(1,3));
      plot(sax(:, Llapse), ye(:, Llapse), '-', 'Color', 0.8*ones(1,3));
      plot(ses(Llapse), ys(Llapse), 'k.');
      axis([1 ses(end) 0 ylm(mm,ii)]);
      Lf = Llapse & isfinite(ys) & isfinite(es);
      if ii < 3
         [bw,bwint] = regressW_mike(ys(Lf), es(Lf), A(Lf,:));
         [b,bint]   = regress(ys(Lf), A(Lf,:));
         if sign(bwint(1,1)) == sign(bint(1,2))
            plot(ses, A*b, 'k--', 'LineWidth', 2);
         end
      else
         f = getFIT_exp1([A(Lf,1) ys(Lf) es(Lf)]);
         plot(ses, valFIT_exp1(f, ses), 'k--', 'LineWidth', 2)
         disp(sprintf('mm=%d,tau=%.2f',mm,f(1)))
      end
      axis([1 ses(end) 0 ylm(mm,ii)]);
   end
end


return

%% OLD STUFF

pmfdat = FS_loadProjectFile('2008_Bias', 'figBIAS_pmfParams');

if isempty(pmfdat)
   
   pmfdat = cell(monkn, 4);
   fcs    = getBIAS_fcs;
   funs   = {@ddExp2_L; @ddExp4fz_LF};
   nfuns  = length(funs);
   
   for mm = 1:monkn
      
      dat          = FS_getDotsTrainingData(monks{mm});
      Lgood        = dat(:,2) <= 2 & dat(:,3)>=0 & isfinite(fcs{mm});% & dat(:,6)<0.8;
      L51          = dat(:,5) == 0.5120;
      L99          = dat(:,5) >  0.9;
      Lcor         = dat(:,3) == 1;
      sessions     = unique(dat(:,1));
      num_sessions = length(sessions);
      chc          = dat(:,[9 9]);
      chc(~Lcor,1) = 0;
      chc( Lcor,2) = 0;
      fdat         = {[], fcs{mm}, chc, chc};
      
      pmfdat{mm,1} = nans(num_sessions, 3);           % corrcoeff, CIs
      pmfdat{mm,2} = nans(num_sessions, 7, 2, nfuns); % fits/sems; last arg is thresh
      pmfdat{mm,3} = nans(num_sessions, 3, nfuns);    % bs mean/sem, bs/dv
      pmfdat{mm,4} = nans(size(dat,1),  3, nfuns);    % preds, resids, biases
      
      for ss = 1:num_sessions
         
         Lses  = Lgood & dat(:,1) == sessions(ss);
         disp([ss sessions(ss) sum(Lses)])
         
         % check 99 & 51 coh
         Ltm = Lses & dat(:,6) > nanmedian(dat(Lses,6));
         if sum(Ltm&L99&Lcor)/sum(Ltm&L99) < sum(Ltm&L51&Lcor)/sum(Ltm&L51)
            Lses = Lses & ~L99;
         end
         
         if sum(Lses) > 100
            
            % Get lapse rate
            lapse  = getSLOPE_lapse(dat(Lses, [3 5 6]));
            ltofit = min(0.2, max(lapse, 0.001));
            
            % fit to simple model
            [f1,s1,t1,p1,r1] = ctPsych_fit(funs{1}, dat(Lses, [5 6 8 9]), ...
               dat(Lses, 3), [], [], [], [], ltofit);
            % ctPsych_plot(fun1, f1, dat(Lses, [5 6 8 9]), dat(Lses, 3))
            
            pmfdat{mm,2}(ss,1:7,1,1) = [f1' nans(1,3) lapse ctPsych_thresh(funs{1}, f1)];
            pmfdat{mm,2}(ss,1:2,2,1) = s1;
            pmfdat{mm,4}(Lses,1:2,1) = [p1, r1(:,4)];
            
            % correlation coef between wk-filtered choices and resids
            [R, P, RLO, RUP]    = corrcoef(fcs{mm}(Lses), r1(:,4));
            pmfdat{mm,1}(ss, :) = [R(2,1) RLO(2,1) RUP(2,1)];
            
            % fit to funs, using A and lapse from ddExp3
            for ff = 2:nfuns
               if ff == 2
                  va = [];
               else
                  va = fdat{ff-1}(Lses,:);
               end
               [fits,sems,ts,ps,rs,d]  = ctPsych_fit(funs{ff}, dat(Lses,[5 6 8 9]), ...
                  dat(Lses, 3), [], [], f1(1), [], ltofit, va);
               bs                      = feval(funs{ff}, fits, d, nan, va);
               pmfdat{mm,2}(ss,:,1,ff) = [fits' nans(1,7-length(fits)-1) ctPsych_thresh(funs{ff}, fits)];
               pmfdat{mm,2}(ss,:,2,ff) = [sems' nans(1,7-length(fits))];
               pmfdat{mm,3}(ss,:,ff)   = [nanmean(abs(bs(:,1))), ...
                  nanmean(abs(bs(:,1)))./nanmean(abs(bs(:,2))), nanse(abs(bs(:,1)))];
               pmfdat{mm,4}(Lses,:,ff) = [ps rs(:,4) bs(:,1)];
            end
         end
      end
   end
   
   % save it to disk
   FS_saveProjectFile('2008_Bias', 'figBIAS_pmfParams', pmfdat);
end

% for mm = 1:4
%     ses    = (1:mse(mm))';
%     Llapse = pmfdat{mm,2}(ses,4,1,1) < 0.1;
%     ses    = ses(Llapse);
%     t1     = pmfdat{mm,2}(ses,5,1,1);
%     t2     = pmfdat{mm,2}(ses,5,2,1);
%     ts     = t1./t2;
%     Lf     = isfinite(t1) & isfinite(t2);
%     disp([nanmedian(ts) iqr(ts(Lf)) signrank(t1(Lf), t2(Lf))])
% %    disp([nanmean(pmfdat{mm,2}(ses,1,1,1)./pmfdat{mm,2}(ses,1,2,1)) ...
% %        nanmedian(pmfdat{mm,2}(ses,5,2,1)./pmfdat{mm,2}(ses,5,1,1))])
% end
%

%% piecewise linear
%         f = getFIT_pwl2([A(Lf,1) ys(Lf) es(Lf)]);
%         plot(ses, valFIT_pwl2(f,ses), 'g--');
%         disp(sprintf('mm=%d,ii=%d: s=%.2f, c=%d', mm, ii, ...
%             f(2), round(f(1)*ses(find(Lf,1,'last')))))

%         if ii == 1
%             [f,s] = getFIT_exp1U([A(Lf,1) ys(Lf) es(Lf)]);
%             plot(ses, valFIT_exp1U(f, ses), 'r--')
%         else
%             [f,s] = getFIT_exp1([A(Lf,1) ys(Lf) es(Lf)]);
%             plot(ses, valFIT_exp1(f, ses), 'r--')
%         end
%       [bw,bwint] = regressW(ys(Lf), es(Lf), A(Lf,:));
%       [b,bint]   = regress(ys(Lf), A(Lf,:));
%       if sign(bwint(1,1)) == sign(bint(1,2))
%           plot(ses, A*b, 'k--', 'LineWidth', 2);
%       end


%% OLD STUFF BELOW
%%
% for each monkey
fi    = 2;
li    = [5 3 1]; %  bias, lapse, A (no local bias=4)
lims  = [0 1; 0 .2; 0 75]; % 0 10 for lb
nfuns = size(pmfdat{1,2},3)-1;
for mm = 1:monkn
   
   % lbias vs thresh
   for ll = 1:length(li)
      axes(axs((mm-1)*(length(li)+1)+ll)); cla reset; hold on;
      xs = pmfdat{mm,1}(:,1);
      if li(ll) == 5
         ys = abs(pmfdat{mm,2}(:,li(ll),fi));
      else
         ys = pmfdat{mm,2}(:,li(ll),fi);
      end
      plot(xs,ys,'k.');
      Lgood = isfinite(xs) & isfinite(ys);
      [R,P] = corrcoef(xs(Lgood), ys(Lgood));
      if P(2,1)<0.05
         title(sprintf('%.2f*', R(2,1)))
         h=lsline;
         set(h, 'LineWidth', 2, 'LineStyle', '--');
      else
         title(sprintf('%.2f', R(2,1)))
      end
      plot([-0.1 0.42], [0 0], 'k:');
      plot([0 0], lims(ll,:), 'k:');
      axis([-0.1 0.42 lims(ll,:)]);
   end
   
   % thresh (ddExp3) vs thresh (ddExp5fzz)
   axes(axs((mm-1)*(length(li)+1)+(length(li)+1))); cla reset; hold on;
   plot([0 1], [0 1], 'k:');
   plot(pmfdat{mm,2}(:,6,1), pmfdat{mm,2}(:,6,fi), 'k.')
   axis([0 1 0 1]);
   if mm == 4
      xlabel('Threshold');
   end
   
   % display geometric mean threshold fractions
   tfs = nans(1,nfuns);
   for ff = 1:nfuns
      vals = pmfdat{mm,2}(:,6,ff+1)./pmfdat{mm,2}(:,6,1);
      tfs(ff) = geomean(vals(isfinite(vals)));
   end
   disp(tfs)
end
axes(axs(13));
ylabel('Bias (zscore)');
xlabel('r_behav');
axes(axs(14));
ylabel('lapse');
xlabel('r_behav');
axes(axs(15));
ylabel('A');
xlabel('r_behav');

