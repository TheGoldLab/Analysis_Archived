function fig_ = figBIAS_pmfCompare(num)
% function fig_ = figBIAS_pmfCompare(num)
%

if nargin < 1 || isempty(num)
   num = 7;
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
pmfdat = FS_loadProjectFile('2008_Bias', 'figBIAS_pmfCompare');

if isempty(pmfdat)
   
   pmfdat = cell(monkn, 5);
   fcs    = getBIAS_fcs;
   funs   = {@ddExp2_L; @ddExp3z_L; @ddExp41B_LF; @ddExp5RL_LF; @ddExp3fz_LF; @ddExp4fz_LF; @ddExp6RL_LF};
   nfuns  = length(funs);
   acsz   = 100;
   
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
      
      pmfdat{mm,1} = nans(num_sessions, 3);           % corrcoeff, CIs
      pmfdat{mm,2} = nans(num_sessions, 7, 2, nfuns); % fits/sems; last arg is thresh
      pmfdat{mm,3} = nans(num_sessions, 3, nfuns);    % bs mean/sem, bs/dv
      pmfdat{mm,4} = nans(size(dat,1),  3, nfuns);    % preds, resids, biases
      pmfdat{mm,5} = nans(num_sessions, 3, 2, nfuns); % mean/median/sem; resid/chc/dir; DC/AC
      
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
            lapse  = getBIAS_lapse(dat(Lses, [3 5 6]));
            ltofit = min(0.2, max(lapse, 0.001));
            
            % fit to simple model
            [f1,s1,t1,p1,r1] = ctPsych_fit(funs{1}, dat(Lses, [5 6 8 9]), ...
               dat(Lses, 3), [], [], [], [], false, ltofit);
            % ctPsych_plot(fun1, f1, dat(Lses, [5 6 8 9]), dat(Lses, 3))
            
            pmfdat{mm,2}(ss,1:7,1,1) = [f1' nans(1,3) lapse ctPsych_thresh(funs{1}, f1)];
            pmfdat{mm,2}(ss,1:2,2,1) = s1;
            pmfdat{mm,4}(Lses,1:2,1) = [p1, r1(:,4)];
            
            pmfdat{mm,5}(ss,:,1,1) = [nanmean(r1(:,4)) nanmean(r1(:,4).^2) nanse(r1(:,4))];
            ac = xcorr(r1(:,4), acsz, 'coeff');
            pmfdat{mm,5}(ss,:,2,1) = [ac(acsz+2) abs(mean(ac(acsz+(2:51)))) ...
               abs(mean(ac(acsz+(2:101))))];
            
            % correlation coef between wk-filtered choices and resids
            [R, P, RLO, RUP]    = corrcoef(fcs{mm}(Lses), r1(:,4));
            pmfdat{mm,1}(ss, :) = [R(2,1) RLO(2,1) RUP(2,1)];
            
            % fit to funs, using A and lapse from ddExp3
            for ff = 2:nfuns
               if ff == 2
                  va = {};
               elseif ff==3 || ff==4 || ff==7
                  va = {chc(Lses,:)};
               elseif ff==5 || ff==6
                  va = {fcs{mm}(Lses)};
               end
               % [fits,sems,ts,ps,rs,d]  = ctPsych_fit(funs{ff}, dat(Lses,[5 6 8 9]), ...
               %     dat(Lses, 3), [], [], f1(1), [], ltofit, va);
               % updated by jig 1/5/18
               [fits,sems,~,ps,rs,d]   = ctPsych_fit(funs{ff}, dat(Lses,[5 6 8 9]), ...
                  dat(Lses, 3), [], [], f1(1), [], true, ltofit, va{:});
               bs                      = feval(funs{ff}, fits, d, nan, va{:});
               pmfdat{mm,2}(ss,:,1,ff) = [fits' nans(1,7-length(fits)-1) ctPsych_thresh(funs{ff}, fits)];
               pmfdat{mm,2}(ss,:,2,ff) = [sems' nans(1,7-length(fits))];
               nma                     = nanmean(abs(bs(:,1)));
               pmfdat{mm,3}(ss,:,ff)   = [nma, nma./nanmean(abs(bs(:,2))), nanse(abs(bs(:,1)))];
               pmfdat{mm,4}(Lses,:,ff) = [ps rs(:,4) bs(:,1)];
               
               pmfdat{mm,5}(ss,:,1,ff) = [nanmean(rs(:,4)) nanmean(rs(:,4).^2) nanse(rs(:,4))];
               ac = xcorr(rs(:,4), acsz, 'coeff');
               pmfdat{mm,5}(ss,:,2,ff) = [ac(acsz+2) abs(mean(ac(acsz+(2:51)))) ...
                  abs(mean(ac(acsz+(2:101))))];
            end
         end
      end
   end
   
   % save it to disk
   FS_saveProjectFile('2008_Bias', 'figBIAS_pmfCompare', pmfdat);
end

rsdat = FS_loadProjectFile('2008_Bias', 'figBIAS_residSummary');
c3 = ones(1,3);
co = linspace(0.8,0,5)'*c3;
sy = {'>' 'v' 'd' '^' 's'};
fi = 2:6;
ms = 7;
for mm = 1:monkn

   % ses = (1:mse(mm))';
   ses = (1:size(pmfdat{mm,1},1))';
   
   % Boxplots of residual fraction
   axes(axs(mm)); cla reset; hold on;
   mdat = [];
   for ff = 1:5
      mdat = cat(1, mdat, [pmfdat{mm,5}(ses,2,1,fi(ff))./pmfdat{mm,5}(ses,2,1,1), ...
         (ff)*ones(length(ses),1)]);
   end
   h=boxplot(mdat(:,1), mdat(:,2), 'notch', 'on');
   set(h(end,:), 'Visible', 'off');
   set(h(1:2,:), 'LineStyle', '-');
   for ff = 1:5
      set(h(:,ff),'Color',co(ff,:));
      plot(ff,nanmedian(mdat(mdat(:,2)==ff),1), sy{ff}, 'Color', co(ff,:), ...
         'MarkerFaceColor', co(ff,:), 'MarkerSize', ms);
   end
   ylim([0.77 1.01]);
   
   % weirdo plot
   axes(axs(mm+4)); cla reset; hold on;
   plot([0 0], [0 1], 'k:');
   for ff = 1:5
      xs = abs(pmfdat{mm,5}(ses,1,2,fi(ff)));%.*sign(pmfdat{mm,5}(ses,2,2,fi(ff)));
      %        ys = abs(pmfdat{mm,5}(ses,3,2,fi(ff)));
      ys = abs(pmfdat{mm,5}(ses,1,1,fi(ff)));
      % plot(xs, ys, '.', 'Color', co{ff}, 'MarkerSize', 5);
      pcx = prctile(xs(isfinite(xs)), [25 50 75]);
      pcy = prctile(ys(isfinite(ys)), [25 50 75]);
      plot(pcx([2 2]), pcy([1 3]), '-', 'Color', co(ff,:));
      plot(pcx([1 3]), pcy([2 2]), '-', 'Color', co(ff,:));
      plot(pcx(2),     pcy(2),     sy{ff}, 'Color', co(ff,:), ...
         'MarkerFaceColor', co(ff,:), 'MarkerSize', ms);
   end
   
%    % dir, for ref
%    xs = abs(rsdat{mm,2}(ses,1,2,2));
%    ys = abs(rsdat{mm,2}(ses,2,2,2));
%    pcx = prctile(xs(isfinite(xs)), [25 50 75]);
%    pcy = prctile(ys(isfinite(ys)), [25 50 75]);
%    %    plot(pcx([2 2]), pcy([1 3]), 'k--', 'LineWidth', 1);
%    %    plot(pcx([1 3]), pcy([2 2]), 'k--', 'LineWidth', 1);
%    plot(pcx(2), pcy(2), 'r*', 'MarkerSize', 9);
   
   axis([0 .1 0 .08]);
end
return


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

%% Plotit
%inds = [1 1; 6 2; 8 2];
ylm  = [75 1.5 1.2; 50 1 1.2; 90 1.5 1.2; 30 2 1.2];
for mm = 1:monkn
   
   ses = (1:mse(mm))';
   sax = [ses ses]';
   A   = [ses ones(size(ses))];
   Llapse = pmfdat{mm,2}(ses,6,1,1) < 0.1;
   
   for ii = 1:3
      axes(axs(mm+(ii-1)*4)); cla reset; hold on;
      if ii == 1
         ys = min(ylm(mm,ii), pmfdat{mm,2}(ses,1,1,1));
         es = max(0.01,       pmfdat{mm,2}(ses,1,2,1));
      else
         ys = min(ylm(mm,ii), pmfdat{mm,3}(ses,ii-1,3));
         es = max(0.01,       pmfdat{mm,3}(ses,3   ,3));
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
         [bw,bwint] = regressW(ys(Lf), es(Lf), A(Lf,:));
         [b,bint]   = regress(ys(Lf), A(Lf,:));
         if sign(bwint(1,1)) == sign(bint(1,2))
            plot(ses, A*b, 'k--', 'LineWidth', 2);
         end
      else
         [f,s] = getFIT_exp1([A(Lf,1) ys(Lf) es(Lf)]);
         plot(ses, valFIT_exp1(f, ses), 'k--', 'LineWidth', 2)
         disp(sprintf('mm=%d,tau=%.2f',mm,f(1)))
      end
      
      
      %% piecewise linear
      %         f = getFIT_pwl2([A(Lf,1) ys(Lf) es(Lf)]);
      %         plot(ses, valFIT_pwl2(f,ses), 'g--');
      %         disp(sprintf('mm=%d,ii=%d: s=%.2f, c=%d', mm, ii, ...
      %             f(2), round(f(1)*ses(find(Lf,1,'last')))))
      axis([1 ses(end) 0 ylm(mm,ii)]);
   end
end

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

return

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

