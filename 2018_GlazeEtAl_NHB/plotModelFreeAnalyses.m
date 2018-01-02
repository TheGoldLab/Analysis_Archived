function plotModelFreeAnalyses(dataType, axes_vs_trials, axes_vs_evidence)
% function plotModelFreeAnalyses(dataType, axes_vs_trials, axes_vs_evidence)
%
% dataType is 'data' or 'simulation'
% other arguments are vectors of axis handles
%
% Utility function to plot "model-free" analyses of real and
%  simulated data
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line
%     learning in an unpredictable environment"

% Where to find the data
[file_list, analysis_data_dir, raw_data_dir] = getDataInfo;

%% Get adaptivity model fits.
% These are generated using: fitAdaptivityModel
% Parameters are:
%  1: logJ0 = log(Hdefault/(1?Hdefault)
%  2: adaptivity (mH)
%  3: decision noise
%  4: lapse
load(fullfile(analysis_data_dir, 'adaptivityModelFits.mat'), 'fits');

% for simulations, load particle-filter fits
if ~strcmp(dataType, 'data')
   load(fullfile(analysis_data_dir, ...
      'allparam_fixedvarfits_partfilt_sigma_partMfixed20.mat'),...
      'paramstrct_sigma_partMfixed');
   part20fits = [paramstrct_sigma_partMfixed.params]';
end

%% Some constants
MU_DIST         = 30;   % distance between targets
INITIAL_CUTOFF  = 2.5;  % LogLR of initial evidence for switch
EVIDENCE_CUTOFF = 2.5;  % LogLR of subsequent evidence for switch
MIN_TRIALS      = 100;  % per block
HAZARD_BINS     = [0 .2 .8 1]; % Binned by hazard rate
RUN_LENGTH_BINS = [0 1 2 3 4 inf]; % Show sequences of runs after change-points
ADAPTIVITY_BINS = [-inf .1 .3 inf]; % Adaptivity bins
EVIDENCE_BINS   = [-30 -5:1:5 30];

% some useful stuff
num_hazard_bins     = length(HAZARD_BINS)-1;
num_run_length_bins = length(RUN_LENGTH_BINS)-1;
num_adaptivity_bins = length(ADAPTIVITY_BINS)-1;
num_evidence_bins   = length(EVIDENCE_BINS)-1;

%% Loop through the subjects
% collecting data -- see below for details
num_subjects = length(file_list);
all_data     = [];
for ss = 1:num_subjects
   
   % get the data file
   data_filename = fullfile(raw_data_dir, file_list{ss});
   load(data_filename);
   
   % Loop through the sessions
   num_sessions = session_ind - 1;
   for ii = 1:num_sessions
      
      % get the appropriate data
      eval(['data=data' num2str(ii) ';'])
      
      % useful stuff
      midpt = mean(data.muall(:,1));
      musgn = sign(data.muall(data.muinds,1)-midpt);
      musgn(data.signaled==0) = 0;
      musgn = [0;musgn(1:end-1)];
      if data.H(1)>1
         data.H = 1./data.H;
      end
      
      if strcmp(dataType, 'data')

         % REAL DATA
         % Was the choice correct or not? needed to compute switch...
         correct_vect = data.pred==data.muinds;
         
         % Did the subject switch choices?
         switch_vect = [0; abs(diff(data.pred))];     
      else
         
         % SIMULATIONS         
         % collect useful parameters from fits, etc
         params  = part20fits(ss,:);
         x       = data.X(:,1)-midpt;
         F       = data.sigma^2/30;
         H0samps = betarnd( ...
            exp(params(6))*params(7),exp(params(6))*(1-params(7)),20000,1);
         eta     = params(2)*randn(2000,1);
         
         % simulate choices
         [Lexp,~] = particle_filter_learnH6( ...
            exp(params(8)),H0samps,F,x,musgn,params(1),1,eta);
         
         % compute choices, accounting for lapses
         choice = double(Lexp>0);
         whoops = logical(rand(2000,1)<params(3));
         choice(whoops) = 1 - choice(whoops);
         correctchoice = double(data.muinds==2);
         correct_vect = double(choice==correctchoice);
         
         % switch?
         switch_vect = [0; abs(diff(choice))];
      end

      % get LLR wrt right triangle
      LLR_right = (data.X(:,1)-midpt) * MU_DIST / data.sigma^2;
      
      % get LLR wrt switch
      LprevRight = (data.muinds==2 & data.cp==0) | ...
         (data.muinds==1 & data.cp==1);
      LLR_switch = LLR_right;
      LLR_switch(LprevRight) = -LLR_switch(LprevRight);
      
      % get LLR wrt correct side
      LLR_correct = LLR_right;
      LLR_correct(data.muinds==1) = -LLR_correct(data.muinds==1);
      
      % indices of change-point trials (except first)
      cp_trial_inds = [setdiff(find(data.cp), 1); data.n];
      
      % make binary array of indices of trials before and after
      %  change_points following strong evidence
      cp_array = zeros(data.n, 2);
      for cc = 1:length(cp_trial_inds)-1
         if LLR_correct(cp_trial_inds(cc)-1)>=INITIAL_CUTOFF
            cp_array(cp_trial_inds(cc)-1,1) = 1;
            cp_array(cp_trial_inds(cc):cp_trial_inds(cc+1)-1,2) = 1;
         end
      end
      
      % H-block indices
      H_blocks = [1; find(diff(data.H(:,1))); data.n+1];
      H_inds   = nans(data.n, 1);
      for hh = 1:length(H_blocks)-1
         H_inds(H_blocks(hh):H_blocks(hh+1)-1) = ...
            H_blocks(hh):H_blocks(hh+1)-1;
      end
      
      % probability of switching to the correct side
      % collect data from this session, columns are:
      %  1. trial number in block
      %  2. trial before "good" switch
      %  3. trials after "good" switch
      %  4. prevous LLR correct
      %  5. current LLR for switch
      %  6. run length
      %  7. correct?
      %  8. switch?
      %  9. H
      %  10. Signaled mean from previous trial (if given)
      %  11. adaptivity (mH from adaptivity fits per subject)
      all_data = cat(1, all_data, [ ...
         H_inds, ...
         cp_array, ...
         [0; LLR_correct(1:end-1)], ...
         LLR_switch, ...
         data.r, ...
         correct_vect, ...
         switch_vect, ...
         data.H(:,1), ...
         musgn, ...
         fits(ss,2)*ones(data.n,1)]);
   end
end

%% FIRST COMPUTE VS TRIAL
if nargin > 1 && ~isempty(axes_vs_trials)
   
   % make data matrix
   pswitch_vs_trial_data = nans(num_run_length_bins, 3, num_hazard_bins, num_adaptivity_bins);
   
   % select only trials with no signaled answer from previous trial
   Lgood = all_data(:,10)==0;
   
   % for each adaptivity bin
   for aa = 1:num_adaptivity_bins
      
      % select appropriate adaptivities
      Ladaptivity = Lgood & all_data(:,11)>=ADAPTIVITY_BINS(aa) & ...
         all_data(:,11)<ADAPTIVITY_BINS(aa+1);
      
      % for each hazard bin
      for hh = 1:num_hazard_bins
         
         % select appropriate hazards
         Lhazard = Ladaptivity & ...
            all_data(:,9)>=HAZARD_BINS(hh) & ...
            all_data(:,9)<HAZARD_BINS(hh+1);
         
         % for each run length
         for rr = 1:num_run_length_bins
            
            if rr==1 % first just checks for previous trial was an error
               tmp_array = 1-all_data(Lhazard&all_data(:,2)==1,7);
            else
               % get good trials
               Larr = Lhazard & all_data(:,3)==1 & ...
                  abs(all_data(:,5))<=EVIDENCE_CUTOFF & ...
                  all_data(:,6)>=RUN_LENGTH_BINS(rr) & ...
                  all_data(:,6)< RUN_LENGTH_BINS(rr+1);
               tmp_array = all_data(Larr,7);
            end
            
            if length(tmp_array) > 12
               pswitch_vs_trial_data(rr, :, hh, aa) = [ ...
                  mean(tmp_array), std(bootstrp(200,@mean, tmp_array)), ...
                  length(tmp_array)];
            end
         end
      end
   end
   
   for aa = 1:num_adaptivity_bins
      axes(axes_vs_trials(aa)); cla reset; hold on;
      
      errorbar( ...
         repmat(1:num_run_length_bins,num_hazard_bins,1)', ...
         squeeze(pswitch_vs_trial_data(:,1,:,aa)), ...
         squeeze(pswitch_vs_trial_data(:,2,:,aa)), ...
         '-','linewidth',1.5);
      
      axis([0.5 num_run_length_bins+0.5 -0.05 1.05])
      set(gca,'box','off','ticklength',[0 0], ...
         'fontsize',12,'xtick',1:num_run_length_bins,'xticklabel',{'pre','1','2','3','4+'})
   end
end

%% SECOND COMPUTE VS EVIDENCE
if nargin > 2 && ~isempty(axes_vs_evidence)
   
   % make data matrix
   pswitch_vs_evidence_data = nans(num_evidence_bins, 4, num_hazard_bins, num_adaptivity_bins);
   
   % select trials with feedback or strong evidence on previous trial, after
   % first 100 trials of block
   Lgood = all_data(:,1) > MIN_TRIALS & (all_data(:,4)>=EVIDENCE_CUTOFF | all_data(:,10)~=0);
   
   % for each adaptivity bin
   for aa = 1:num_adaptivity_bins
      
      % select appropriate adaptivities
      Ladaptivity = Lgood & all_data(:,11)>=ADAPTIVITY_BINS(aa) & ...
         all_data(:,11)<ADAPTIVITY_BINS(aa+1);
      
      % for each hazard bin
      for hh = 1:num_hazard_bins
         
         % select appropriate hazards
         Lhazard = Ladaptivity & ...
            all_data(:,9)>=HAZARD_BINS(hh) & ...
            all_data(:,9)<HAZARD_BINS(hh+1);
         
         % for each evidence bin
         for ee = 1:num_evidence_bins
            Levidence = Lhazard & all_data(:,5)>=EVIDENCE_BINS(ee) & ...
               all_data(:,5)<EVIDENCE_BINS(ee+1);
            
            if sum(Levidence) > 10
               pswitch_vs_evidence_data(ee,:,hh,aa) = [ ...
                  mean(all_data(Levidence,5)), ...
                  mean(all_data(Levidence,8)), ...
                  std(bootstrp(200,@mean, all_data(Levidence,8))), ...
                  sum(Levidence)];
            end
         end
      end
   end
   
   %% Plot it, with logistic fits
   fits = nans(3,3,3);
   LLRs = nans(3,3);
   ops  = optimoptions(@fmincon,'Algorithm','interior-point', 'Display', 'off');
   xax  = (-10:0.1:10)';
   xaxv = [ones(size(xax)) xax];
   for aa = 1:num_adaptivity_bins
      axes(axes_vs_evidence(aa)); cla reset; hold on;
      
      % collect the data
      xs = squeeze(pswitch_vs_evidence_data(:,1,:,aa));
      ys = squeeze(pswitch_vs_evidence_data(:,2,:,aa));
      es = squeeze(pswitch_vs_evidence_data(:,3,:,aa));
      ns = squeeze(pswitch_vs_evidence_data(:,4,:,aa));
      
      % plot as errorbar
      hs=errorbar(xs,ys,es,'o','linewidth',1.5);
      
      % Fit each hazard-specific set to a logistic function
      for hh = 1:3
         
         % get the fits using fitLogistic_err
         xv = [ones(size(xs,1),1) xs(:,hh)];
         yv = [ys(:,hh) ns(:,hh)];
         myFun = @(fits)fitLogistic_err(fits, xv, yv);
         [fits(ii,:,hh), LLRs(ii,hh)] = fmincon(myFun, [(xv\(6*yv(:,1)-3))' 0], [], [], [], ...
            [], [-100 -100 -100], [100 100 100], [], ops);
         
         % plot the fits
         plot(xax, fits(ii,3,hh)+(1-2.*fits(ii,3,hh)).*1./(1+exp(-(xaxv*fits(ii,1:2,hh)'))), ...
            '-', 'LineWidth', 1.5, 'Color', get(hs(hh), 'Color'));
         set(hs(hh), 'MarkerFaceColor', get(hs(hh), 'Color'));
      end
      axis([-11 11 -.05 1.05]);
      set(gca,'box','off','fontsize',12)
   end
   
   %% Logistic fits comparing paired adaptivity groups
   % per hazard
   for hh = 1:3
      
      % two pairwise comparisons:
      %  1. low vs mid adaptivity
      %  2. mid vs high adaptivity
      for xx = 1:2
         
         % get y values (p_switch) and ns for the two conditions
         yv = reshape(permute( ...
            pswitch_vs_evidence_data(:,[2 4],hh,xx:xx+1),[1 4 2 3]),[],2,1);
         
         % Full model
         xv = cat(2, ...
            [ones(num_evidence_bins,1); zeros(num_evidence_bins,1)], ...
            [zeros(num_evidence_bins,1); ones(num_evidence_bins,1)], ...
            [pswitch_vs_evidence_data(:,1,hh,xx); zeros(num_evidence_bins,1)], ...
            [zeros(num_evidence_bins,1); pswitch_vs_evidence_data(:,1,hh,xx+1)]);
         myFun = @(fits)fitLogistic_err(fits, xv, yv);
         [ffF,errF] = fmincon(myFun, [(xv\(6*yv(:,1)-3))' 0], [], [], [], ...
            [], -100.*ones(size(xv,2)+1,1), 100.*ones(size(xv,2)+1,1), [], ops);
         
         % Reduced model: one slope
         xv = cat(2, ...
            [ones(num_evidence_bins,1); zeros(num_evidence_bins,1)], ...
            [zeros(num_evidence_bins,1); ones(num_evidence_bins,1)], ...
            reshape(pswitch_vs_evidence_data(:,1,hh,xx:xx+1),[],1));
         myFun = @(fits)fitLogistic_err(fits, xv, yv);
         [ffR,errR] = fmincon(myFun, [(xv\(6*yv(:,1)-3))' 0], [], [], [], ...
            [], -100.*ones(size(xv,2)+1,1), 100.*ones(size(xv,2)+1,1), [], ops);
         
         % report the fits (slopes) and some stats
         P = 1-chi2cdf(2*(errR-errF),length(ffF)-length(ffR));
         disp(sprintf('h=%d,xx=%d: P=%.5f, slope1=%.2f, slope2=%.2f', hh,xx,P,ffF(3),ffF(4)))
      end
   end
end
