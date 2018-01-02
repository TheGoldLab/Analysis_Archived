function fitAdaptivityModel(filename)
% function fitAdaptivityModel
%
% Fit testing data to the Adaptivity Model
%
% This takes <1 min on a MacBook Pro 2.8 GHz Intel Core i7
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"

%% Initialize values, params are:
%   1. J_default
%   2. J_slope
%   3. noise in the DV
%   4. lapse rate
MU_DIST = 30;
SIGMAS  = round([7.0711 10.0000 12.2474]);
nsigmas = length(SIGMAS);
ub      = [ 10  10 10 .02];
lb      = [-10 -10  0 .0001];
inits   = [0.45 0.4 0.3 .0025];

%% Get information about data files
[file_list, analysis_data_dir, raw_data_dir] = getDataInfo;

%% Loop through the subjects
num_subjects = length(file_list);
fits         = nans(num_subjects, size(inits, 2));
LLRs         = nans(num_subjects, 1);
pcts         = nans(num_subjects, nsigmas); % per sigma

for ss = 1:num_subjects
   
   disp(ss)
   
   % get the data file
   data_filename = fullfile(raw_data_dir, file_list{ss});
   load(data_filename)
   
   % collect some data for each session:
   %   1. likelihood of x given left
   %   2. likelihood of x given right
   %   3. choice
   %   4. logJ0 = log[H/(1-H)]
   %   5. Signaled mean from previous trial (if given)
   num_sessions = session_ind - 1;
   session_data = cell(num_sessions, 1);
   ctmp = nans(num_sessions, 3); % sigma, # correct, N
   
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
      
      % collect data in a single matrix in the celery
      session_data{ii} = cat(2, ...
         normpdf(data.X(:,1)-midpt, MU_DIST/2, data.sigma), ...
         normpdf(data.X(:,1)-midpt,-MU_DIST/2, data.sigma), ...
         double(data.pred==2), ...
         log(data.H(:,1)./(1-data.H(:,1))), ... % data.H(:,1), ...
         musgn);
      
      % save stuff to compute % correct per sigma
      ctmp(ii,:) = [ ...
         data.sigma, ...
         sum(data.pred==data.muinds), ...
         size(data.pred,1)];
   end
   
   % do the fit for this subject
   myFun = @(x)fitAdaptivityModel_err(x, session_data);
   
   % now... fit it
   [fits(ss,:), LLRs(ss)] = ...
      fmincon(myFun, inits, [], [], [], [], lb, ub, []);
   % optimoptions(@fmincon,'Algorithm','interior-point'));
   
   % compute overall accuracy per sigma
   for gg = 1:3
      Lsig = round(ctmp(:,1)) == SIGMAS(gg);
      if any(Lsig)
         pcts(ss,gg) = sum(ctmp(Lsig,2))./sum(ctmp(Lsig,3)).*100;
      end
   end
end

% Save it
if nargin<1 || (nargin>=1 && ~ischar(filename))
   filename = 'adaptivityModelFits.mat';
end
save(fullfile(analysis_data_dir, filename), 'fits', 'LLRs', 'pcts');
