function fitTrainingData
% function fitTrainingData
%
% Fit matched training and testing data to the 
%  Adaptivity and Particle Filter Models
%
% NOTE that the particle Filter fits can take hours to run
%  (2017 MacBook Pro)
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"

%% Find the data
[file_names, data_dir, raw_data_dir] = getDataInfo;
num_subjects                         = length(file_names);

%% Define task values
MU_DIST = 30;

%% Set initial values, lower- and upper-bounds for both models
% -- Adaptivity Model fit parameters
%   1. H_default
%   2. H_slope
%   3. noise in the DV
%   4. lapse rate
models(1).name    = 'Adaptivity';
models(1).func    = @fitAdaptivityModel_err;
models(1).params0 = [-0.2 0.2 0.4 0.004];
models(1).lb      = [-5 -2 0 0];
models(1).ub      = [5 2 10 0.03];
models(1).ic      = 1; % for indexing data
models(1).is      = [1 2 4 5 6];
models(1).options = optimoptions(@fmincon,'Algorithm','interior-point');

% -- Particle Filter Model fit parameters
%   1. H0       ... default H
%   2. logR0    ... precision of prior on H
%   3. logK     ... volatility
%   4. vx_scale ... scale factor on standard normal sensory noise
%   5. lapse rate
models(2).name    = 'ParticleFilter';
models(2).func    = @fitParticleFilterModel_err;
models(2).params0 = [0.4 2 -5 1.5 0.0025];
models(2).lb      = [0.001 log(0.1) log(.0001) 0.0001 0.0001];
models(2).ub      = [0.999 log(2000) log(.1) 5 0.05];
models(2).ic      = 1:2;
models(2).is      = [3 4 6];
models(2).options = psoptimset(...
   'MaxIter',5000,'Display','iter','CompletePoll','on', ...
   'TolMesh',1e-6,'UseParallel',true,'ScaleMesh','on','InitialMeshSize',1, ...
   'TimeLimit',900,'CompleteSearch','off');

% initialize objective function for particle filter model
numberOfParticles = 20;
simN              = 10000;
eta0              = randn(20000,1);
feval(models(2).func, [], [], numberOfParticles, simN, eta0);

% add fields for LLR, fits
num_models = length(models);
for mm = 1:num_models
   models(mm).LLRs   = nans(num_subjects, 2);
   models(mm).params = nans(num_subjects, length(models(mm).params0), 2);
end

%% Loop through the subjects
for ss = 1:num_subjects
   
   % get the data file
   data_filename = fullfile(raw_data_dir, file_names{ss});
   load(data_filename)
   inds         = [];
   num_sessions = session_ind - 1;
   
   for ii = 1:num_sessions
      
      % look for practice data set with matched task parameters
      if exist(['data' num2str(ii) '_practice'], 'var')
         
         eval(['data=data' num2str(ii) ';']);
         eval(['datap=data' num2str(ii) '_practice;']);
         
         uH  = unique(data.H(:,1));
         uHp = unique(datap.H(:,1));
         
         if isequal(uH, uHp) && length(uH)==2 && ...
               mode(data.signaled)==mode(datap.signaled)
            inds = cat(2, inds, ii);
         end
      end
   end
   
   % found at least two matched sequencies
   if length(inds) > 2
      
      % see below for data collected
      session_data = cell(length(inds), 2, 2);
      ind          = 1;
      for ii = inds
         for dd = 1:2
            
            % practice data first
            if dd==1
               eval(['data=data' num2str(ii) '_practice;']);
               swp = find(diff(data.H(:,1))~=0,1);
               trs = 1:data.N;
            else
               eval(['data=data' num2str(ii) ';']);
               swd = find(diff(data.H(:,1))~=0,1);
               trs = trs-swp+swd;
            end
            
            % useful stuff
            midpt = mean(data.muall(:,1));
            musgn = sign(data.muall(data.muinds,1)-midpt);
            musgn(data.signaled==0) = 0;
            musgn = [0;musgn(1:end-1)];
            if data.H(1)>1
               data.H = 1./data.H;
            end
            
            % Collect some data for each session
            % FIRST CELL is per trial, columns are:
            %   1. likelihood of x given left
            %   2. likelihood of x given right
            %   3. x location wrt midpoint
            %   4. choice
            %   5. H (was logJ0 = log[H/(1-H)])
            %   6. Signaled mean from previous trial (if given)
            session_data{ind,1,dd} = cat(2, ...
               normpdf(data.X(:,1)-midpt, MU_DIST/2, data.sigma), ...
               normpdf(data.X(:,1)-midpt,-MU_DIST/2, data.sigma), ...
               data.X(:,1)-midpt, ...
               double(data.pred==2), ...
               log(data.H(:,1)./(1-data.H(:,1))), ... % data.H(:,1)
               musgn);

            % include just a subset of trials
            session_data{ind,1,dd} = session_data{ind,1,dd}(trs,:);

            % SECOND CELL is scalar:
            %   F = sigma^2/mu
            session_data{ind,2,dd} = data.sigma.^2./MU_DIST;            
         end
         ind = ind + 1;
      end
      
      % Do the fits for this subject, per model
      for mm = 1:num_models
         for dd = 1:2
            
            % get arguments for specific model, using weird indices
            args = session_data(:,models(mm).ic,dd);
            for aa = 1:size(args,1)
               args{aa,1} = args{aa,1}(:,models(mm).is);
            end
            
            % do the fit for this subject
            myFun = @(x)models(mm).func(x, args);
            
            % now... fit it
            [models(mm).params(ss,:,dd), models(mm).LLRs(ss,dd)] = ...
               fmincon(myFun, models(mm).params0, [], [], [], [], ...
               models(mm).lb, models(mm).ub, [], models(mm).options);
         end
      end
   end
end

% save 'em
save(fullfile(data_dir, 'trainingDataFits.mat'), 'models');
