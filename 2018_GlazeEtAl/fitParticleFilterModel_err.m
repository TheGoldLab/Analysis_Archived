function err_ = fitParticleFilterModel_err(params, session_data, i_numberOfParticles, i_simN, i_eta0, i_H0N)
% function err_ = fitParticleFilterModel_err(params, session_data, i_numberOfParticles, i_simN, i_eta0, i_H0N)
%
% session_data is cell array (per session)
%  first cell is per trial (rows), columns are:
%   1. x location wrt midpoint
%   2. choice
%   3. Signaled mean from previous trial (if given)
%  second cell is scalar:
%   F = sigma^2/mu

persistent numberOfParticles simN eta0 H0N

% reset persistent variables
if nargin > 2 || isempty(numberOfParticles)
    
    if nargin < 3 || isempty(i_numberOfParticles)
        numberOfParticles = 20;
    else
        numberOfParticles = i_numberOfParticles;        
    end
    if nargin < 4 || isempty(i_simN)
        simN = 10000;
    else
        simN = i_simN;
    end    
    if nargin < 5 || isempty(i_eta0)
        eta0 = randn(20000,1);
    else
        eta0 = i_eta0;
    end
    if nargin < 6 || isempty(i_H0N)
        H0N = 10000;
    else
        H0N = i_H0N;
    end
    
    if nargout == 0
        return
    end
end

% free parameters
H0       = params(1);
r0       = exp(params(2));
K        = exp(params(3));
vx_scale = params(4);
lapse    = params(5);

% sample H0 from a beta distribution defined by free parameters r0, H0
H0_samples = betarnd(r0*H0, r0*(1-H0), H0N, 1);

% scale the sensory noise
eta = vx_scale*eta0;

% loop through the sessions, computing Log Likelihood
num_sessions  = size(session_data,1);
LLvc          = zeros(num_sessions,1);

parfor ii = 1:num_sessions
        
    % get the estimate
    [DV, Lvar] = particle_filter_learnH6(...
        K,...
        H0_samples,...
        session_data{ii,2},...
        session_data{ii,1}(:,1),...
        session_data{ii,1}(:,3), ...
        numberOfParticles,...
        simN,...
        eta);

    % get p_choice estimates, corrected for lapses
    choicehat = lapse + (1 - 2*lapse).*(0.5 + 0.5*erf(DV./sqrt(2*max(0.01,Lvar))));

    % compute -sum log-likelihood
    LLvc(ii) = -sum(session_data{ii,1}(:,2).*log(choicehat)+ ...
        (1-session_data{ii,1}(:,2)).*log(1-choicehat)); 
end

err_ = sum(LLvc);
