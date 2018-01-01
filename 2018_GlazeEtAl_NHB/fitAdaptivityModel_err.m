function err_ = fitAdaptivityModel_err(params, session_data)
% function err_ = fitAdaptivityModel_err(params, session_data)
%
% params are:
%   1. J_default
%   2. J_slope
%   3. noise in the DV
%   4. lapse rate
%
% session_data is cell array (per session), each cell has a single
%   matrix with rows as trials, columns are:
%   1. likelihood of x given left
%   2. likelihood of x given right
%   3. choice
%   4. logJ0 = log[H/(1-H)]
%   5. Signaled mean from previous trial (if given)

% fit parameters
J_default     = params(1);
J_slope       = params(2);
eta           = params(3);
lapse         = params(4);
num_sessions  = size(session_data,1);
evc           = zeros(num_sessions,1);

% compute separately for each session
for ii = 1:num_sessions
    
    % Subjective hazard is computed as a linear function of objective hazard,
    %    expressed as log-prior-odds (J)
    J   = J_default + J_slope*session_data{ii}(:,4);
    Hvc = 1./(1+exp(-J));
    
    % get posteriors
    [q1,q2] = dsprt_mixed_Hvc2(session_data{ii}(:,1), ...
        session_data{ii}(:,2), Hvc, session_data{ii}(:,5));

     % L = Log-posterior-odds
    L = log(q1)-log(q2);
    
    % LogLR from predicted and actual choices
    choicehat = lapse+(1-2*lapse)./(1+exp(-L/eta));
    choice    = session_data{ii}(:,3);
    evc(ii)   = -sum((1-choice).*log(1-choicehat)) - sum(choice.*log(choicehat));
end
err_ = sum(evc);
