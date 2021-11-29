function vals_ = GSDV(fits, data, spec)
% function vals_ = GSDV(fits, data, spec)
% 
% Fit behavior to Gold and Shadlen model to extract decision variables
%
%
%   Computed at values in "data":
%       data(1)   ... coh [0 ... 1]
%       data(2)   ... OPTIONAL time (sec)
%       data(3)   ... OPTIONAL dot dir (-1/1)
%
%   Given parameters in "fits":
%       1     parameters for coherence dependence term, A
%       0 - 1 parameters for exponent for coherence
%       0 - 1 parameters for exponent for time
%       0 - 1 parameters for fano factor, default = 0.3
% 
% Determined by spec ... see below
%

% Modified by jcl on 2/2/07 from quickTs.m
% Copyright 2007 by Joshua I. Gold
%   University of Pennsylvania

persistent func

% return initial values (init, min, max)
if nargin < 2 || isempty(fits)
    
    % check data matrix.. should include columns:
    %   Coherence
    %   OPTIONAL Time
    %   OPTIONAL Direction
    %   % correct
    %   OPTIONAL n
    if isempty(data)
        error('quickTs: no data')
    end
        
    % build function and update guesses
    % A-term
    A     = 'fits(1)';
    vals_ = [20 0 100];   
    
    % Coherence
    switch spec(1)
        case 0
            C     = 'data(:,1)';
        
        case 1
            C     = 'data(:,1).^fits(2)';
            vals_ = cat(1, vals_, [1 0 3]);
    end
    
    % Time
    fiti = size(vals_, 1) + 1;
    switch spec(2)
        case 0
            T     = 'data(:,2)';
        
        case 1
            T     = sprintf('data(:,2).^fits(%d)', fiti);
            vals_ = cat(1, vals_, [1 0 3]);
    end
    
    % Fano factor
    fiti = size(vals_, 1) + 1;
    switch spec(3)
        case 0
            f     = '0.3';
        
        case 1
            f     = sprintf('fits(%d)', fiti);
            vals_ = cat(1, vals_, [0.3 0 20]);
    end
    
    % baseline firing
    R0 = '10';
    
    
    % create mean and sd for gaussian
    m  = [A '.*' C '.*' T];
    sd = ['(' f '.*(2*' R0 '.*' T '+' m ')).^0.5'];
    
    % make the func
    func = ['vals_= 1-normcdf(0,' m ',' sd ');'];

else
    
    % just eval the string
    eval(func);
end
