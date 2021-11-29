function vals_ = GSDVses(fits, data, ses)
% function vals_ = GSDVses(fits, data, spec)
% 
% Fit behavior to Gold and Shadlen model to extract decision variables,
% assuming that the decision variable, a, is a function of session.
%
%
%   Computed at values in "data":
%       data(1)   ... coh [0 ... 1]
%       data(2)   ... viewing time (in second)
%
%
%   Given parameters in "fits":
%       the only parameters are m and b, given by A=m*session+b
% 
%

% Modified by jcl on 2/2/07 from quickTs.m
% Copyright 2007 by Joshua I. Gold
%   University of Pennsylvania

persistent func

% return initial values (init, min, max)
if nargin < 2 || isempty(fits)
    % check data matrix.. should include columns:
    %   Coherence
    if isempty(data)
        error('quickTs: no data')
    end
        
    % build function and update guesses
    % A-term
    A     = '(fits(1)+fits(2).*ses)';
    vals_ = [0 0 100; 0.5 0 1];   
    
    % Coherence
    C     = 'data(:,1)';
    
    T     = 'data(:,2)';
    
    % Fano factor
    f     = '0.3';
    
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
