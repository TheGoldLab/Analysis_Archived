function [output]= fitLinMod_wHet (input)

% goal: fit a weighted regression where weights depend linearly on some
% factor.

% inputs:

% input.yVar   = y variable... for example perceptual errors
% input.xVar   = x matrix... for example including intercept, PE, and modulated PE terms
% input.resVar = variable on which residuals depend... e.g. prediction errors
% input.fitSlope      = if true, we fit the slope, if false, we don't

% output.coefficients  = one coefficient for each column of xVar
% output.weightCoeffs  = intercept and linear coefficient of weighting function
% output.weights       = weights given to each point in the fitting process


% params:
% 1) = intercept of weighting function
% 2) slope of weighting function
% 3-inf) regression coefficients

% limits:

nParams=size(input.xVar, 2)
lb=[ .5, -10e2 ];
ub=[ 10e10, 10e2 ];


%% Call fmincon to evaluate squred error in fitting a curve
%  with a parameters initially set to starting point.
input.start_point = [100, 0];

wp=true(size(lb))
if ~input.fitSlope % if we are not fitting the slope, set it to zero.
    wp(2)=false;
end

model = @linMod_wHetFun;                % this is the function handle to the function that takes the parameters and outputs the thing we want to minimize

    oldOpts = optimset('fmincon');
    options=optimset(oldOpts, 'maxFunEvals', 1000000000000, 'MaxIter', 10e10);
    [estimates, sse, ef, o, l, g, h] = fmincon(model, input.start_point, [], [], [], [], lb(wp), ub(wp), [], options);
    
    % output.coefficients  = one coefficient for each column of xVar
    % output.weightCoeffs  = intercept and linear coefficient of weighting function
    % output.weights       = weights given to each point in the fitting process
    % output.fitSlope      = if true, we fit the slope, if false, we don't
    
    
    % get info about fit to spit out:
    [negLogLike, weights, modCoeffs]=linMod_wHetFun(estimates);
    tParams=input.start_point;
    tParams(wp)=estimates;
    output=struct;
    output.weightCoeffs=tParams;
    output.coefficients=modCoeffs;
    output.weights=weights;
    output.logLike=negLogLike.*-1;




% frugFun accepts some inputs (params) and an array that tells it which
% inputs they are (wp)
    function [negLogLike, wt, modCoeffs] = linMod_wHetFun(params)                %input comes from fmincon (initially start_point)
    
        tParams=input.start_point;
        tParams(wp)=params;
        tParams=real(tParams);
        tParams=max([tParams; lb]);
        tParams=min([tParams; ub]);
        
        
        %COMPUTE THE WEIGHTS:
        wt=tParams(1)+abs(input.resVar).*tParams(2);
        % ok, we wont let weights go below zero. this is fudgy.
        wt=max([zeros(size(wt))+.5, wt],[], 2);
        % keyboard

        [modCoeffs, ~, ~, ~, logLike] = regressW_mike(input.yVar,wt, input.xVar);
    
        negLogLike=-1.*logLike;
        
        if ~isfinite(negLogLike)
           %disp('potential problem')
            % keyboard
        end
        
    end

end

