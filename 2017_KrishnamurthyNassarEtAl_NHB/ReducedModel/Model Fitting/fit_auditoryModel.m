function [output]= fit_auditoryModel(input)


% input = structure with following fields:
% whichParams is a logical array, length=6, specifying which parameters
% you want to fit.

% order for startPoint and whichParams:

      
% 1) hazard rate
% 2) motor std
% 3) width of auditory likelihood distribution

% just in case, i guess
wp=logical(input.whichParams);

% limits:
lb=[ 10e-10 , 10e-10, 10e-10 ];
ub=[1-10e-10,  100,   100];


%% Call fmincon to evaluate squred error in fitting a curve
%  with a parameters initially set to starting point.
start_point = input.startPoint(wp);
model = @audModel;                % this is the function handle to the function that takes the parameters and outputs the thing we want to minimize


if ~isfield(input, 'computeNLLflag')|~input.computeNLLflag
    oldOpts = optimset('fmincon');
    options=optimset(oldOpts, 'maxFunEvals', 100000000, 'MaxIter', 10e10);
    [estimates, sse, ef, o, l, g, h] = fmincon(model, start_point, [], [], [], [], lb(wp), ub(wp), [], options);
    output=struct;
    output.params=estimates;
    output.negLogLike=audModel(estimates);

else
    [negLogLike] = audModel(start_point);
    output.negLogLike=negLogLike;
end


% frugFun accepts some inputs (params) and an array that tells it which
% inputs they are (wp)
    function [negLLik] = audModel(params)                %input comes from fmincon (initially start_point)
       
        input.params(wp)=params;
        input.params(~wp)=input.startPoint(~wp);
        
        
        [negLLik, output] = auditoryModel_gen(input);
        
        if ~isfinite(negLLik)
            keyboard
        end
        
    end

end

