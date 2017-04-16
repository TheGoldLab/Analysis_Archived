function [estimates, estimateserror, sse] = fitEricDeanSigmoidFunction(xdata, ydata, SP)
% Call fminsearch with a random starting point.
start_point = SP;
model = @sigfun;
[estimates, sse, ef, o, l, g, h] = fmincon(model, start_point, [], [], [], [], [.0000000001 .0000000001], [100 1000]);
estimateserror = sqrt(diag(-((-h)^(-1))));
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = sigfun(params)
       
        %% PART 1:  MAke predictions about y's given x's.
sp=params(1);    %% SHAPE PARAMETER... large negative values = step function, small negs are line.
bp=params(2);    %% Bound parameter... bp = 1 means that 0 x = 0y, 180 x = 180 y
x=xdata;   


    miny=1./(1+exp(sp.*(0-90))).*180;
    maxy=1./(1+exp(sp.*(180-90))).*180;
    predY=((1./(1+exp(sp.*(x-90))).*180)  - miny) .*180.*bp./(maxy-miny)+.5.*(180-bp.*180);

  
    
    
    
        %% PART 2:  compute the error... ie how different are predicted Y's
        %% from actual y's... 
        
        
        ErrorVector = ydata-predY;
        
        %% PART 3:  How bad is the error on average (or total)
        
        sse = sum(ErrorVector .^ 2);
        
    end
end