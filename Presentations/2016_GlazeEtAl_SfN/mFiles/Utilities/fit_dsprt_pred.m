function [params,fval] = fit_dsprt_pred(LLR,choice)
% FITS CHOICE DATA TO A SINGLE-HAZARD RATE MODEL WITH LIKELIHOOD WEIGHT AS
% A FREE PARAMETER AS WELL

ub = [1 20];
params0 = [.4 2];
lb = [0 0];

model = @fit_dsprt_min;

opts = optimoptions(@fmincon,'Algorithm','interior-point');
ms = MultiStart('UseParallel','always','Display','off');

problem = createOptimProblem('fmincon','objective',...
    model,'x0',params0,'lb',lb,'ub',ub,'options',opts);
[params,fval] = run(ms,problem,4);


    function e = fit_dsprt_min(params)

        H = params(1);
        psi = params(2);

        % C-MEX ROUTINE IMPLEMENTING FAST VERSION
        L = dsprt_c(LLR,H);
        
        % THIS IS FOR FITTING PREDICTIVE BEHAVIOR SO WE'RE GOING TO WRITE
        % CHOICES AS A FUNCTION OF LOG PRIOR PREDICTIVE FOR THAT TRIAL, WITH IS
        % SIMPLY LOG POSTERIOR - LOG LIKELIHOOD RATIO PER BAYES THEOREM IN
        % LOG SPACE
        Psi = L-LLR;
       
        choicehat = .5 + .5*erf(Psi/(sqrt(2)*psi));
        e = -sum((1-choice).*log(max(1-choicehat,eps))) - sum(choice.*log(max(choicehat,eps)));
      
    end

end