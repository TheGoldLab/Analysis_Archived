function e = getfit_Hregress(l1arr,l2arr,choicearr,musgnarr,logJarr,params)

arrN = numel(l1arr);

logJ0 = params(1);
alpha = params(2);
psi = params(3);
whoops = params(4);

evc = zeros(arrN,1);
parfor arrind = 1:arrN
    
    logJvc = logJ0 + alpha*(logJarr{arrind}-logJ0);
    Hvc = 1./(1+exp(-logJvc));
    l1 = l1arr{arrind};
    l2 = l2arr{arrind};
    [q1,q2] = dsprt_mixed_Hvc2(l1,l2,Hvc,musgnarr{arrind});
    L = log(q1)-log(q2);
    
    choicehat = 1./(1+exp(-L/psi));
    choicehat = (1-whoops)*choicehat+whoops*(1-choicehat);
    
    choicehat = choicehat(:);
    choice = choicearr{arrind};
    evc(arrind) = -sum((1-choice).*log(1-choicehat)) - sum(choice.*log(choicehat));
    
end
e = sum(evc);
% 
% e = e - lognormpdf(logJ0,0,2);