function [posteriorDist,  CPP]=combinePriorAndAllLikes_KK(priorMu, priorStd, likeMu, likeStd, haz)

%% This code is based on combinePriorAndLike but extended for use in model fitting.
% In particular, since we do not know the particular likelihood
% representation, we'd like to consider all possible internal
% representations.


% IDEA:
% we need this quantity summed over all lambda in [0 180]
% [P(O_hat |S=0,mu_hat,sig_hat,lambda)*(1-pCha) + pCha*P(O_hat
% |S=1,mu_hat,sig_hat,lambda)]*p(lambda |O)

% where pCha = (H*U[0,180])/(H*U[0,180] + (1-H) P(lambda | mu_hat,sig_hat))
% dimension 1: internal representation
% dimension 2: sound location
% dimension 3: change-point variable


gridTicks = 0:180;
priorPdf = normpdf(gridTicks,priorMu,priorStd);
priorPdf = priorPdf/sum(priorPdf);

normFac = normcdf(180,priorMu,priorStd) - normcdf(0,priorMu,priorStd);


p_change = (1.0/180)./(1.0/180 + ...
                   normpdf(likeMu,priorMu,priorStd)/normFac);
               


p_Lambda_given_O = normpdf(gridTicks,likeMu,likeStd);
p_Lambda_given_O = p_Lambda_given_O/sum(p_Lambda_given_O);


%prior_no_change = priorPdf.*p_Lambda_given_O;
%prior_no_change = prior_no_change/sum(prior_no_change);


%posteriorDist = haz*p_change*p_Lambda_given_O + (1-haz)*(1-p_change)*prior_no_change;

new_prior = haz*1.0/181 + (1-haz)*priorPdf;

posteriorDist = new_prior.*p_Lambda_given_O;
posteriorDist = posteriorDist/sum(posteriorDist);


% 
% % % ---------------------------------- MRN debug --------------------------
% prior = normpdf(gridTicks, priorMu, priorStd);  %prior over the grid
% prior = prior./(sum(prior)); %normalise
% likelihood = normpdf(gridTicks, likeMu, likeStd); %likelihood over the grid
% likelihood = likelihood./sum(likelihood);
% 
% % p(out|data,s)= p(s,data|outcome)p(data,s)   (UNNORMALIZED)
% 
% 
% pS = repmat([1-haz; haz],1,length(gridTicks));
% pS = pS./sum(sum(pS));
% pDat = repmat(likelihood, 2,1);
% pDat = pDat./sum(sum(pDat));
% 
% % probability of data and S given the observed outcome
% pDatAndS=pDat.*pS;
% pDatAndS=pDatAndS./sum(sum(pDatAndS));
% 
% %prior on outcome is stored prior if s=0, flat if s=1;
% flat=ones(1, length(gridTicks))./length(gridTicks);
% pOutcome=[prior; flat]; %prior over means
% 
% jointDist=pOutcome.*pDatAndS;
% jointDist=jointDist./sum(sum(jointDist));
% posteriorDist = nansum(jointDist);
% % % sDist=nansum(jointDist, 2);
% % % 
% % % posteriorMu=nansum(outcomeDist.*gridTicks);
% % % postS=sDist(2);
% % 
% % %------------------------------------------------------------------------
% % 
% 









%
%
%
% % aux variable for Bayesian combo of N(lamda,sig_like) & N(mu_hat,sig_hat)
%
% aux_internal_pdf = nan(length(gridTicks),length(gridTicks));
% for i = 1:length(gridTicks)
%     aux_internal_pdf(i,:)=normpdf(gridTicks, gridTicks(i), likeStd);
%     aux_internal_pdf(i,:)= aux_internal_pdf(i,:)./(sum(aux_internal_pdf(i,:)));
% end
%
% p_Ohat_noChange = bsxfun(@times,priorPdf,aux_internal_pdf);
% aux_norm = sum(p_Ohat_noChange,2);
% p_Ohat_noChange = bsxfun(@times,1./aux_norm,p_Ohat_noChange);
%
% % Note each row of p_Ohat_noChange corresponds to a pdf for a value of
% % lambda
%
% p_Ohat_change = aux_internal_pdf;
%
% p_Lambda_given_O = normpdf(gridTicks,likeMu,likeStd);
% p_Lambda_given_O = p_Lambda_given_O/sum(p_Lambda_given_O);
%
%
% aux_mult_CP = repmat(pCha'.*p_Lambda_given_O',1,length(gridTicks));
% aux_mult_noCP = repmat((1-pCha').*p_Lambda_given_O',1,length(gridTicks));
%
%
% posteriorDist = p_Ohat_change.*aux_mult_CP + p_Ohat_noChange.*aux_mult_noCP;
% posteriorDist = posteriorDist./sum(posteriorDist(:));








%
% if nargout>2
%     likeFromPrior=internalRepDist.* pInternalRepGivMu.*repmat(priorNoChange, 1, length(gridTicks));
%     likeFromUniform=internalRepDist.* pInternalRepGivMu.*repmat(priorChange, 1, length(gridTicks));
%     pCha= haz.*sum(likeFromUniform(:))./( sum(likeFromPrior(:)).*(1-haz)+ sum(likeFromUniform(:)).*haz);
% end



%CPP = p_change;
CPP = 1.0;

%imagesc(posteriorDist)
%imagesc(pDat)
%imagesc(jointDist)
%plot(gridTicks,pDist)
end
