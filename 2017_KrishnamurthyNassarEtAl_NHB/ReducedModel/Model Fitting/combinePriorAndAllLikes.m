function [posteriorDist, gridTicks, pCha]=combinePriorAndAllLikes(priorMu, priorStd, likeMu, likeStd, haz)

%% This code is based on combinePriorAndLike but extended for use in model fitting.
% In particular, since we do not know the particular likelihood
% representation, we'd like to consider all possible internal
% representations.


% dimension 1: internal representation
% dimension 2: sound location
% dimension 3: change-point variable

gridTicks = 0:180;

priorNoChange = normpdf(gridTicks, priorMu, priorStd)';  %prior over the grid
priorNoChange=priorNoChange./sum(priorNoChange);



priorChange=  ones(length(gridTicks), 1)./length(gridTicks);
prior=priorNoChange.*(1-haz) + priorChange.*(haz);
prior=repmat(prior, 1, length(gridTicks));
prior = prior./(sum(prior(:))); %normalise



% create matrix of all possible internal representations (replicated across
% means)
internalRepDist=normpdf(gridTicks, likeMu, likeStd);
internalRepDist=repmat(internalRepDist, length(gridTicks), 1);
internalRepDist=internalRepDist./sum(internalRepDist(:));

% get conditional distribution of internal representations given means:
for i = 1:length(gridTicks)
    pInternalRepGivMu(i,:)=normpdf(gridTicks, gridTicks(i), likeStd);
    pInternalRepGivMu(i,:)= pInternalRepGivMu(i,:)./(sum(pInternalRepGivMu(i,:)));
end

posteriorDist=  internalRepDist.* pInternalRepGivMu.*prior;
posteriorDist=posteriorDist./sum(posteriorDist(:));

if nargout>2
    likeFromPrior=internalRepDist.* pInternalRepGivMu.*repmat(priorNoChange, 1, length(gridTicks));
    likeFromUniform=internalRepDist.* pInternalRepGivMu.*repmat(priorChange, 1, length(gridTicks));
    pCha= haz.*sum(likeFromUniform(:))./( sum(likeFromPrior(:)).*(1-haz)+ sum(likeFromUniform(:)).*haz);
end

    
    
 


%imagesc(posteriorDist)
%imagesc(pDat)
%imagesc(jointDist)
%plot(gridTicks,pDist)
end
