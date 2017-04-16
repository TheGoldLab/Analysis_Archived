function output=clusterBasedPermutationTesting(input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input  -- structure containing the following fields:

% data   : M X N matrix containing coefficients for M subjects at N timepoints
% clt    : cluster forming threshold p value... (value you need to beat to
% be considered as part of a cluster.
% numReps: number of permutations
% corrType: 't' = equivalent to voxelwise correction, 'mass' = cluster
% mass, 'size' cluster size. 



% 

% output -- structure containing output

% clusterMass: cluster mass for each timepoint
% clusterSize: cluster size for each timepoint
%       absT : absolute T value for each timepoint
% permutationMass= numReps array of max permutation cluster mass (size*mean t statistic)
% permutationSize= numReps array of max permutation cluster size
% permutationT   = numReps array of max abs(T) value. 
% corrP          = corrected p value (based on correction typ




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=size(input.data, 1);
[~, p, ~, STATS]=ttest(input.data);
[output.clusterSize, output.clusterMass]=getNumNeigbors(p<input.clt, abs(STATS.tstat)) ;
output.absT=abs(STATS.tstat);

% create permutation distributions for potential test statistics:
output.permutationMass=nan(input.numReps, 1);
output.permutationSize=nan(input.numReps, 1);
output.permutationT=nan(input.numReps, 1);

% get null distribution based on permuted maxima:

for j = 1:input.numReps
    
    permMat=repmat((binornd(1, .5, n, 1)-.5).*2, 1, length(p));
    [~,sP,~,sSTATS] = ttest(permMat.*input.data) ;
    simSig=sP<input.clt;
    output.permutationT(j)=max(abs(sSTATS.tstat));
    % permutationMass= numReps array of max permutation cluster mass (size*mean t statistic)
    % permutationSize= numReps array of max permutation cluster size
    % permutationT   = numReps array of max abs(T) value.
    
    if sum(simSig)==0
        output.permutationMass(j)=0;
        output.permutationSize(j)=0;
    else
        [runLength, runMass]=getNumNeigbors(simSig, abs(sSTATS.tstat));
        output.permutationMass(j)=max(runMass);
        output.permutationSize(j)=max(runLength);
    end
end

% compute corrected p values for each time bin:
output.corrP=nan(size(output.absT));
switch input.corrType
    case  't'
        for i=1:length(output.absT)
            output.corrP(i)=nanmean(output.absT(i) <= output.permutationT);
        end
    case 'mass'
        for i=1:length(output.absT)
            output.corrP(i)=nanmean(output.clusterMass(i) <= output.permutationMass);
        end
    case 'size'
        for i=1:length(output.absT)
            output.corrP(i)=nanmean(output.clusterSize(i) <= output.permutationSize);
        end
end

