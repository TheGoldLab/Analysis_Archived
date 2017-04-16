function [numNeighbors neighborhoodMass]=getNumNeigbors(logicalArray, mass)



numNeighbors=nan(size(logicalArray));
neighborhoodMass=nan(size(logicalArray));
ct=0;
for i = 1:length(logicalArray)
    if ~logicalArray(i)
        numNeighbors(i)=0;
        neighborhoodMass(i)=0;
        ct=0;
        
    else
        ct=ct+1;
        numNeighbors(i-ct+1:i)=ct;
        neighborhoodMass(i-ct+1:i)=nansum(mass(i-ct+1:i));
    end
end
