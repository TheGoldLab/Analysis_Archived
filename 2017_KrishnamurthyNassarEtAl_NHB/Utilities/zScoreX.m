function xMat=zScoreX(xMat)

if nanstd(xMat(:,1))==0
    xMat=xMat(:,2:end);
    
    
    xMat=(xMat-repmat(nanmean(xMat, 1), length(xMat), 1))./...
        repmat(nanstd(xMat), length(xMat), 1);
    
    xMat=[ones(length(xMat),1) xMat];
else
    
    xMat=(xMat-repmat(nanmean(xMat, 1), length(xMat), 1))./...
        repmat(nanstd(xMat), length(xMat), 1);
    
    
end
