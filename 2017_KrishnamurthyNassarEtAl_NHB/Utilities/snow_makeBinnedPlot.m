function [a B p R points]=snow_makeBinnedPlot(x,y, binSize, XLABEL, YLABEL, markSize, markColor, inCurrFig, lineOn, lc)
% bin size is in fractional units... ie binSize=.01 will make a hundred bins, each
% of which includes 1 percent of the date


 


if nargin<10|isempty(lc)
    lc=markColor
end



if nargin<9|isempty(lineOn)
    lineOn=1
end



if nargin<6|isempty(markSize)
    markSize=3
end

if nargin<7|isempty(markColor)
    markColor='k'
end

if nargin<8|isempty(inCurrFig)|inCurrFig==0
a=figure;
else
a=gcf;
end

bins=ceil(1/binSize);


siz=length(x);
binEl=binSize.*siz;
[X,I] = sort(x);
Y=(y(I));


for i = 1:bins
    sel= floor((i-1).*binEl+1):ceil(i.*binEl);
    if max(sel)>length(X)
        sel=sel(1:end-1);
    end
    
    
    meanX(i)=nanmean(X(sel));
    devX(i)=nanstd(X(sel))./sqrt(length(sel));
    meanY(i)=nanmean(Y(sel));
    devY(i)=nanstd(Y(sel))./sqrt(length(sel));
end

ErrLinesX=repmat(meanX, 2, 1)
ErrLinesY=cat(1, meanY+devY, meanY-devY);

bErrLinesX=cat(1, meanX+devX, meanX-devX); 
bErrLinesY=repmat(meanY, 2, 1)






if size(y, 2)>size(y, 1)
    y=y'
end


if size(x, 2)>size(x, 1)
    x=x'
end
%keyboard





[B,BINT,R,RINT,STATS] = regress(y, [ones(length(x), 1) x]);





if lineOn==1
plot([min(meanX) max(meanX)], [B(1)+min(meanX).*B(2) B(1)+max(meanX).*B(2)], 'color', lc)
end


R = corr(y, x)
p=STATS(3);





hold on
plot(ErrLinesX, ErrLinesY, 'k', 'lineWidth', 1)
plot(bErrLinesX, bErrLinesY, 'k', 'lineWidth', 1)
points=plot(meanX, meanY, 'o', 'markersize', markSize, 'markerEdgeColor', 'k', 'markerFaceColor', markColor, 'lineWidth', 1)
xlabel(XLABEL)
ylabel(YLABEL)










   




