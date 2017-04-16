function [a points] = makeTACfig(yDat, xDat, cutoff, inCurFig, MC, ms, offset,LC)

if nargin<4|inCurFig==0
a=figure;
else
a=gcf;
end

if nargin<5|isempty(LC)
    LC='b'
end

if nargin<6|isempty(ms)
    ms=6;
end


if nargin<7|isempty(offset)
    offset=0;
end




hold on


xDat(xDat>cutoff)=cutoff;

for i = 1:cutoff
sel=xDat==i;
anyT(i)=nanmean(yDat(sel));
anyTs(i)=nanstd(yDat(sel))./sqrt(sum(isfinite(yDat(sel))));
end

xes=[1:cutoff; 1:cutoff];
yes=[anyT+anyTs; anyT-anyTs];

points=plot((xes)+offset, yes, 'color', LC);
plot((xes(1,:))+offset, anyT, 'o', 'color', MC, 'markerFaceColor', MC, 'markerEdgeColor', LC, 'markerSize', ms);


