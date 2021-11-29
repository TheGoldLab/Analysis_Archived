function [correctmatmp,Nmatmp,correctmp] = calculate_correct_mat(R,correct,Rbins,Rmin,method)

cpanl = find(R(1:end-1)>=Rmin & R(2:end)==1)+1;

correctmatmp = zeros(numel(cpanl),numel(Rbins)-1);
Nmatmp = correctmatmp;

N = numel(correct);

for cpind = 1:numel(cpanl)
    nextcp = find(R(cpanl(cpind)+1:end)==1,1) + cpanl(cpind);
    if isempty(nextcp)
        nextcp = N + 1;
    end
    
    correctmp = correct(cpanl(cpind):nextcp-1);
    Rtmp = R(cpanl(cpind):nextcp-1);
    
    for rind = 1:numel(Rbins)-1
        inds = Rtmp>=Rbins(rind) & Rtmp<Rbins(rind+1);
        correctmatmp(cpind,rind) = sum(correctmp(inds));
        Nmatmp(cpind,rind) = sum(inds);
    end
    
end

if nargin<5
    method = 1;
end

switch method
    case 1
        correctmp = sum(correctmatmp) ./ sum(Nmatmp);
    case 2
        correctmp = sum(correctmatmp) ./ sum(Nmatmp);
end
