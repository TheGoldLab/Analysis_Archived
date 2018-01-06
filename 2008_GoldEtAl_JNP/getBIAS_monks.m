function [monks_, monkn_, monkse_] = getBIAS_monks(exp)
% function [monks_, monkn_, monkse_] = getBIAS_monks(exp)

inds = 1:4;
if nargin >= 1
    if strncmpi(exp, 'FEF', 2)
        inds = 1:2;
    elseif strncmpi(exp, 'MTLIP', 2)
        inds = 3:4;
    end
end

monks_ = {'Atticus', 'Ava', 'Cyrus', 'ZZ'};
monks_ = monks_(inds);

if nargout > 1
    monkn_ = length(monks_);
end

if nargout > 2
    se      = [161 222 160 130];
    monkse_ = se(inds);
end
