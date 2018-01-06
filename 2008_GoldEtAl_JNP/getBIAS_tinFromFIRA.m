function [tin_, tq_] = getBIAS_tinFromFIRA(sid, task)
% function [tin_, tq_] = getBIAS_tinFromFIRA(sid, task)

if nargin < 2 || isempty(task)
    task = 6;
end

global FIRA

tin_ = nan;
tq_  = nan;
rmax = nan;

if isempty(FIRA)
    return
end

Ltst = getFIRA_ecodesByName('task')==task;
if task == 6
    Ltst = Ltst & getFIRA_ecodesByName('correct')==1;
end
tdat = getFIRA_ecodesByName({'dot_dir', 'fp_off', 'dot_on', 'dot_off'});
dirs = nonanunique(tdat(:,1));
ns   = nans(length(dirs),1);
as   = cell(length(dirs),1);

for dd = 1:length(dirs)
    Ld = Ltst&tdat(:,1)==dirs(dd);
    ns(dd) = sum(Ld);
    if sum(Ld) > 15
        % r  = getFIRA_rateByBin(Ld, sid, tdat(Ld,2)+50, tdat(Ld,2)+350);
        [r,z,as{dd}] = getFIRA_rateByBin(Ld, sid, tdat(Ld,3)+100, tdat(Ld,4));
        %disp([dirs(dd) r(1) ns(dd)])
        if ~isempty(r) && (isnan(rmax) || r(1) > rmax)
            rmax = r(1);
            tin_ = dirs(dd);
            tii  = dd;
        end
    end
end

% tuning quality ...err, uhh...
if isfinite(rmax)
    [Y,I] = sort(ns);
    toi   = I(find(I~=tii,1,'last'));
    if ns(tii) > 15
        tq_   = rocN(as{tii}, as{toi});
    end
end

if 0 %length(dirs) == 2
    cla reset; hold on;
    Fd = find(Ltst&tdat(:,1)==dirs(1));
    for yy = 1:length(Fd)
        xs = FIRA.spikes.data{Fd(yy),sid}-tdat(Fd(yy),3);
        plot(xs, yy*ones(size(xs)),'kx');
    end
    Fd = find(Ltst&tdat(:,1)==dirs(2));
    yst = yy;
    for yy = 1:length(Fd)
        xs = FIRA.spikes.data{Fd(yy),sid}-tdat(Fd(yy),3);
        plot(xs, yst+yy*ones(size(xs)),'rx');
    end
    disp([tin_==dirs(1) tq_])
    r = input('next')
end
