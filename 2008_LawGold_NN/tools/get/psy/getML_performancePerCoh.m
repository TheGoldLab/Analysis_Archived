function [p n pci] = getML_performancePerCoh(fn, tlim, recompute)
% get average performance at each coherence for each session

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_performancePerCoh_' fn '.mat'];

if nargin<2
    tlim      = [-inf inf];
    recompute = 1;
elseif nargin<3
    recompute = 1;
end


if recompute
    % load session infos
    a     = getML_txt(fn);
    fn    = a.data{strcmp(a.name,'dat_fn')};
    ses   = a.data{strcmp(a.name,'session')};


    % get performance
    coh    = [0 3.2 6.4 12.8 25.6 51.2 99.9];
    p      = nans(length(coh),length(fn));
    pci    = nans(length(coh),2,length(fn));
    n      = nans(length(coh),length(fn));
   
    global FIRA
    for i = 1:length(fn)
        openFIRA(fn{i})
        fprintf('%d: %s\n', i, fn{i})
        co = getFIRA_ecodesByName('dot_coh');
        cr = getFIRA_ecodesByName('correct');
        tt = getFIRA_ecodesByName('dot_dur');

        % remove nan trials
        Lbd     = isnan(co)|isnan(cr);
        co(Lbd) = [];
        cr(Lbd) = [];
        tt(Lbd) = [];
        
        % remove broken fix and no choice trials
        Lbd     = cr<0;
        co(Lbd) = [];
        cr(Lbd) = [];
        tt(Lbd) = [];
        
        % select for time
        Lbd     = tt<=tlim(1) & tt>=tlim(2);
        if ~isempty(Lbd)
            co(Lbd) = [];
            cr(Lbd) = [];
            tt(Lbd) = [];
        end

        % get performance per coh
        for j = 1:length(coh)
            Lcoh = co==coh(j);
            if nansum(Lcoh&cr>=0)>=10   % only calculate percent correct if trials>=10
                n(j,i)           = nansum(Lcoh&cr>=0);
                p(j,i)           = nansum(Lcoh&cr==1)/nansum(Lcoh&cr>=0);
                [xxx pci(j,:,i)] = binofit(p(j,i)*n(j,i),n(j,i),1-0.68); % 68% ci
            end
        end
    end

    save(savepath, 'p', 'n', 'pci')
else

    load(savepath)
end

