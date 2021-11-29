function [p n pci ratio] = getML_performancePerCohAdjusted(fn, tlim, recompute)
% get average performance at each coherence for each session for lapse

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_performancePerCohAdjusted_' fn '.mat'];

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
    coh    = [0 3.2 6.4 12.8 25.6 51.2 99.9]';
    p      = nans(length(coh),length(fn));
    pci    = nans(length(coh),2,length(fn));
    n      = nans(length(coh),length(fn));
    ratio  = nans(length(coh),length(fn));
    
    
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
                nc(j,1) = nansum(Lcoh&cr>=0);
                pc(j,1) = nansum(Lcoh&cr==1)/nansum(Lcoh&cr>=0);
            end
        end
        % fit weibull to performance
        f = ctPsych_fit('quick3',coh/100,[pc nc]);
        
        % compute the ratio between fitted performance with and without
        % forcing lapse to zero, and use this ratio to adjust the raw
        % performance
        ratio(:,i) = quick3([f(1:end-1);0],coh/100)./quick3(f,coh/100);
        p(:,i)     = pc.*ratio(:,i);
        n(:,i)     = nc;
    end

    save(savepath, 'p', 'n', 'pci', 'ratio')
else

    load(savepath)
end

