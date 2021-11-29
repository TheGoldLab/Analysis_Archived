function [p n pci] = getML_performanceAsso(fn, recompute)
% get average performance at high coh for each session
[home_dir, lab_dir, current_dir, tmat_dir] = dirnames;
savepath = [tmat_dir '/getML_performanceAsso_' fn '.mat'];

if recompute
    % load session infos
    a     = getML_txt(fn);
    fn    = a.data{strcmp(a.name,'dat_fn')};
    ses   = a.data{strcmp(a.name,'session')};


    % get performance
    p      = nans(1,length(fn));
    pci    = nans(2,length(fn));
    n      = nans(1,length(fn));
    pfit   = nans(1,length(fn)); % performance at 1s obtained by fitting an exponential to behavioral data
    bfit   = nans(2,length(fn)); % performance at 1s obtained by fitting an exponential to behavioral data
    pfitci = nans(2,length(fn)); % performance at 1s obtained by fitting an exponential to behavioral data

    global FIRA
    for i = 1:length(fn)
        openFIRA(fn{i})
        fprintf('%d: %s\n', i, fn{i})
        co = getFIRA_ecodesByName('dot_coh');
        cr = getFIRA_ecodesByName('correct');
        tt = getFIRA_ecodesByName('dot_dur');

        Lbd     = isnan(co)|isnan(cr);
        co(Lbd) = [];
        cr(Lbd) = [];
        tt(Lbd) = [];

        Lbd     = cr<0;
        co(Lbd) = [];
        cr(Lbd) = [];
        tt(Lbd) = [];


        % compute performance for only long viewing time trials (trials
        % longer than median vt)
        Lbd     = tt<myprctile(tt,20);
        %Lbd      = [];
        if ~isempty(Lbd)
            co(Lbd) = [];
            cr(Lbd) = [];
            tt(Lbd) = [];
        end
        nanmin(tt)
        Lcoh = co==99.9;
        n(i)      = nansum(Lcoh&cr>=0);
        [p(i) ci] = binofit(nansum(Lcoh&cr==1), n(i), 1-0.68);
        pci(:,i)  = ci;
        %             p(i)     = nansum(Lcoh&cr==1)/nansum(Lcoh&cr>=0);
        %             pci(:,i) = [p(i)-sqrt(p(i)*(1-p(i))/n(i)) ...  % 68% ci using normal approximation of binomial distribution
        %                          p(i)+sqrt(p(i)*(1-p(i))/n(i))];

        
        
        % if 99% coh was not used in that session, use 51% instead
        Lcoh = co==51.2;
        if nansum(Lcoh&cr==1)/nansum(Lcoh&cr>=0) >= p(i)  & ~(strcmp(fn,'ZZTRain_asso.txt') & i==9)
            n(i)      = nansum(Lcoh&cr>=0);
            [p(i) ci] = binofit(nansum(Lcoh&cr==1), n(i), 1-0.68);
            pci(:,i)  = ci;
        end

        %             xxx             = unique(round(tt));
        %             xxx(isnan(xxx)) = [];
        %             if length(xxx)>10       % if time is a variable, fit expontial function to get performance
        %                 [b, bsem, gof, ym] = exp_fit3W(tt(Lcoh), cr(Lcoh), ones(size(cr(Lcoh))), [sum(cr(Lcoh)==1)/length(cr(Lcoh)); myprctile(tt(Lcoh),36)], [0 1; 0 1000]);
        %                 pfit(j,i)          = b(1)+(1-b(1))*(1-exp(-1000/b(2)));
        %                 bfit(j,:,i)        = b;
        %
        %                 % estimate ci with monte carlo stimulation
        %                 mcn    = 2;
        %                 mcfits = zeros(1,mcn);
        %                 for ii = 1:mcn
        %                     if ~mod(ii,10)
        %                         disp(sprintf('Bootstrap CIs, set %d', ii))
        %                     end
        %                     % compute fit on simulated data set
        %                     bmc = exp_fit3W(tt(Lcoh), binornd(ones(size(cr(Lcoh))), ym), ones(size(cr(Lcoh))), ...
        %                         [sum(cr(Lcoh)==1)/length(cr(Lcoh)); myprctile(tt(Lcoh),36)], [0 1; 0 1000]);
        %                     mcfits(ii) = bmc(1)+(1-bmc(1))*(1-exp(-1000/bmc(2)));
        %                 end
        %                 CI            = 90;
        %                 pfitci(j,:,i) = [myprctile(mcfits,50-CI/2)' myprctile(mcfits,50+CI/2)'];
        %
        %             else                    % else, just calculate the average
        %                 pfit(j,i)     = nansum(Lcoh&cr==1)/nansum(Lcoh&cr>=0);
        %                 pfitci(j,:,i) = [p(j,i)-1.96*sqrt(p(j,i)*(1-p(j,i))/n(j,i)) ...  % 68% ci using normal approximation of binomial distribution
        %                                     p(j,i)+1.96*sqrt(p(j,i)*(1-p(j,i))/n(j,i))];
        %             end
    end

    save(savepath, 'p', 'n', 'pci')
else

    load(savepath)
end

