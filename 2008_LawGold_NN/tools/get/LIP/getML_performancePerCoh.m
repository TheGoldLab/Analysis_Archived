function [pc, pce, b, bi, pv] = getML_performancePerCoh(Monk, recompute)
% performance at each coh

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_performancePerCoh_' Monk '.mat'];


if recompute
    fn    = [Monk 'TRain_psy.txt'];
    a     = getML_txt(fn);
    fname = a.data{strcmp(a.name,'dat_fn')};
    pc     = nans(length(fname),7);
    pce    = nans(length(fname),7);


    % get fits
    coh = [0 3.2 6.4 12.8 25.6 51.2 99.9];
    global FIRA
    for i = 1:length(fname)
        fprintf('%d: %s\n', i, fname{i})
        openFIRA(fname{i})
        co  = getFIRA_ecodesByName('dot_coh');
        crt = getFIRA_ecodesByName('correct');

        for j = 1:7
            Lgd = co==coh(j);
            c   = crt(Lgd);
            [pc(i,j) pce_]  = binofit(sum(c==1),sum(c>=0));
            pce(i,j)        = pce_(2)-pc(i,j);
        end
    end


    % get performance improvement vs coh
    ses = 1:length(fname);
    for i = 1:7
        Lgd = ~isnan(pc(:,i)) & ses'<100;
        [b_, bi_, h, p] = nestedFW(pc(Lgd,i), ones(sum(Lgd),1), [ones(size(pc(Lgd,i))) ses(Lgd)']);
        b(i)  = b_(2);
        bi(i) = b_(2)-bi_(2,1);
        pv(i) = p;
    end


    save(savepath, 'pc', 'pce', 'b', 'bi', 'pv')

else
    load(savepath)
end



