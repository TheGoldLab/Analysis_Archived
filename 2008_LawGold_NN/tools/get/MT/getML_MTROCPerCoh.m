function [b, bi, pv] = getML_MTROCPerCoh(Monk, recompute)
% MT ROC area at each coh

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_MTROCPerCoh_' Monk '.mat'];


if recompute
    [r rn re] = getML_ROCs([Monk 'TRain_MT.txt'], [0 100], [], [], [], [], 0);

    r  = shiftdim(mynanmedian(r,1),2);
    rn = shiftdim(nanmean(rn,1),2);
    
    
    % get performance improvement vs coh
    ses = 1:length(r);
    for i = 1:7
        Lgd = ~isnan(r(:,i));
        [b_, bi_, h, p] = nestedFW(r(Lgd,i), ones(sum(Lgd),1), [ones(size(r(Lgd,i))) ses(Lgd)']);
        b(i)  = b_(2);
        bi(i) = b_(2)-bi_(2,1);
        pv(i) = p;
    end


    save(savepath, 'b', 'bi', 'pv')

else
    load(savepath)
end



