function [rp, rpci, rpp, rpm, rth, rthci, rthp, rthm] = getML_learningRate(monk, wlen, recompute)
% get learning rate of behavioral thresholds and lapse


savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_learningRate_' monk '_' int2str(wlen) '.mat'];

if recompute
    % get data
    [p n pci]        = getML_performanceAsso([monk 'TRain_asso.txt'], 0);
    [pt nt ptci]     = getML_performanceAsso([monk 'TRain_psy.txt'], 0);
    [fits, sems, th] = getML_psyPerformanceCTFixLapse([monk 'TRain_psy.txt'], 0, 1-pt);

    
    p    = 1-[p pt];
    pci  = [pci ptci];
    pd   = (pci(2,:)-pci(1,:))/2;
    th   = 100*fits(1,:);
    thci = shiftdim(sems(1,:,:),1)*100;
    thd  = (thci(2,:)-thci(1,:))/2;

   
    % get learning rate for lapse
    rp    = nans(length(p),1);
    rpci  = nans(length(p),1);
    rpp   = nans(length(p),1);
    rpm   = nans(length(p),1);

    rth   = nans(length(th),1);
    rthci = nans(length(th),1);
    rthp  = nans(length(th),1);
    rthm  = nans(length(th),1);

    dlen  = wlen/2; % half window size
    [b bi gof ym] = exp_fit2W(1:length(p), p, pd, [0.5; 5], [0.4 0.6;0.01 100]);
    for i = 1:length(rp)
        if i-dlen-1>=1 & i+dlen-1<=length(rp)
            % data
            y               = p(i-dlen-1:i+dlen-1)';
            yd              = pd(i-dlen-1:i+dlen-1)';
            [b,bi,h,rpp(i)] = nestedFW(y, yd, [ones(size(y)) [1:2*dlen+1]'], 1-0.68);
            rp(i)           = b(2);
            rpci(i)         = bi(2,2)-b(2);

            % model
            y               = ym(i-dlen-1:i+dlen-1);
            yd              = ones(size(y));
            b               = nestedFW(y, yd, [ones(size(y)) [1:2*dlen+1]']);
            rpm(i)          = b(2);
        elseif i-dlen-1<1
            % data
            y               = p(1:i+dlen-1)';
            yd              = pd(1:i+dlen-1)';
            [b,bi,h,rpp(i)] = nestedFW(y, yd, [ones(size(y)) [1:i+dlen-1]'], 1-0.68);
            rp(i)           = b(2);
            rpci(i)         = bi(2,2)-b(2);

            % model
            y               = ym(1:i+dlen-1);
            yd              = ones(size(y));
            b               = nestedFW(y, yd, [ones(size(y)) [1:i+dlen-1]']);
            rpm(i)          = b(2);
        elseif i+dlen-1>length(rp)
            % data
            y               = p(i-dlen-1:length(rp))';
            yd              = pd(i-dlen-1:length(rp))';
            [b,bi,h,rpp(i)] = nestedFW(y, yd, [ones(size(y)) [1:length(rp)-i+dlen+2]'], 1-0.68);
            rp(i)           = b(2);
            rpci(i)         = bi(2,2)-b(2);

            % model
            y               = ym(i-dlen-1:length(rp));
            yd              = ones(size(y));
            b               = nestedFW(y, yd, [ones(size(y)) [1:length(rp)-i+dlen+2]']);
            rpm(i)          = b(2);
        end
    end

  
    
    
    
    % get learning rate for log thresholds
    if strcmp(monk, 'Cy')
        [b bci gof ym] = exp_fitW(1:length(th), th, ones(size(thd)), [10; 80; 40], [0 20; 0 100; 0 160]);
    else
        [b bci gof ym] = exp_fitW(1:length(th), th, ones(size(thd)), [10; 80; 40], [0 20; 0 100; 0 160]);
    end
    
    
    thd(thd<1) = nanmean(thd(thd<1)); % otherwise, log(thd) will be negative and will be troublesome for weighted fit
    for i = 1:length(rth)
        if i-dlen-1>=1 & i+dlen-1<=length(rth)
            % data
            y               = log(th(i-dlen-1:i+dlen-1)');
            yd              = log(thd(i-dlen-1:i+dlen-1)');
            [b,bi,h,rthp(i)] = nestedFW(y, yd, [ones(size(y)) [1:2*dlen+1]'], 1-0.68);
            rth(i)           = b(2);
            rthci(i)         = bi(2,2)-b(2);

            % model
            y               = log(ym(i-dlen-1:i+dlen-1));
            yd              = ones(size(y));
            b               = nestedFW(y, yd, [ones(size(y)) [1:2*dlen+1]']);
            rthm(i)          = b(2);
        elseif i-dlen-1<1
            % data
            y               = log(th(1:i+dlen-1)');
            yd              = log(thd(1:i+dlen-1)');
            [b,bi,h,rthp(i)] = nestedFW(y, yd, [ones(size(y)) [1:i+dlen-1]'], 1-0.68);
            rth(i)           = b(2);
            rthci(i)         = bi(2,2)-b(2);

            % model
            y               = log(ym(1:i+dlen-1));
            yd              = ones(size(y));
            b               = nestedFW(y, yd, [ones(size(y)) [1:i+dlen-1]']);
            rthm(i)          = b(2);
        elseif i+dlen-1>length(rth)
            % data
            y               = log(th(i-dlen-1:length(rth))');
            yd              = log(thd(i-dlen-1:length(rth))');
            [b,bi,h,rthp(i)] = nestedFW(y, yd, [ones(size(y)) [1:length(rth)-i+dlen+2]'], 1-0.68);
            rth(i)           = b(2);
            rthci(i)         = bi(2,2)-b(2);

            % model
            y               = log(ym(i-dlen-1:length(rth)));
            yd              = ones(size(y));
            b               = nestedFW(y, yd, [ones(size(y)) [1:length(rth)-i+dlen+2]']);
            rthm(i)          = b(2);
        end
    end
  
    
    
    save(savepath, 'rp', 'rpci', 'rpp', 'rpm', 'rth', 'rthci', 'rthp', 'rthm')
else
    load(savepath)
end

