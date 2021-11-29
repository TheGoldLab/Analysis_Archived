function [fits, sems, p] = getML_MTModelCT(Monk, recompute)

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_MTModelCT_' Monk '.mat'];

if recompute
    a      = getML_txt([Monk 'TRain_MT.txt']);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    ses    = a.data{strcmp(a.name,'session')};
    usable = a.data{strcmp(a.name,'usable')};

    crtf = [1];
    nvf  = [0 1];
    dirf = 0;

    bb     = 0;
    be     = 1000;
    bs     = 25;
    bw     = 50;
    bins   = [zeros(size([bb+bw:bs:be]')) [bb+bw:bs:be]'];    
    mbins  = bins(:,2)./1000;
    coh    = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
    
    [rM, rSD, rN, rMn] = getML_mRate([Monk 'TRain_MT.txt'], bins, crtf, nvf, dirf, 0);
  
    m  = rM([1:7],:,:)-rM([8:14], :, :);   % pref-null responses
    se = rSD./sqrt(rN);
    se = se(1:7,:,:)+se(8:14,:,:);
    
    t  = repmat(mbins',7,1);
    c  = repmat(coh,1,length(mbins));
    
    
    % outputs
    fits = nans(3,length(fn));
    sems = nans(3,length(fn));
    p    = nans(3,length(fn));
    
    
    for i = 1:length(fn)
        % preferred
        if usable(i)==1
            fprintf(sprintf('%d: %s\n', i, fn{i}))
            warning off
            
            m_  = m(:,:,i);
            se_ = se(:,:,i); 
            [fits(:,i), sems(:,i), p(:,i)] = getML_neuroACT([c(:) t(:) m_(:) se_(:)]);
     
            
            % plot and check fit
            if 0
                clf
                lc  = [];
                for j = 1:7
                    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
                end
                
                for j = 1:7
                   hold on
                   plot(mbins, m(j,:,i), 'o', 'MarkerSize', 10, ...
                                         'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
                   plot(mbins, fits(1)*coh(j)^fits(2)*mbins.^fits(3), '-', 'LineWidth', 3, 'Color', lc(j,:))
                   xlim([0 1]);              
                end
                title(sprintf('%.3f, %.3f, %.3f', fits(1), fits(2), fits(3)))
                pause
            end
                
               
        end
    end

    save(savepath, 'fits', 'sems', 'p')
else

    load(savepath)
end
