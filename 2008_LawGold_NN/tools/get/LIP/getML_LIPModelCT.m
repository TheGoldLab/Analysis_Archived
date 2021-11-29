function [fits, sems, p] = getML_LIPModelCT(Monk, recompute)

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_LIPModel_' Monk '.mat'];

if recompute
    a      = getML_txt([Monk 'TRain_LIP.txt']);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    ses    = a.data{strcmp(a.name,'session')};
    usable = a.data{strcmp(a.name,'usable')};

    crtf = [1];
    nvf  = [0 1];
    dirf = 0;

    bb     = -200;
    be     = 1000;
    bs     = 25;
    bw     = 50;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    mbins  = mean(bins,2)./1000;

    [rM, rSD, rN, rMn] = getML_mRate([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, 0);


    on   = nans(length(fn),1);
    off  = nans(length(fn),1);
    fits = nans(4,length(fn));
    sems = nans(4,length(fn));
    p    = nans(4,length(fn));

    m      = rM(1:7,:,:);                           
    sd     = rSD(1:7,:,:)./(rN(1:7,:,:).^0.5);      % standard error
    
    coh = [0 0.032 0.064 0.128 0.256 0.512 0.999]';

    for i = 1:length(fn)
        if usable(i)==1
            coh = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
            fprintf(sprintf('%d: %s\n', i, fn{i}))
            warning off

            % fit raw spikes with s-ramp to estimate onset time
            if strcmp(Monk, 'Cy')
                L = mbins>=0.05 & mbins<=0.7;    % cyrus's time
            elseif strcmp(Monk, 'ZZ')
                L = mbins>=0.25 & mbins<=1.0;      % zsa zsa's time
            end

            % use either 99.9 or 51.2% to estimate onset time, depending on
            % which coh has more trials (later sessions 51.2% has more trials)
            if max(rN(7,:,i)) > max(rN(6,:,i))
                bl  = nanmean(nanmean(m(7,mbins>=0 & mbins<=0.2,i)));
                [bHi, bHise, gof, y]= sramp_fitW(mbins(L)', m(7,L,i), sd(7,L,i), bl, [],...
                    [min(mbins(L)) 0 min(mbins(L)); max(mbins(L)) 1000 max(mbins(L))]);
                on(i)  = bHi(1);
                off(i) = bHi(3); 
            else
                bl  = nanmean(nanmean(m(6,mbins>=0 & mbins<=2,i)));
                [bHi, bHise, gof, y]= sramp_fitW(mbins(L)', m(6,L,i), sd(6,L,i), bl, [],...
                    [min(mbins(L)) 0 min(mbins(L)); max(mbins(L)) 1000 max(mbins(L))]);
                on(i)  = bHi(1);
                off(i) = bHi(3);
            end

            if on(i)>off(i)
                xxx    = off(i);
                off(i) = on(i);
                on(i)  = xxx;
            end
%             if on(i)>0.9
%                 on(i)=0.9;
%             end
%              
%             % fit LIP model
%             if strcmp(Monk, 'Cy')
%                 fitLength = 275;    % cyrus's time
%             elseif strcmp(Monk, 'ZZ')
%                 fitLength = 325;    % zsa zsa's time
%             end
% 
%             
%             if strcmp(Monk,'Cy') & i==12 | strcmp(Monk,'Cy') & i==46
%                 % set onset time of these two cells because the algorithm to detect 
%                 % onset time didn't work very well for this cell
%                 init  = [];
%                 on(i) = on(i)+200;
%             else
%                 init  = [];
%             end

            init  = [];
            
            coh   = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
            ot    = on(i)  - mod(on(i),0.025);
            ft    = off(i) + mod(off(i),0.025);
            LP    = mbins>=ot & mbins<=ft;
            mbP   = mbins(LP);
            mbP   = repmat(mbP',length(coh),1);
            mbP   = mbP(:);
            c     = repmat(coh,1,sum(LP));
            c     = c(:);
            r     = m(1:7,LP,i);
            rd    = sd(1:7,LP,i);
            r     = r(:);
            rd    = rd(:);

            
            [fits(:,i), sems(:,i)] = getML_neuroACT([c(:) mbP(:) r(:) rd(:)]);
     
            fits(:,i);
            % plot and check fits
            if 0
                clf
                lc  = [];
                for j = 1:7
                    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
                end
                
                for j = 1:7
                   hold on
                   plot(mbins, m(j,:,i)', 'o', 'MarkerSize', 10, ...
                                         'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
                   x  = mbins(LP);
                   ym = fits(1,i) + fits(2,i) .* (coh(j).^fits(3,i)) .* ((x-x(1)).^fits(4,i));
                   plot(x, ym, '-', 'LineWidth', 3, 'Color', lc(j,:))
                   xlim([0 0.8]);              
                end
               % title(sprintf('%.3f - %.3f', a(3,i), p(i)))
                pause
            end
                
               
        end
    end

    save(savepath, 'fits', 'sems', 'p')
else

    load(savepath)
end
