function [a, ase, p] = getML_LIPModel(Monk, recompute)

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
    mbins  = mean(bins,2);

    [rM, rSD, rN, rMn] = getML_mRate([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, 0);


    % find onset time like before, fit a line separately for each coh, plot
    % each and see which function (linear coh, log coh or power function) fit
    % the best
    on  = nans(length(fn),1);
    a   = nans(4,length(fn));
    ase = nans(4,length(fn));
    p   = nans(length(fn),1);

    m      = rM;%-rM([8:14 1:7], :, :);    % select between normalized and unnormalized spikes
    sd     = rSD./(rN.^0.5);               % standard error
    %sd(sd==Inf) = 100;
    coh = [0 0.032 0.064 0.128 0.256 0.512 0.999]';

    for i = length(fn):-1:1
        % preferred
        if usable(i)==1
            coh = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
            fprintf(sprintf('%d: %s\n', i, fn{i}))
            warning off

            % fit raw spikes with s-ramp to estimate onset time
            if strcmp(Monk, 'Cy')
                L = mbins>=50 & mbins<=700;     % cyrus's time
            elseif strcmp(Monk, 'ZZ')
                L = mbins>=400 & mbins<=1000;   % zsa zsa's time was 250
            end

            % use either 99.9 or 51.2% to estimate onset time, depending on
            % which coh has more trials (later sessions 51.2% has more trials)
            if max(rN(7,:,i)) > max(rN(6,:,i))
                bl  = nanmean(nanmean(m(7,mbins>=0 & mbins<=200,i)));
                [bHi, bHise, gof, y]= sramp_fitW(mbins(L)', m(7,L,i), sd(7,L,i), bl, [],...
                    [min(mbins(L)) 0 min(mbins(L)); max(mbins(L)) 3 max(mbins(L))]);
                on(i) = bHi(1);
            else
                bl  = nanmean(nanmean(m(6,mbins>=0 & mbins<=200,i)));
                [bHi, bHise, gof, y]= sramp_fitW(mbins(L)', m(6,L,i), sd(6,L,i), bl, [],...
                    [min(mbins(L)) 0 min(mbins(L)); max(mbins(L)) 3 max(mbins(L))]);
                on(i) = bHi(1);
            end

%             if on(i)>0.9
%                 on(i)=0.9;
%             end
             
            % fit LIP model
            if strcmp(Monk, 'Cy')
                fitLength = 275;    % cyrus's parameters
                bcon      = [0 -0.020 -0.050 0; 100 0.2 1.5 100];
            elseif strcmp(Monk, 'ZZ')
                fitLength = 325;    % zsa zsa's time
                bcon      = [0 -0.020 -0.025 0; 100 0.2 1.5 100];
            end

            
            if strcmp(Monk,'Cy') & i==12 | strcmp(Monk,'Cy') & i==46
                % set onset time of these two cells because the algorithm to detect 
                % onset time didn't work very well for this cell
                init  = [];
                on(i) = on(i)+200;
            else
                init  = [];
            end
           
            if strcmp(Monk,'ZZ') & i==6 | strcmp(Monk,'ZZ') & i==16 | strcmp(Monk,'ZZ') & i==25 ...
                | strcmp(Monk,'ZZ') & i==86
                % set initial condition for these cells
                init  = [];
                on(i) = on(i)+200;
                
            elseif strcmp(Monk,'ZZ') & i==52
                init  = [];
                on(i) = on(i)+300;
                
            else
                init  = [];
            end
           
            %on(i) = on(i)+100;
            
            coh   = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
            ot    = on(i)-mod(on(i),25);
            LP    = mbins>=ot & mbins<=ot+fitLength;
            mbP   = mbins(LP);
            mbP   = repmat(mbP',length(coh),1);
            mbP   = mbP(:);
            c     = repmat(coh,1,sum(LP));
            c     = c(:);
            r     = m(1:7,LP,i);
            rd    = sd(1:7,LP,i);
            r     = r(:);
            rd    = rd(:);

            
            [a(:,i), ase(:,i), p(i), yp] = getML_fitLIPslope3(mbP, c, r, rd, init,...
                [0 -0.020 -0.025 0; 100 0.2 1.5 100]);
            
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
                   s  = a(2,i)+coh(j)*a(3,i);
                   x = mbins(LP);
                   ym = s.*(x-x(1))+a(1,i);
                   if any(ym>a(4,i))
                       ym(ym>a(4,i)) = a(4,i);
                   end
                   plot(x, ym, '-', 'LineWidth', 3, 'Color', lc(j,:))
                   xlim([0 1000]);              
                end
                title(sprintf('%.3f - %.3f', a(3,i), p(i)))
                pause
            end
                
               
        end
    end

    save(savepath, 'a', 'ase', 'p')
else

    load(savepath)
end
