function [a, ae, p] = getML_LIPModel(Monk, recompute, on)

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_LIPModel_' Monk '.mat'];

if recompute
    a      = getML_txt([Monk 'TRain_LIP.txt']);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    ses    = a.data{strcmp(a.name,'session')};
    usable = a.data{strcmp(a.name,'usable')};

    pfn   = [Monk 'TRain_psy.txt'];
    xx    = getML_txt(pfn);
    pses  = xx.data{strcmp(xx.name,'session')};
    pdir = xx.data{strcmp(xx.name,'ddir')};
    dd = nans(size(ses));
    for i = 1:length(ses)      % get dots direction
        dd(i) = pdir(pses==ses(i));
    end

    
    crtf = [1];
    nvf  = [0 1];
    dirf = 0;

    bb     = -200;
    be     = 1500;
    bs     = 25;
    bw     = 50;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    mbins  = mean(bins,2);

    [rM, rSD, rN, rMn] = getML_mRate([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, 0);


    % find onset time like before, fit a line separately for each coh, plot
    % each and see which function (linear coh, log coh or power function) fit
    % the best
    if nargin<3 | isempty(on)
        on  = nans(length(fn),1);
    end
    a   = nans(4,length(fn));
    ae  = nans(4,2,length(fn));
    p   = nans(length(fn),1);

    m      = rM-rM([8:14 1:7], :, :);    % select between normalized and unnormalized spikes
    sd     = rSD./(rN.^0.5);                       
    %sd(sd==Inf) = 100;
    coh  = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
    Lcoh = logical([1 1 1 1 1 1 1]); 
    
    for i = 1:length(fn)
        % preferred
        if usable(i)==1
            fprintf(sprintf('%d: %s\n', i, fn{i}))
            warning off

            if nargin<3 | isempty(on)
                % fit raw spikes with s-ramp to estimate onset time
                if strcmp(Monk, 'Cy')
                    L = mbins>=50 & mbins<=700;    % cyrus's time
                elseif strcmp(Monk, 'ZZ')
                    L = mbins>=250 & mbins<=1000;      % zsa zsa's time
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
            end
            
            
            % fit LIP model
            if strcmp(Monk, 'Cy')
                if i==126 | i==112 | i==110 | i==106 | i==104 | i==103 | i==69
                    fitLength = 200;    % cyrus's time
                else
                    fitLength = 275;    % cyrus's time
                end    
            elseif strcmp(Monk, 'ZZ')
                if i==87 | i==90 | i==91
                    fitLength = 175;    % zsa zsa's time
                else
                    fitLength = 200;    % zsa zsa's time
                end
            end
            
            
            if strcmp(Monk,'Cy') % adjust onset time for these cells
                if     i==89 | i==128 | i==122 | i==85 | i==24 ...
                        | i==138 | i==125 | i==116 ...
                         | i==113 | i==84 | i==49
                    init  = []; 
                    on(i) =  on(i)+25;   
                   
                elseif  i==64  | i==109 | i==70 | i==46 | i==1 | i==136 ...
                         | i==58 | i==146 | i==137  | i==130 | i==127 | i==59
                    init  = [];
                    on(i) =  on(i)+50;   
                    
                elseif  i==96 | i==53 | i==114 | i==103 | i==97 | i==56 | i==30
                    init = [];
                    on(i) = on(i)+75;
                    
                elseif  i==76 | i==12 | i==142 ...
                        | i==124 | i==57 | i==19 | i==16
                    init  = [];
                    on(i) =  on(i)+100; 
          
                elseif i==52 | i==14 
                    init = [];
                    on(i) = on(i)+150;
                    
                elseif i==140  | i==28 | i==123 | i==60 | i==45
                    init  = [];
                    on(i) =  on(i)+200;

                elseif i==93 | i==43
                    init  = [];
                    on(i) =  on(i)+250;
                
                elseif i==78 | i==90 | i==91 | i==107
                    init  = [];
                    on(i) =  on(i)+275;

                elseif i==81
                    init  = [];
                    on(i) =  on(i)+300;

                elseif i==104
                    init  = [];
                    on(i) =  on(i)+450;

                elseif i==41 | i==131 | i==77 
                    init = [];
                    on(i) = on(i)+500;
                
                elseif i==51
                    init = [];
                    on(i) = on(i)+550;
              
                elseif i==42 | i==35 | i==94 ...
                        | i==149 | i==82 | i==83 | i==126 | i==69
                    init  = [];
                    on(i) =  on(i)-25;
                    
                elseif i==147 | i==141 | i==110 | i==79
                    init  = [];
                    on(i) =  on(i)-50;
                    
                elseif i==22 | i==55 | i==66
                    init  = [];
                    on(i) =  on(i)-100;
                    
                elseif i==47 | i==4
                    init  = [];
                    on(i) =  on(i)-200;   
                    
                else
                    init = [];
                end
            end
           
             
            if strcmp(Monk,'ZZ') 
                if  i==6
                    init = [];
                    on(i) = on(i)-500;
                
                elseif i==89 | i==35
                    init = [];
                    on(i) = on(i)-400;
                
                elseif i==12 | i==62 | i==59 | i==51
                    init = [];
                    on(i) = on(i)-300;
                
                elseif i==95
                    init = [];
                    on(i) = on(i)-250;
                
                elseif i==32 | i==30 | i==74 | i==24 | i==43
                    init = [];
                    on(i) = on(i)-200;
                
                elseif i==71 | i==25
                    init = [];
                    on(i) = on(i)-150;
                    
                elseif i==47 | i==26 | i==11
                    init = []; 
                    on(i) = on(i)-100;
                    
                elseif i==70 | i==58 | i==52 | i==36 | i==22 | i==13 | i==90 
                    init = [];
                    on(i) = on(i)-50;
                
                elseif i==69 
                    init = [];
                    on(i) = on(i)-25;
                
                elseif i==49
                    init = [];
                    on(i) = on(i)+650;
                
                elseif i==1 | i==53
                    init = [];
                    on(i) = on(i)+600;
                
                elseif i==48 
                    init = [];
                    on(i) = on(i)+400;

                elseif i==64 | i==56 | i==42
                    init = [];
                    on(i) = on(i)+300;
               
                elseif i==45 | i==20 | i==57 | i==10 | i==82
                    init = [];
                    on(i) = on(i)+200;
                    
                elseif i==9
                    init = [];
                    on(i) = on(i)+325;
                    
                elseif i==77 | i==86 | i==37
                    init = [];
                    on(i) = on(i)+150;
               
                elseif i==86
                    init = [];
                    on(i) = on(i)+175;
                
                elseif i==16 | i==67 | i==39 | i==91 | i==92
                    init  = [];
                    on(i) = on(i)+100;
                     
                elseif i==79 | i==88 ...
                        | i==72 | i==66 | i==8 | i==87 | i==93
                    init  = [];
                    on(i) = on(i)+50;
                
                elseif i==76  | i==81 ...
                        | i==50 | i==85 
                    init  = [];
                    on(i) = on(i)+25;
                
                else
                    init = [];
                end
            end
           
            ot    = on(i)-mod(on(i),25);
            LP    = mbins>=ot & mbins<=ot+fitLength;
            mbP   = mbins(LP);
            mbP   = repmat(mbP',sum(Lcoh),1);
            mbP   = mbP(:);
            c     = repmat(coh(Lcoh),1,sum(LP));
            c     = c(:);
            r     = m(Lcoh,LP,i);
            rd    = sd(Lcoh,LP,i);
            r     = r(:);
            rd    = rd(:);

            if strcmp('ZZ', Monk)
                acon = [0 -1 -0.050 0; 100 0.2 1.5 100];
            else
                acon = [0 -0.020 -0.050 0; 100 0.2 1.5 100];
            end
            
            citype = {100,68};
            [a(:,i), ae(:,:,i), p(i), yp] = getML_fitLIPslope3(mbP, c, r, rd, init,...
                acon, citype);
            
            % plot and check fit
            if 0
                clf
                lc  = [];
                for j = 1:7
                    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
                end
                  
                for j = find(Lcoh)
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
                   xlim([0 1500]);              
                end
                title(sprintf('%.f (%.f,%.f)', 1000*a(3,i), 1000*ae(3,1,i), 1000*ae(3,2,i)))
                fprintf(sprintf('fits: %.f (%.f,%.f)\n', 1000*a(3,i), 1000*ae(3,1,i), 1000*ae(3,2,i)))
                pause
            end
                
               
        end
    end

    a  = real(a);
    ae = real(ae); 
    save(savepath, 'a', 'ae', 'p')
else

    load(savepath)
end
