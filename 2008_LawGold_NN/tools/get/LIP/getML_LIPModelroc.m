function [a, ae, l] = getML_LIPModelroc(Monk, recompute)
% fit LIP roc areas to a piecewise-linear slope model

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_LIPModelroc_' Monk '.mat'];

if recompute
    a      = getML_txt([Monk 'TRain_LIP.txt']);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    usable = a.data{strcmp(a.name,'usable')};

    % get ROC areas
    crtf = [1];
    nvf  = [0 1];
    dirf = 0;

    bb     = -200;
    be     = 1500;
    bs     = 50;
    bw     = 100;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    mbins  = mean(bins,2)/1000;

    [rs rn re] = getML_ROCs([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, [], 0);

    a  = nans(4,length(fn));
    ae = nans(4,2,length(fn));
    l  = nans(1,length(fn));

    m    = rs;   % select between normalized and unnormalized spikes
    sd   = re./(rn.^0.5);   % standard error
    coh  = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
    if strcmp(Monk,'Cy')
        Lcoh = logical([0 1 1 1 1 1 1]);
    else
        Lcoh = logical([0 1 1 1 1 1 1]);
    end
    
    for i = 1:length(fn)
        % preferred
        if usable(i)==1
            fprintf(sprintf('%d: %s\n', i, fn{i}))
            warning off

            a0 = [];
            
            if strcmp(Monk, 'ZZ')
                if      i==93    L = mbins>=0.55  & mbins<=0.7;
                elseif  i==92    L = mbins>=0.4   & mbins<=0.75;
                elseif  i==91    L = mbins>=0.35  & mbins<=0.6;
                elseif  i==90    L = mbins>=0.45  & mbins<=0.6;
                elseif  i==89    L = mbins>=0.65  & mbins<=0.8;
                elseif  i==88    L = mbins>=0.4   & mbins<=0.6;
                elseif  i==86    L = mbins>=0.3   & mbins<=0.4;
                elseif  i==85    L = mbins>=0.3   & mbins<=0.5;
                elseif  i==82    L = mbins>=0.65  & mbins<=0.8;
                elseif  i==81    L = mbins>=0.5   & mbins<=0.65;
                elseif  i==79    L = mbins>=0.6   & mbins<=1.5;
                elseif  i==78    L = mbins>=0.4   & mbins<=0.5;
                elseif  i==77    L = mbins>=0.45  & mbins<=0.6;
                elseif  i==76    L = mbins>=0.55  & mbins<=0.75;
                elseif  i==74    L = mbins>=0.55  & mbins<=0.7;
                elseif  i==72    L = mbins>=0.35  & mbins<=0.55;
                elseif  i==71    L = mbins>=0.5   & mbins<=0.9;
                elseif  i==70    L = mbins>=0.9   & mbins<=1.2;
                elseif  i==69    L = mbins>=0.5   & mbins<=0.8;
                elseif  i==68    L = mbins>=0.5   & mbins<=0.7;
                elseif  i==67    L = mbins>=0.35  & mbins<=0.6;
                elseif  i==66    L = mbins>=0.45  & mbins<=0.7;
                elseif  i==65    L = mbins>=0.45  & mbins<=0.65;
                elseif  i==64    L = mbins>=0     & mbins<=1.5;
                elseif  i==63    L = mbins>=0.7   & mbins<=0.9;
                elseif  i==62    L = mbins>=0.6   & mbins<=1;
                elseif  i==59    L = mbins>=0.7   & mbins<=1.3;
                elseif  i==58    L = mbins>=0.5   & mbins<=0.9;
                elseif  i==57    L = mbins>=0.5   & mbins<=0.7;
                elseif  i==56    L = mbins>=0.5   & mbins<=0.65;
                elseif  i==53    L = mbins>=0.1   & mbins<=0.6;
                elseif  i==52    L = mbins>=0.9   & mbins<=1;
                elseif  i==51    L = mbins>=0.8   & mbins<=1.2;
                elseif  i==50    L = mbins>=0.9   & mbins<=1;
                elseif  i==49    L = mbins>=0.2   & mbins<=0.5;
                elseif  i==48    L = mbins>=0     & mbins<=1;
                elseif  i==47    L = mbins>=0.55  & mbins<=0.7;
                elseif  i==45    L = mbins>=0.6   & mbins<=1;
                elseif  i==44    L = mbins>=0.3   & mbins<=0.5;
                elseif  i==43    L = mbins>=1     & mbins<=1.2;
                elseif  i==42    L = mbins>=0.4   & mbins<=0.8;
                elseif  i==40    L = mbins>=0.75  & mbins<=0.95;
                elseif  i==39    L = mbins>=0.45  & mbins<=1;
                elseif  i==38    L = mbins>=0.75  & mbins<=0.9;
                elseif  i==37    L = mbins>=0.4   & mbins<=1.2;
                elseif  i==36    L = mbins>=0.6   & mbins<=1.2;
                elseif  i==35    L = mbins>=0.6   & mbins<=1.2;
                elseif  i==34    L = mbins>=0.9   & mbins<=1.3;
                elseif  i==33    L = mbins>=0.75  & mbins<=1.15;
                elseif  i==32    L = mbins>=0.5   & mbins<=1.5;
                elseif  i==31    L = mbins>=0.1   & mbins<=1.3;
                elseif  i==30    L = mbins>=0.35  & mbins<=0.7;
                elseif  i==29    L = mbins>=0.9   & mbins<=1.3;
                elseif  i==28    L = mbins>=0.7   & mbins<=1.2;
                elseif  i==27    L = mbins>=0.4   & mbins<=0.9;
                elseif  i==26    L = mbins>=0.3   & mbins<=1;
                elseif  i==25    L = mbins>=0.2   & mbins<=1;
                elseif  i==24    L = mbins>=0.7   & mbins<=1.3;
                elseif  i==23    L = mbins>=0.3   & mbins<=0.8;
                elseif  i==22    L = mbins>=0.55  & mbins<=1.5;
                elseif  i==20    L = mbins>=0.3   & mbins<=1;
                elseif  i==18    L = mbins>=0.1   & mbins<=0.7;
                elseif  i==17    L = mbins>=0.2   & mbins<=0.8;
                elseif  i==16    L = mbins>=1     & mbins<=1.5;
                elseif  i==15    L = mbins>=0.2   & mbins<=1.1;
                elseif  i==14    L = mbins>=0.2   & mbins<=1.3;
                elseif  i==13    L = mbins>=0.5   & mbins<=1.4;
                elseif  i==12    L = mbins>=0.5   & mbins<=1.5;
                elseif  i==11    L = mbins>=0.2   & mbins<=0.6;
                elseif  i==10    L = mbins>=0.1   & mbins<=0.8;
                elseif  i==9     L = mbins>=0.4   & mbins<=0.7;
                elseif  i==8     L = mbins>=0.4   & mbins<=0.8;
                elseif  i==7     L = mbins>=1     & mbins<=1.5;
                elseif  i==6     L = mbins>=0.7   & mbins<=1.2;
                elseif  i==5     L = mbins>=0.7   & mbins<=1.5;
                elseif  i==4     L = mbins>=0.9   & mbins<=1.5;
                elseif  i==3     L = mbins>=0.5   & mbins<=1.2;
                elseif  i==2     L = mbins>=0     & mbins<=1.5;    
                elseif  i==1     L = mbins>=0     & mbins<=1.5;
                    
                else
                    L = mbins>=0.4 & mbins<=0.6;  % Zsa Zsa's time
                end
            end

            if strcmp(Monk, 'Cy')   % pick the best time period to fit data
                if     i==149   L = mbins>=0.125  & mbins<=0.2;
                elseif i==147   L = mbins>=0.2    & mbins<=0.3;
                elseif i==146   L = mbins>=0.2    & mbins<=0.35;
                elseif i==145   L = mbins>=0.15   & mbins<=0.35;
                elseif i==143   L = mbins>=0.2    & mbins<=0.4;
                elseif i==142   L = mbins>=0.1    & mbins<=0.25;
                elseif i==141   L = mbins>=0.15   & mbins<=0.25; 
                elseif i==140   L = mbins>=0.2    & mbins<=0.4; 
                elseif i==139   L = mbins>=0.15   & mbins<=0.3; 
                elseif i==138   L = mbins>=0.1    & mbins<=0.3;
                elseif i==137   L = mbins>=0.15   & mbins<=0.35;
                elseif i==136   L = mbins>=0.15   & mbins<=0.35;
                elseif i==135   L = mbins>=0.1    & mbins<=0.25;
                elseif i==131   L = mbins>=0.25   & mbins<=0.45; 
                elseif i==130   L = mbins>=0.1    & mbins<=0.4;
                elseif i==129   L = mbins>=0.2    & mbins<=0.4;
                elseif i==128   L = mbins>=0.15   & mbins<=0.4;
                elseif i==127   L = mbins>=0.15   & mbins<=0.3;
                elseif i==126   L = mbins>=0.15   & mbins<=0.3;
                elseif i==125   L = mbins>=0.15   & mbins<=0.35;
                elseif i==124   L = mbins>=0.15   & mbins<=0.35; 
                elseif i==123   L = mbins>=0.2    & mbins<=0.5; 
                elseif i==122   L = mbins>=0.15   & mbins<=0.4; 
                elseif i==120   L = mbins>=0.2    & mbins<=0.35;  
                elseif i==119   L = mbins>=0.25   & mbins<=0.35;  
                elseif i==117   L = mbins>=0.2    & mbins<=0.5;  
                elseif i==116   L = mbins>=0.15   & mbins<=0.4;  
                elseif i==114   L = mbins>=0.1    & mbins<=0.3;  
                elseif i==113   L = mbins>=0.05   & mbins<=0.3;  
                elseif i==112   L = mbins>=0.2    & mbins<=0.5;  
                elseif i==111   L = mbins>=0.2    & mbins<=0.4; 
                elseif i==110   L = mbins>=0.2    & mbins<=0.4; 
                elseif i==109   L = mbins>=0.1    & mbins<=0.3; 
                elseif i==108   L = mbins>=0.8    & mbins<=1;   
                elseif i==107   L = mbins>=0.55   & mbins<=0.8; 
                elseif i==106   L = mbins>=0.15   & mbins<=0.425;
                elseif i==104   L = mbins>=0.2    & mbins<=0.3; 
                elseif i==103   L = mbins>=0.2    & mbins<=0.5; 
                elseif i==102   L = mbins>=0.2    & mbins<=0.3; 
                elseif i==101   L = mbins>=0.2    & mbins<=0.35;
                elseif i==100   L = mbins>=0.15   & mbins<=0.25;
                elseif i==98    L = mbins>=0.15   & mbins<=0.25;
                elseif i==97    L = mbins>=0.2    & mbins<=0.4;
                elseif i==96    L = mbins>=0.15   & mbins<=0.3;
                elseif i==95    L = mbins>=0.15   & mbins<=0.35;
                elseif i==94    L = mbins>=0.15   & mbins<=0.3;
                elseif i==93    L = mbins>=0.075  & mbins<=0.35;
                elseif i==92    L = mbins>=0.15   & mbins<=0.4;
                elseif i==91    L = mbins>=0.15   & mbins<=0.3;
                elseif i==90    L = mbins>=0.2    & mbins<=0.35;
                elseif i==89    L = mbins>=0.15   & mbins<=0.4;
                elseif i==88    L = mbins>=0.15   & mbins<=0.35;
                elseif i==87    L = mbins>=0.15   & mbins<=0.3;
                elseif i==86    L = mbins>=0.15   & mbins<=0.25;
                elseif i==85    L = mbins>=0.15   & mbins<=0.25;
                elseif i==84    L = mbins>=0.1    & mbins<=0.25;
                elseif i==83    L = mbins>=0.15   & mbins<=0.3;
                elseif i==82    L = mbins>=0.2    & mbins<=0.3;
                elseif i==81    L = mbins>=0.15   & mbins<=0.3;
                elseif i==79    L = mbins>=0.15   & mbins<=0.4;
                elseif i==78    L = mbins>=0.3    & mbins<=0.5;
                elseif i==77    L = mbins>=0.2    & mbins<=0.6;
                elseif i==76    L = mbins>=0.15   & mbins<=0.35;
                elseif i==75    L = mbins>=0.15   & mbins<=0.3;
                elseif i==74    L = mbins>=0.1    & mbins<=0.4;
                elseif i==70    L = mbins>=0.1    & mbins<=0.4;
                elseif i==69    L = mbins>=0.15   & mbins<=0.35;
                elseif i==66    L = mbins>=0.2    & mbins<=0.4;
                elseif i==65    L = mbins>=0.5    & mbins<=0.7;
                elseif i==64    L = mbins>=0.15   & mbins<=0.35;
                elseif i==61    L = mbins>=0.15   & mbins<=0.35;
                elseif i==60    L = mbins>=0.35   & mbins<=1;
                elseif i==59    L = mbins>=0.15   & mbins<=0.35;
                elseif i==58    L = mbins>=0.15   & mbins<=0.35;
                elseif i==57    L = mbins>=0      & mbins<=1.5;
                elseif i==56    L = mbins>=0      & mbins<=0.7;
                elseif i==55    L = mbins>=0.6    & mbins<=0.9;
                elseif i==53    L = mbins>=0.3    & mbins<=0.7;
                elseif i==52    L = mbins>=0      & mbins<=0.5;
                elseif i==51    L = mbins>=0.5    & mbins<=1;
                elseif i==49    L = mbins>=0      & mbins<=0.6;
                elseif i==48    L = mbins>=0.15   & mbins<=0.4;
                elseif i==47    L = mbins>=0.4    & mbins<=0.65;
                elseif i==46    L = mbins>=0.1    & mbins<=1.5;
                elseif i==45    L = mbins>=0.15   & mbins<=0.5;
                elseif i==44    L = mbins>=0.3    & mbins<=0.45;
                elseif i==43    L = mbins>=0.1    & mbins<=0.8;
                elseif i==42    L = mbins>=0.1    & mbins<=0.7;
                elseif i==41    L = mbins>=0.2    & mbins<=0.6;
                elseif i==40    L = mbins>=0      & mbins<=0.5;
                elseif i==38    L = mbins>=0.2    & mbins<=0.6;
                elseif i==37    L = mbins>=0.4    & mbins<=1;
                elseif i==35    L = mbins>=0.1    & mbins<=0.5;
                elseif i==34    L = mbins>=0.8    & mbins<=1.4;
                elseif i==33    L = mbins>=0.2    & mbins<=0.4;
                elseif i==32    L = mbins>=0.8    & mbins<=1;
                elseif i==30    L = mbins>=0.3    & mbins<=0.6;
                elseif i==28    L = mbins>=0.25   & mbins<=0.6;
                elseif i==24    L = mbins>=0.15   & mbins<=0.4;
                elseif i==23    L = mbins>=0      & mbins<=0.5;
                elseif i==22    L = mbins>=0.1    & mbins<=0.4;
                elseif i==21    L = mbins>=0.1    & mbins<=0.4;
                elseif i==19    L = mbins>=0.2    & mbins<=0.5;
                elseif i==18    L = mbins>=0.2    & mbins<=0.8;
                elseif i==17    L = mbins>=0.4    & mbins<=1;
                elseif i==16    L = mbins>=0      & mbins<=0.4;
                elseif i==14    L = mbins>=0      & mbins<=0.5;
                elseif i==13    L = mbins>=0.1    & mbins<=0.4;
                elseif i==12    L = mbins>=0.2    & mbins<=0.45;
                elseif i==8     L = mbins>=0.25   & mbins<=0.55;
                elseif i==7     L = mbins>=0      & mbins<=0.8;
                elseif i==6     L = mbins>=0.1    & mbins<=0.5;
                elseif i==5     L = mbins>=0.3    & mbins<=0.8;
                elseif i==4     L = mbins>=0      & mbins<=1.5;
                elseif i==3     L = mbins>=0.35   & mbins<=0.7;
                elseif i==2     L = mbins>=0.1    & mbins<=0.5;
                elseif i==1     L = mbins>=0.25   & mbins<=0.7;

                else
                    L  = mbins>=0.075 & mbins<=0.35;  % cyrus's time
                end
            end
            if isempty(a0)
                a0 = nan;
            end


            
            if strcmp(Monk, 'Cy')
                bcon = [-0.1 -0.5 min(mbins(L)) 0.5; 30 0.4 max(mbins(L)) 1];
            else
                bcon = [-0.1 -0.2 min(mbins(L)) 0.5; 20 0.4 max(mbins(L)) 1];
            end

            init = [a0 nan nan nan];
            
            coh   = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
            mbP   = mbins(L);
            mbP   = repmat(mbP',sum(Lcoh),1);
            mbP   = mbP(:);
            c     = repmat(coh(Lcoh),1,sum(L));
            c     = c(:);
            r     = m(L,Lcoh,i)';
            rd    = sd(L,Lcoh,i)';
            r     = r(:);
            rd    = rd(:);

            citype = {10, 68};
            [a(:,i), ae(:,:,i), yp] = getML_fitLIPslope7(mbP, c, r, rd, init, bcon, citype);
            l(i)                    = yp(mbP==nanmax(mbP)&c==nanmax(c)); % 'lapse rate' predicted by model

            a(:,i)
            ae(:,:,i)
            % plot and check fit
            if 1
                clf
                lc  = [];
                for j = 1:7
                    lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
                end

                for j = find(Lcoh)
                    subplot(4,2,j)
                    hold on
                    plot(mbins, m(:,j,i), 'o', 'MarkerSize', 6, ...
                        'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
                    Lc = c==coh(j);
                    plot(mbP(Lc), yp(Lc), '-', 'LineWidth', 3, 'Color', lc(j,:))
                    xlim([0 1.5]);
                    ylim([0.4 1]);
                end
                title(sprintf('%.2f (%.2f,%.2f)', a(1,i), ae(1,1,i), ae(1,2,i)))

                %                 subplot(4,2,8)
                %                 plot(coh, slope(:,1,i), 'ro', 'MarkerSize', 6, ...
                %                                         'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
                %                 line(repmat(coh',2,1), [(slope(:,1,i)-slope(:,2,i))'; (slope(:,1,i)+slope(:,2,i))'],...
                %                                 'Color', 'r', 'LineWidth', 1)
                %                 line(coh, a(1,i)+coh.*a(2,i), 'Color', 'k', 'LineWidth', 2)
                %                 title(sprintf('%.3f', a(2,i)))
                %
                %                 ylim([min(slope(:,1,i)), max(slope(:,1,i))])
                %
                pause
            end
        end
    end


    save(savepath, 'a', 'ae', 'l')
else

    load(savepath)
end

