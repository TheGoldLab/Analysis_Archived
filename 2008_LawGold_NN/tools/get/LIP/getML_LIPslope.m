function [f fi s] = getML_LIPslope(Monk, recompute)
% getML_LIPslope
% get slope of LIP responses for each coherence

[hdir, ldir, cdir, tdir] = dirnames;
savepath = [tdir '/' Monk '.mat'];


if recompute
    %% load data
    a      = getML_txt([Monk 'TRain_LIP.txt']);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    ses    = a.data{strcmp(a.name,'session')};
    usable = a.data{strcmp(a.name,'usable')};

    % get ROC areas
    crtf = [1];
    nvf  = [0 1];
    dirf = 0;

    bb     = -200;
    be     = 1500;
    bs     = 25;
    bw     = 50;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    mbins  = mean(bins,2)/1000;

    [rM, rSD, rN, rMn] = getML_mRate([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, 0);

  
    b  = nans(7,length(fn));
    bs = nans(7,length(fn));

    m    = rMn(1:7,:,:)-rMn(8:14,:,:);            % select between normalized and unnormalized spikes
    sd   = (1./sqrt(rN(1:7,:,:)))+(1./sqrt(rN(8:14,:,:)));         % weight by # of trials because I don't have sd for normalized spike rate
    coh  = [0 0.032 0.064 0.128 0.256 0.512 0.999]';

    lc  = [];
    for i = 1:7
        lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
    end


    for i = 1:length(fn)
        if usable(i)==1
            fprintf(sprintf('%d: %s\n', i, fn{i}))

            % ramp on ramp off time estimated by the piecewise-linear
            % function
              if strcmp(Monk, 'ZZ')
                if      i==95    L = mbins>=0.175  & mbins<=0.3;
                elseif  i==93    L = mbins>=0.55  & mbins<=0.9;
                elseif  i==92    L = mbins>=0.5   & mbins<=0.7;
                elseif  i==91    L = mbins>=0.3   & mbins<=0.6;
                elseif  i==90    L = mbins>=0.4   & mbins<=0.55;
                elseif  i==89    L = mbins>=0.55  & mbins<=0.75;
                elseif  i==88    L = mbins>=0.4   & mbins<=0.7;
                elseif  i==87    L = mbins>=0.5   & mbins<=0.775;
                elseif  i==86    L = mbins>=0.25  & mbins<=0.8;
                elseif  i==85    L = mbins>=0.3   & mbins<=0.45;
                elseif  i==82    L = mbins>=1     & mbins<=1.1;
                elseif  i==81    L = mbins>=0.4   & mbins<=0.95;
                elseif  i==79    L = mbins>=0.5   & mbins<=0.75;
                elseif  i==78    L = mbins>=0.4   & mbins<=0.8;
                elseif  i==77    L = mbins>=0.4   & mbins<=1.2;
                elseif  i==76    L = mbins>=0.55  & mbins<=0.9;
                elseif  i==74    L = mbins>=0.2   & mbins<=0.45;
                elseif  i==72    L = mbins>=0.3   & mbins<=0.55;
                elseif  i==71    L = mbins>=0.45  & mbins<=0.8;
                elseif  i==70    L = mbins>=0.35  & mbins<=0.8;
                elseif  i==69    L = mbins>=0.3   & mbins<=0.7;
                elseif  i==68    L = mbins>=0.5   & mbins<=0.7;
                elseif  i==67    L = mbins>=0.35  & mbins<=0.65;
                elseif  i==66    L = mbins>=0.2   & mbins<=0.7;
                elseif  i==65    L = mbins>=0.35  & mbins<=0.65;
                elseif  i==64    L = mbins>=0.4   & mbins<=1.3;
                elseif  i==63    L = mbins>=0.6   & mbins<=0.95;
                elseif  i==62    L = mbins>=0.375 & mbins<=0.55;
                elseif  i==59    L = mbins>=0.1   & mbins<=0.35;
                elseif  i==58    L = mbins>=0.5   & mbins<=0.8;
                elseif  i==57    L = mbins>=0.675 & mbins<=0.9;
                elseif  i==56    L = mbins>=0.5   & mbins<=0.675;
                elseif  i==53    L = mbins>=0.7   & mbins<=1.1;
                elseif  i==52    L = mbins>=0.7   & mbins<=1.1;
                elseif  i==51    L = mbins>=0.4   & mbins<=1;
                elseif  i==50    L = mbins>=0.6   & mbins<=0.8;
                elseif  i==49    L = mbins>=0.7   & mbins<=1;
                elseif  i==48    L = mbins>=0.2   & mbins<=0.5;
                elseif  i==47    L = mbins>=0.35  & mbins<=0.7;
                elseif  i==45    L = mbins>=0.7   & mbins<=1;
                elseif  i==44    L = mbins>=0.3   & mbins<=0.5;
                elseif  i==43    L = mbins>=0.7   & mbins<=1.2;
                elseif  i==42    L = mbins>=0.8   & mbins<=1.05;
                elseif  i==40    L = mbins>=0.4   & mbins<=0.95;
                elseif  i==39    L = mbins>=0.45  & mbins<=1;
                elseif  i==38    L = mbins>=0.4   & mbins<=0.9;
                elseif  i==37    L = mbins>=0.4   & mbins<=1.2;
                elseif  i==36    L = mbins>=0.1   & mbins<=0.3;
                elseif  i==35    L = mbins>=0.4   & mbins<=1.5;
                elseif  i==34    L = mbins>=0     & mbins<=1.5;
                elseif  i==33    L = mbins>=0.1   & mbins<=1.15;
                elseif  i==32    L = mbins>=0.1   & mbins<=1.2;
                elseif  i==31    L = mbins>=0.1   & mbins<=1.5;
                elseif  i==30    L = mbins>=0     & mbins<=1.5;
                elseif  i==29    L = mbins>=0.3   & mbins<=0.5;
                elseif  i==28    L = mbins>=0     & mbins<=1.2;
                elseif  i==27    L = mbins>=0     & mbins<=1.5;
                elseif  i==26    L = mbins>=0     & mbins<=1.5;
                elseif  i==25    L = mbins>=0     & mbins<=1.5;
                elseif  i==24    L = mbins>=0     & mbins<=1.5;
                elseif  i==23    L = mbins>=0     & mbins<=1.5;
                elseif  i==22    L = mbins>=0     & mbins<=1.5;
                elseif  i==20    L = mbins>=0     & mbins<=1;
                elseif  i==18    L = mbins>=0     & mbins<=1.5;
                elseif  i==17    L = mbins>=0.2   & mbins<=1;
                elseif  i==16    L = mbins>=0.1   & mbins<=1.5;
                elseif  i==15    L = mbins>=0     & mbins<=1;
                elseif  i==14    L = mbins>=0     & mbins<=1.3;
                elseif  i==13    L = mbins>=0.1   & mbins<=1;
                elseif  i==12    L = mbins>=0.1   & mbins<=0.9;
                elseif  i==11    L = mbins>=0.1   & mbins<=1.5;
                elseif  i==10    L = mbins>=0     & mbins<=1.5;
                elseif  i==9     L = mbins>=0     & mbins<=1.5;
                elseif  i==8     L = mbins>=0     & mbins<=1.5;
                elseif  i==7     L = mbins>=0.6   & mbins<=0.9;
                elseif  i==6     L = mbins>=0.1   & mbins<=1.5;
                elseif  i==5     L = mbins>=0     & mbins<=1.2;
                elseif  i==4     L = mbins>=0.1   & mbins<=1.5;
                elseif  i==3     L = mbins>=0.8   & mbins<=1.5;
                elseif  i==2     L = mbins>=0     & mbins<=1.5;    
                elseif  i==1     L = mbins>=0.7   & mbins<=1.5;

                else
                    L = mbins>=0.4 & mbins<=0.6;  % Zsa Zsa's time
                end
            end

            if strcmp(Monk, 'Cy')   % pick the best time period to fit data
                if     i==149   L = mbins>=0.15   & mbins<=0.25;
                elseif i==147   L = mbins>=0.15   & mbins<=0.225;
                elseif i==146   L = mbins>=0.225  & mbins<=0.375;
                elseif i==145   L = mbins>=0.2    & mbins<=0.4;
                elseif i==143   L = mbins>=0.15   & mbins<=0.4;
                elseif i==142   L = mbins>=0.125  & mbins<=0.325;
                elseif i==141   L = mbins>=0.1    & mbins<=0.2;
                elseif i==140   L = mbins>=0.15   & mbins<=0.35;
                elseif i==139   L = mbins>=0.2    & mbins<=0.375;
                elseif i==138   L = mbins>=0.125  & mbins<=0.225;
                elseif i==137   L = mbins>=0.15   & mbins<=0.25;
                elseif i==136   L = mbins>=0.15   & mbins<=0.3;
                elseif i==135   L = mbins>=0.125  & mbins<=0.275;
                elseif i==131   L = mbins>=0.25   & mbins<=0.45;
                elseif i==130   L = mbins>=0.175  & mbins<=0.3;
                elseif i==129   L = mbins>=0.2    & mbins<=0.4;
                elseif i==128   L = mbins>=0.15   & mbins<=0.25;
                elseif i==127   L = mbins>=0.15   & mbins<=0.275;
                elseif i==126   L = mbins>=0.15   & mbins<=0.25;
                elseif i==125   L = mbins>=0.15   & mbins<=0.375;
                elseif i==124   L = mbins>=0.2    & mbins<=0.4;
                elseif i==123   L = mbins>=0.3    & mbins<=0.7;
                elseif i==122   L = mbins>=0.1    & mbins<=0.25;
                elseif i==120   L = mbins>=0.125  & mbins<=0.45;
                elseif i==119   L = mbins>=0.225  & mbins<=0.4;
                elseif i==117   L = mbins>=0.1    & mbins<=0.8;
                elseif i==116   L = mbins>=0.2    & mbins<=0.3;
                elseif i==114   L = mbins>=0.1    & mbins<=0.25;
                elseif i==113   L = mbins>=0.05   & mbins<=0.3;
                elseif i==112   L = mbins>=0.1    & mbins<=0.325;
                elseif i==111   L = mbins>=0.2    & mbins<=0.4;
                elseif i==110   L = mbins>=0.15   & mbins<=0.3;
                elseif i==109   L = mbins>=0.1    & mbins<=0.3;
                elseif i==108   L = mbins>=0.75   & mbins<=1;
                elseif i==107   L = mbins>=0.55   & mbins<=0.8;
                elseif i==106   L = mbins>=0.15   & mbins<=0.35;
                elseif i==104   L = mbins>=0      & mbins<=0.6;
                elseif i==103   L = mbins>=0.2    & mbins<=0.4;
                elseif i==102   L = mbins>=0.15   & mbins<=0.4;
                elseif i==101   L = mbins>=0.2    & mbins<=0.4;
                elseif i==100   L = mbins>=0.15   & mbins<=0.3;
                elseif i==98    L = mbins>=0.15   & mbins<=0.25;
                elseif i==97    L = mbins>=0.2    & mbins<=0.4;
                elseif i==96    L = mbins>=0.15   & mbins<=0.3;
                elseif i==95    L = mbins>=0.15   & mbins<=0.35;
                elseif i==94    L = mbins>=0.15   & mbins<=0.3;
                elseif i==93    L = mbins>=0.2    & mbins<=0.6;
                elseif i==92    L = mbins>=0.15   & mbins<=0.4;
                elseif i==91    L = mbins>=0.15   & mbins<=0.25;
                elseif i==90    L = mbins>=0.2    & mbins<=0.7;
                elseif i==89    L = mbins>=0.1    & mbins<=0.4;
                elseif i==88    L = mbins>=0.125   & mbins<=0.35;
                elseif i==87    L = mbins>=0.15   & mbins<=0.3;
                elseif i==86    L = mbins>=0.15   & mbins<=0.25;
                elseif i==85    L = mbins>=0.15   & mbins<=0.5;
                elseif i==84    L = mbins>=0.1    & mbins<=0.25;
                elseif i==83    L = mbins>=0.15   & mbins<=0.3;
                elseif i==82    L = mbins>=0.2    & mbins<=0.3;
                elseif i==81    L = mbins>=0.15   & mbins<=0.4;
                elseif i==79    L = mbins>=0.125  & mbins<=0.35;
                elseif i==78    L = mbins>=0.35   & mbins<=0.5;
                elseif i==77    L = mbins>=0.2    & mbins<=0.6;
                elseif i==76    L = mbins>=0.15   & mbins<=0.3;
                elseif i==75    L = mbins>=0.15   & mbins<=0.25;
                elseif i==74    L = mbins>=0.1    & mbins<=0.4;
                elseif i==70    L = mbins>=0.1    & mbins<=0.4;
                elseif i==69    L = mbins>=0.1    & mbins<=0.4;
                elseif i==66    L = mbins>=0.15   & mbins<=0.35;
                elseif i==65    L = mbins>=0.5    & mbins<=0.7;
                elseif i==64    L = mbins>=0.15   & mbins<=0.35;
                elseif i==61    L = mbins>=0.1    & mbins<=0.275;
                elseif i==60    L = mbins>=0      & mbins<=1;
                elseif i==59    L = mbins>=0.15   & mbins<=0.4;
                elseif i==58    L = mbins>=0.15   & mbins<=0.35;
                elseif i==57    L = mbins>=0      & mbins<=1;
                elseif i==56    L = mbins>=0      & mbins<=0.7;
                elseif i==55    L = mbins>=0.6    & mbins<=0.9;
                elseif i==53    L = mbins>=0.3    & mbins<=0.7;
                elseif i==52    L = mbins>=0.2    & mbins<=0.5;
                elseif i==51    L = mbins>=0.5    & mbins<=1;
                elseif i==49    L = mbins>=0      & mbins<=0.4;
                elseif i==48    L = mbins>=0.1   & mbins<=0.4;
                elseif i==47    L = mbins>=0.4    & mbins<=0.65;
                elseif i==46    L = mbins>=0.1    & mbins<=1.5;
                elseif i==45    L = mbins>=0.15   & mbins<=0.5;
                elseif i==44    L = mbins>=0.3    & mbins<=0.45;
                elseif i==43    L = mbins>=0      & mbins<=0.8;
                elseif i==42    L = mbins>=0.1    & mbins<=0.7;
                elseif i==41    L = mbins>=0.2    & mbins<=0.6;
                elseif i==40    L = mbins>=0      & mbins<=0.5;
                elseif i==38    L = mbins>=0.2    & mbins<=0.6;
                elseif i==37    L = mbins>=0.4    & mbins<=1;
                elseif i==35    L = mbins>=0.075  & mbins<=0.5;
                elseif i==34    L = mbins>=0      & mbins<=1.4;
                elseif i==33    L = mbins>=0.1    & mbins<=0.6;
                elseif i==32    L = mbins>=0.6    & mbins<=1;
                elseif i==31    L = mbins>=0.1    & mbins<=0.6;
                elseif i==30    L = mbins>=0.3    & mbins<=0.6;
                elseif i==28    L = mbins>=0.25   & mbins<=0.6;
                elseif i==24    L = mbins>=0.15   & mbins<=0.4;
                elseif i==23    L = mbins>=0      & mbins<=0.5;
                elseif i==22    L = mbins>=0.1    & mbins<=0.4;
                elseif i==21    L = mbins>=0.1    & mbins<=0.4;
                elseif i==19    L = mbins>=0.2    & mbins<=0.5;
                elseif i==18    L = mbins>=0.2    & mbins<=0.8;
                elseif i==17    L = mbins>=0.4    & mbins<=1;
                elseif i==16    L = mbins>=0.125  & mbins<=0.45;
                elseif i==14    L = mbins>=0      & mbins<=0.5;
                elseif i==13    L = mbins>=0.1    & mbins<=0.4;
                elseif i==12    L = mbins>=0.2    & mbins<=0.6;
                elseif i==8     L = mbins>=0.25   & mbins<=0.55;
                elseif i==7     L = mbins>=0      & mbins<=0.8;
                elseif i==6     L = mbins>=0.1    & mbins<=0.5;
                elseif i==5     L = mbins>=0.3    & mbins<=0.8;
                elseif i==4     L = mbins>=0      & mbins<=1.5;
                elseif i==3     L = mbins>=0.0   & mbins<=0.7;
                elseif i==2     L = mbins>=0.2    & mbins<=0.6;
                elseif i==1     L = mbins>=0.3    & mbins<=0.625;

                else
                    L  = mbins>=0.075 & mbins<=0.35;  % cyrus's time
                end
            end
            %
            %         % fit data to high coh to estimate onset and offset time
            %         if nanmax(n(:,7,i))>nanmax(n(:,6,i))
            %             [bb, bsem, gof, ym] = sramp_fitW(mbins(L), m(L,7,i), sd(L,7,i), [], [], ...
            %                 [min(mbins(L)) 0 min(mbins(L)); max(mbins(L)) 10 max(mbins(L))]);
            %             Lt = mbins>=round(20*bb(1))/20 & mbins<=round(20*bb(3))/20;
            %             if 1
            %                 clf
            %                 subplot(2,4,8)
            %                 hold on
            %                 plot(mbins, m(:,7,i), 'o', 'MarkerSize', 6, ...
            %                     'MarkerFaceColor', lc(7,:), 'MarkerEdgeColor', 'none')
            %                 plot(mbins(L),ym,'k', 'LineWidth', 3)
            %                 hold off
            %                 title(sprintf('%.2f - %.2f', bb(1), bb(3)))
            %                 xlim([0 1.5]);
            %                 ylim([0.4 1]);
            %             end
            %
            %         else
            %             [bb, bsem, gof, ym] = sramp_fitW(mbins(L), m(L,6,i), sd(L,6,i), [], [], ...
            %                 [min(mbins(L)) 0 min(mbins(L)); max(mbins(L)) 10 max(mbins(L))]);
            %             Lt = mbins>=round(20*bb(1))/20 & mbins<=round(20*bb(3))/20;
            %             bb = nans(2,7);
            %             if 1
            %                 clf
            %                 subplot(2,4,8)
            %                 hold on
            %                 plot(mbins, m(:,6,i), 'o', 'MarkerSize', 6, ...
            %                     'MarkerFaceColor', lc(7,:), 'MarkerEdgeColor', 'none')
            %                 plot(mbins(L),ym,'k', 'LineWidth', 3)
            %                 hold off
            %                 title(sprintf('%.2f - %.2f', bb(1), bb(3)))
            %                 xlim([0 1.5]);
            %                 ylim([0.4 1]);
            %             end
            %
            %         end



            bb = nans(2,7);
            for j = 1:7
                Lgd = ~isnan(m(j,:,i))& ~isnan(sd(j,:,i)) &sd(j,:,i)~=0;
                Lt  = L';
                if sum(Lgd&Lt)>3
                    sd(sd(j,:,i)==0,j,i) = nanmean(sd(j,:,i));
                    [bb(:,j), bbs] = regressW(m(j,Lt&Lgd,i)', sd(j,Lt&Lgd,i)', [ones(sum(Lt&Lgd),1) mbins(Lt&Lgd)]);
                    b(j,i)    = bb(2,j);
                    bs(j,i)   = bbs(2,2)-bb(2,j);
                end
            end


            if 0
                clf
                for j = 1:7
                    subplot(2,4,j)
                    hold on
                    plot(mbins, m(j,:,i), 'o', 'MarkerSize', 6, ...
                        'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
                    plot(mbins(L), bb(1,j)+bb(2,j)*mbins(L), '-', 'LineWidth', 3, 'Color', 'r')
                    hold off
                    xlim([0 1.5]);
                    ylim([-1 2]);

                    title(sprintf('%.2f',coh(j)))
                end
                pause
            end
        end
    end

    for j = 1:7
        L   = ~isnan(b(j,:));
        
        if strcmp(Monk, 'Cy')   % some of cyrus's data has really big or really small sd, so don't do weight fit
            [fit fiti h p] = nestedFW(b(j,L)',ones(size(bs(j,L)')),[ones(sum(L),1) ses(L)], 1-0.68);
        else
            [fit fiti h p] = nestedFW(b(j,L)',bs(j,L)',[ones(sum(L),1) ses(L)], 1-0.68);
        end
        
        f(j)  = fit(2);
        fi(j) = fiti(2,2)-fit(2);
        s(j)  = p;
    end


    save(savepath, 'f', 'fi', 's')

else
    load(savepath)
end


if 0

    figure
    hold on
    plot([1 3.2 6.4 12.8 25.6 51.2 99.9],f,'o','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
    plot([[1 3.2 6.4 12.8 25.6 51.2 99.9];[1 3.2 6.4 12.8 25.6 51.2 99.9]],[f-fi; f+fi],'r-')
    hold off
    set(gca, 'xscale', 'log', 'xlim', [1 99.9], ...
        'XTick', [1 3.2 6.4 12.8 25.6 51.2 99.9], 'XTickLabel', [0 3.2 6.4 12.8 25.6 51.2 99.9])

end










% %% old code, fit to ROC area
% function [f fi s] = getML_LIPslope(Monk, recompute)
% % getML_LIPslope
% % get slope of LIP responses for each coherence
% 
% 
% savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_LIPslope_' Monk '.mat'];
% 
% 
% if recompute
%     %% load data
%     a      = getML_txt([Monk 'TRain_LIP.txt']);
%     fn     = a.data{strcmp(a.name,'dat_fn')};
%     ses    = a.data{strcmp(a.name,'session')};
%     usable = a.data{strcmp(a.name,'usable')};
% 
%     % get ROC areas
%     crtf = [1];
%     nvf  = [0 1];
%     dirf = 0;
% 
%     bb     = -200;
%     be     = 1500;
%     bs     = 50;
%     bw     = 100;
%     bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
%     mbins  = mean(bins,2)/1000;
% 
%     [rs rn re] = getML_ROCs([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, [], 0);
% 
%     b  = nans(7,length(fn));
%     bs = nans(7,length(fn));
% 
%     m    = rs;
%     sd   = re./(rn.^0.5);   % standard error
%     n    = rn;
%     coh  = [0 0.032 0.064 0.128 0.256 0.512 0.999]';
% 
%     lc  = [];
%     for i = 1:7
%         lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
%     end
% 
% 
%     for i = 1:length(fn)
%         if usable(i)==1
%             fprintf(sprintf('%d: %s\n', i, fn{i}))
% 
%             % ramp on ramp off time estimated by the piecewise-linear
%             % function
%             if strcmp(Monk, 'Cy')  % cyrus's time
%                 if      i==1    L = mbins>=0.1 & mbins<=0.6;
%                 elseif  i==2    L = mbins>=0   & mbins<=1.4;
%                 elseif  i==3    L = mbins>=0.5 & mbins<=1.5;
%                 elseif  i==4    L = mbins>=0.5 & mbins<=1.5;
%                 elseif  i==5    L = mbins>=0   & mbins<=1;
%                 elseif  i==6    L = mbins>=0   & mbins<=0.6;
%                 elseif  i==7    L = mbins>=0   & mbins<=0.8;
%                 elseif  i==8    L = mbins>=0.1 & mbins<=0.55;
%                 elseif  i==12   L = mbins>=0.1 & mbins<=0.55;
%                 elseif  i==13   L = mbins>=0   & mbins<=0.5;
%                 elseif  i==14   L = mbins>=0   & mbins<=1.4;
%                 elseif  i==16   L = mbins>=0   & mbins<=0.6;
%                 elseif  i==17   L = mbins>=0.2 & mbins<=1.4;
%                 elseif  i==18   L = mbins>=0.2 & mbins<=0.9;
%                 elseif  i==19   L = mbins>=0.1 & mbins<=0.75;
%                 elseif  i==21   L = mbins>=0.1 & mbins<=1.5;
%                 elseif  i==22   L = mbins>=0.1 & mbins<=0.5;
%                 elseif  i==23   L = mbins>=0.1 & mbins<=0.5;
%                 elseif  i==24   L = mbins>=0   & mbins<=0.4;
%                 elseif  i==25   L = mbins>=0   & mbins<=0.5;
%                 elseif  i==28   L = mbins>=0   & mbins<=0.5;
%                 elseif  i==30   L = mbins>=0.1 & mbins<=0.8;
%                 elseif  i==31   L = mbins>=0   & mbins<=0.5;
%                 elseif  i==32   L = mbins>=0.3 & mbins<=0.7;
%                 elseif  i==33   L = mbins>=0   & mbins<=0.4;
%                 elseif  i==34   L = mbins>=0   & mbins<=0.6;
%                 elseif  i==35   L = mbins>=0.1 & mbins<=0.4;
%                 elseif  i==37   L = mbins>=0.2 & mbins<=1.2;
%                 elseif  i==38   L = mbins>=0.15& mbins<=0.75;
%                 elseif  i==40   L = mbins>=0.1 & mbins<=0.5;
%                 elseif  i==41   L = mbins>=0.1 & mbins<=1.1;
%                 elseif  i==42   L = mbins>=0   & mbins<=1;
%                 elseif  i==44   L = mbins>=0.1 & mbins<=1;
%                 elseif  i==45   L = mbins>=0.1 & mbins<=0.7;
%                 elseif  i==46   L = mbins>=0   & mbins<=0.3;
%                 elseif  i==47   L = mbins>=0   & mbins<=0.9;
%                 elseif  i==48   L = mbins>=0.1 & mbins<=0.4;
%                 elseif  i==49   L = mbins>=0   & mbins<=0.4;
%                 elseif  i==51   L = mbins>=0.5 & mbins<=1;
%                 elseif  i==52   L = mbins>=0   & mbins<=0.4;
%                 elseif  i==53   L = mbins>=0.4 & mbins<=0.6;
%                 elseif  i==55   L = mbins>=0.6 & mbins<=0.8;
%                 elseif  i==56   L = mbins>=0   & mbins<=0.7;
%                 elseif  i==57   L = mbins>=0   & mbins<=0.7;
%                 elseif  i==58   L = mbins>=0.1 & mbins<=0.3;
%                 elseif  i==59   L = mbins>=0.1 & mbins<=0.3;
%                 elseif  i==60   L = mbins>=0.5 & mbins<=0.9;
%                 elseif  i==61   L = mbins>=0   & mbins<=0.3;
%                 elseif  i==64   L = mbins>=0   & mbins<=0.3;
%                 elseif  i==65   L = mbins>=0.55 & mbins<=0.75;
%                 elseif  i==66   L = mbins>=0   & mbins<=0.3;
%                 elseif  i==69   L = mbins>=0.05& mbins<=0.3;
%                 elseif  i==70   L = mbins>=0.05& mbins<=0.3;
%                 elseif  i==74   L = mbins>=0.05& mbins<=0.4;
%                 elseif  i==75   L = mbins>=0.05& mbins<=0.25;
%                 elseif  i==76   L = mbins>=0.1 & mbins<=0.3;
%                 elseif  i==77   L = mbins>=0   & mbins<=1;
%                 elseif  i==78   L = mbins>=0.1 & mbins<=0.8;
%                 elseif  i==79   L = mbins>=0.1 & mbins<=0.35;
%                 elseif  i==80   L = mbins>=0.05& mbins<=0.3;
%                 elseif  i==81   L = mbins>=0   & mbins<=0.3;
%                 elseif  i==82   L = mbins>=0   & mbins<=0.3;
%                 elseif  i==83   L = mbins>=0   & mbins<=0.3;
%                 elseif  i==84   L = mbins>=0   & mbins<=0.3;
%                 elseif  i==85   L = mbins>=0   & mbins<=0.35;
%                 elseif  i==86   L = mbins>=0   & mbins<=0.35;
%                 elseif  i==87   L = mbins>=0   & mbins<=0.3;
%                 elseif  i==88   L = mbins>=0.05& mbins<=0.3;
%                 elseif  i==89   L = mbins>=0.05& mbins<=0.3;
%                 elseif  i==90   L = mbins>=0.05& mbins<=0.8;
%                 elseif  i==91   L = mbins>=0.05& mbins<=0.3;
%                 elseif  i==92   L = mbins>=0.1 & mbins<=0.35;
%                 elseif  i==93   L = mbins>=0.1 & mbins<=0.35;
%                 elseif  i==94   L = mbins>=0.1 & mbins<=0.3;
%                 elseif  i==95   L = mbins>=0.1 & mbins<=0.3;
%                 elseif  i==96   L = mbins>=0.1 & mbins<=0.3;
%                 elseif  i==97   L = mbins>=0.1 & mbins<=0.4;
%                 elseif  i==98   L = mbins>=0.05& mbins<=0.35;
%                 elseif  i==99   L = mbins>=0   & mbins<=0.4;
%                 elseif  i==100  L = mbins>=0.05& mbins<=0.3;
%                 elseif  i==101  L = mbins>=0.05& mbins<=0.3;
%                 elseif  i==102  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==104  L = mbins>=0   & mbins<=0.35;
%                 elseif  i==106  L = mbins>=0   & mbins<=0.4;
%                 elseif  i==107  L = mbins>=0.1 & mbins<=0.7;
%                 elseif  i==108  L = mbins>=0.1 & mbins<=0.4;
%                 elseif  i==109  L = mbins>=0.1 & mbins<=0.3;
%                 elseif  i==110  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==112  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==113  L = mbins>=0   & mbins<=0.25;
%                 elseif  i==114  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==116  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==119  L = mbins>=0.1 & mbins<=0.55;
%                 elseif  i==120  L = mbins>=0   & mbins<=0.35;
%                 elseif  i==122  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==123  L = mbins>=0.1 & mbins<=0.6;
%                 elseif  i==124  L = mbins>=0   & mbins<=0.4;
%                 elseif  i==125  L = mbins>=0   & mbins<=0.4;
%                 elseif  i==126  L = mbins>=0   & mbins<=0.35;
%                 elseif  i==127  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==128  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==129  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==130  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==131  L = mbins>=0.3 & mbins<=0.6;
%                 elseif  i==135  L = mbins>=0   & mbins<=0.2;
%                 elseif  i==136  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==137  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==138  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==139  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==141  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==142  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==146  L = mbins>=0.05& mbins<=0.35;
%                 elseif  i==147  L = mbins>=0   & mbins<=0.3;
%                 elseif  i==149  L = mbins>=0   & mbins<=0.2;
%                 else
%                     L = mbins>=0 & mbins<=0.6;
%                 end
% 
%             else
% 
%                 if      i==1   L = mbins>=0    & mbins<=1.5;
%                 elseif  i==8   L = mbins>=0    & mbins<=0.7;
%                 elseif  i==9   L = mbins>=0.3  & mbins<=0.8;
%                 elseif  i==15  L = mbins>=0.2  & mbins<=1.3;
%                 elseif  i==20  L = mbins>=0.2  & mbins<=1;
%                 elseif  i==27  L = mbins>=0.2  & mbins<=1.2;
%                 elseif  i==28  L = mbins>=0.4  & mbins<=1.2;
%                 elseif  i==29  L = mbins>=0.5  & mbins<=1.2;
%                 elseif  i==34  L = mbins>=0    & mbins<=0.4;
%                 elseif  i==36  L = mbins>=0.3  & mbins<=1.2;
%                 elseif  i==37  L = mbins>=0.3  & mbins<=1.2;
%                 elseif  i==38  L = mbins>=0.4  & mbins<=1.1;
%                 elseif  i==40  L = mbins>=0.5  & mbins<=1;
%                 elseif  i==42  L = mbins>=0.2  & mbins<=0.7;
%                 elseif  i==43  L = mbins>=0.2  & mbins<=1.5;
%                 elseif  i==44  L = mbins>=0.1  & mbins<=0.5;
%                 elseif  i==47  L = mbins>=0.4  & mbins<=0.75;
%                 elseif  i==49  L = mbins>=0    & mbins<=0.3;
%                 elseif  i==50  L = mbins>=0.5  & mbins<=1.1;
%                 elseif  i==51  L = mbins>=0.7  & mbins<=1.1;
%                 elseif  i==52  L = mbins>=0.4  & mbins<=1.1;
%                 elseif  i==56  L = mbins>=0.45 & mbins<=0.7;
%                 elseif  i==57  L = mbins>=0.1  & mbins<=0.65;
%                 elseif  i==58  L = mbins>=0    & mbins<=0.75;
%                 elseif  i==59  L = mbins>=0.6  & mbins<=1.1;
%                 elseif  i==62  L = mbins>=0.5  & mbins<=1;
%                 elseif  i==63  L = mbins>=0.5  & mbins<=1;
%                 elseif  i==64  L = mbins>=0    & mbins<=0.8;
%                 elseif  i==65  L = mbins>=0.4  & mbins<=0.8;
%                 elseif  i==67  L = mbins>=0.25 & mbins<=0.7;
%                 elseif  i==68  L = mbins>=0.25 & mbins<=0.7;
%                 elseif  i==69  L = mbins>=0.5  & mbins<=0.75;
%                 elseif  i==70  L = mbins>=0.75 & mbins<=1.25;
%                 elseif  i==71  L = mbins>=0.45 & mbins<=0.9;
%                 elseif  i==76  L = mbins>=0.55 & mbins<=0.9;
%                 elseif  i==77  L = mbins>=0.65 & mbins<=1.2;
%                 elseif  i==78  L = mbins>=0.2  & mbins<=0.6;
%                 elseif  i==79  L = mbins>=0.2  & mbins<=0.9;
%                 elseif  i==85  L = mbins>=0.2  & mbins<=0.5;
%                 elseif  i==86  L = mbins>=0.2  & mbins<=0.7;
%                 elseif  i==87  L = mbins>=0.65 & mbins<=1;
%                 elseif  i==89  L = mbins>=0.2  & mbins<=0.8;
%                 elseif  i==90  L = mbins>=0.4  & mbins<=1;
%                 elseif  i==91  L = mbins>=0.4  & mbins<=0.6;
%                 elseif  i==92  L = mbins>=0.4  & mbins<=0.6;
%                 elseif  i==93  L = mbins>=0.45 & mbins<=0.8;
%                 elseif  i==95  L = mbins>=0.2  & mbins<=0.8;
% 
%                 else
%                     L = mbins>=0 & mbins<=1.5;       % Zsa Zsa's time
%                 end
%             end
%             %
%             %         % fit data to high coh to estimate onset and offset time
%             %         if nanmax(n(:,7,i))>nanmax(n(:,6,i))
%             %             [bb, bsem, gof, ym] = sramp_fitW(mbins(L), m(L,7,i), sd(L,7,i), [], [], ...
%             %                 [min(mbins(L)) 0 min(mbins(L)); max(mbins(L)) 10 max(mbins(L))]);
%             %             Lt = mbins>=round(20*bb(1))/20 & mbins<=round(20*bb(3))/20;
%             %             if 1
%             %                 clf
%             %                 subplot(2,4,8)
%             %                 hold on
%             %                 plot(mbins, m(:,7,i), 'o', 'MarkerSize', 6, ...
%             %                     'MarkerFaceColor', lc(7,:), 'MarkerEdgeColor', 'none')
%             %                 plot(mbins(L),ym,'k', 'LineWidth', 3)
%             %                 hold off
%             %                 title(sprintf('%.2f - %.2f', bb(1), bb(3)))
%             %                 xlim([0 1.5]);
%             %                 ylim([0.4 1]);
%             %             end
%             %
%             %         else
%             %             [bb, bsem, gof, ym] = sramp_fitW(mbins(L), m(L,6,i), sd(L,6,i), [], [], ...
%             %                 [min(mbins(L)) 0 min(mbins(L)); max(mbins(L)) 10 max(mbins(L))]);
%             %             Lt = mbins>=round(20*bb(1))/20 & mbins<=round(20*bb(3))/20;
%             %             bb = nans(2,7);
%             %             if 1
%             %                 clf
%             %                 subplot(2,4,8)
%             %                 hold on
%             %                 plot(mbins, m(:,6,i), 'o', 'MarkerSize', 6, ...
%             %                     'MarkerFaceColor', lc(7,:), 'MarkerEdgeColor', 'none')
%             %                 plot(mbins(L),ym,'k', 'LineWidth', 3)
%             %                 hold off
%             %                 title(sprintf('%.2f - %.2f', bb(1), bb(3)))
%             %                 xlim([0 1.5]);
%             %                 ylim([0.4 1]);
%             %             end
%             %
%             %         end
% 
% 
% 
%             bb = nans(2,7);
%             for j = 1:7
%                 Lgd = ~isnan(m(:,j,i)) | ~isnan(sd(:,j,i));
%                 Lt  = L;
%                 if sum(Lgd&Lt)>3
%                     sd(sd(:,j,i)==0,j,i) = nanmean(sd(:,j,i));
%                     [bb(:,j), bbs] = regressW(m(Lt&Lgd,j,i), sd(Lt&Lgd,j,i), [ones(sum(Lt&Lgd),1) mbins(Lt&Lgd)]);
%                     b(j,i)    = bb(2,j);
%                     bs(j,i)   = bbs(2,2)-bb(2,j);
%                 end
%             end
% 
% 
%             if 0
%                 clf
%                 for j = 1:7
%                     subplot(2,4,j)
%                     hold on
%                     plot(mbins, m(:,j,i), 'o', 'MarkerSize', 6, ...
%                         'MarkerFaceColor', lc(j,:), 'MarkerEdgeColor', 'none')
%                     plot(mbins(L), bb(1,j)+bb(2,j)*mbins(L), '-', 'LineWidth', 3, 'Color', lc(j,:))
%                     hold off
%                     xlim([0 1.5]);
%                     ylim([0.4 1]);
% 
%                     title(sprintf('%.2f',coh(j)))
%                 end
%                 pause
%             end
%         end
%     end
% 
%     for j = 1:7
%         L   = ~isnan(b(j,:));
%         
%     %    [fit fiti r ri stat] = regress(b(j,L)',[ones(sum(L),1) ses(L)]);
%         [fit fiti h p] = nestedFW(b(j,L)',bs(j,L)',[ones(sum(L),1) ses(L)]);
%     
%         f(j)  = fit(2);
%         fi(j) = fiti(2,2)-fit(2);
%         s(j)  = p;
%     end
% 
% 
%     save(savepath, 'f', 'fi', 's')
% 
% else
%     load(savepath)
% end
% 
% 
% if 0
% 
%     figure
%     hold on
%     plot([1 3.2 6.4 12.8 25.6 51.2 99.9],f,'o','MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
%     plot([[1 3.2 6.4 12.8 25.6 51.2 99.9];[1 3.2 6.4 12.8 25.6 51.2 99.9]],[f-fi; f+fi],'k-')
%     hold off
%     set(gca, 'xscale', 'log', 'xlim', [1 99.9], ...
%         'XTick', [1 3.2 6.4 12.8 25.6 51.2 99.9], 'XTickLabel', [0 3.2 6.4 12.8 25.6 51.2 99.9])
% 
% end
