function [fits sems] = getML_ddOU(Monk, recompute)
% get coh and time dependence of performance using OU model
savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_ddOU_' Monk '.mat'];


if recompute
    % get session data
    fn     = [Monk 'TRain_psy.txt'];
    a      = getML_txt(fn);
    fname  = a.data{strcmp(a.name,'dat_fn')};
    bb     = 0;
    be     = 1500;
    bw     = 300;
    bs     = 150;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    mbins  = mean(bins,2);
    ftime  = [0:1:1500]'/1000;
    coh    = [0 3.2 6.4 12.8 25.6 51.2 99.9]';

    % get performance
    praw   = nans(7,length(mbins),length(fname));
    nraw   = nans(7,length(mbins),length(fname));
    pou    = nans(7,length(ftime));
    pdE    = nans(7,length(ftime));
    fits   = nans(2,length(fname));
    sems   = nans(2,length(fname));
load(savepath)

    global FIRA
    for i = 1:length(fname)
        fprintf('%d: %s\n', i, fname{i})
        openFIRA(fname{i})

        %%%                                                 %%%
        %       get raw performance for each coh and time     %
        %%%                                                 %%%
        warning off
        Lgd     = ~isnan(FIRA.ecodes.data(:,1));
        bFIRA   = [FIRA.ecodes.data(Lgd,getFIRA_ecodeColumnByName('dot_dir'))...
            FIRA.ecodes.data(Lgd,getFIRA_ecodeColumnByName('dot_coh'))...
            FIRA.ecodes.data(Lgd,getFIRA_ecodeColumnByName('dot_off')) - FIRA.ecodes.data(Lgd,getFIRA_ecodeColumnByName('dot_on'))...
            FIRA.ecodes.data(Lgd,getFIRA_ecodeColumnByName('correct'))];
        for j = 1:length(mbins)
            Lt = bFIRA(:,3)>=bins(j,1) & bFIRA(:,3)<bins(j,2);
            for k = 1:length(coh)
                Lcoh = bFIRA(:,2)==coh(k);
                nraw(k,j,i) = sum(bFIRA(Lt&Lcoh,4)==1 | bFIRA(Lt&Lcoh,4)==0);
                praw(k,j,i) = sum(bFIRA(Lt&Lcoh,4)==1)/nraw(k,j,i);
            end
        end
        warning on




        % get A
        d = [getFIRA_ecodesByName('dot_coh')./100 ...
            (getFIRA_ecodesByName('dot_off')-getFIRA_ecodesByName('dot_on'))./1000 ...
            sign(cos(pi/180.*getFIRA_ecodesByName('dot_dir'))) ...
            getFIRA_ecodesByName('correct')];

        Lgood = d(:,4)>=0;
        d = d(Lgood,:);

        if i==1
            init = 0;
        else
            init = fits(1,i-1);
        end
        
        [f,s] = ctPsych_fit(@ddOU1, d(:,1:2), d(:,4), [], [], init);
        for j = 1:7
            poua(j,:) = ddOU1(f,[repmat(coh(j)/100,size(ftime)) ftime ones(size(ftime))]);
        end
        fits(1,i) = f(1);
        sems(1,i) = s(1);


        
        % get time dependence
        d = [getFIRA_ecodesByName('dot_coh')./100 ...
            (getFIRA_ecodesByName('dot_off')-getFIRA_ecodesByName('dot_on'))./1000 ...
            sign(cos(pi/180.*getFIRA_ecodesByName('dot_dir'))) ...
            getFIRA_ecodesByName('correct')];

        % select trials
        if findstr(fn,'Cy')
            if i==16 | i==18 | i==33 | i==38 | i==62 | i==97
                Lgood = d(:,4)>=0 & d(:,1)<0.6;
            elseif i==22 | i==23
                Lgood = d(:,4)>=0;
            elseif i==106
                Lgood = d(:,4)>=0 & d(:,1)<0.6;
            elseif i==109 | i==131 | i==133 | i==139 | i==141 | i==150
                Lgood = d(:,4)>=0 & d(:,1)<0.2;
            elseif i==110
                Lgood = d(:,4)>=0 & d(:,1)<0.5;
            elseif i==113
                Lgood = d(:,4)>=0 & d(:,1)<0.5 & d(:,2)<0.8;
            elseif i==134 | i==135 | i==138
                Lgood = d(:,4)>=0 & d(:,1)<0.5;
            elseif i>=99 & i<=127 | i==149 | i==151 | i==154
                Lgood = d(:,4)>=0 & d(:,1)<0.6 & d(:,2)<0.8;
            elseif i>=127
                Lgood = d(:,4)>=0 & d(:,1)<0.5 & d(:,2)<0.6;
            else
                Lgood = d(:,4)>=0 & d(:,1)<0.5;
            end
        else
            if    i==16 | i==44
                Lgood = d(:,4)>=0 & d(:,1)<0.6 & d(:,2)<0.8;
            elseif i==106 | i==111 | i==123 | i==125 | i==128 | i==130
                Lgood = d(:,4)>=0 & d(:,1)<0.5;
            elseif i==109
                Lgood = d(:,4)>=0 & d(:,1)<0.2;
            else
                Lgood = d(:,4)>=0 & d(:,1)<0.6;
            end
        end
        d = d(Lgood,:);


        % note:
        %  For Zsa Zsa - use only 99% coh doesn't work (no trend of learning)
        %              - only 51.2% coh doesn't work (better trend, but still
        %              noisy
        %              - use only short viewing time doesn't work, need to play
        %              with time more


        init = [0 -20]';
        [f,s] = ctPsych_fit(@ddOU2, d(:,1:2), d(:,4), [], [], init);
        for j = 1:7
            pou(j,:) = ddOU2(f,[repmat(coh(j)/100,size(ftime)) ftime ones(size(ftime))]);
        end
        fits(2,i) = f(2);
        sems(2,i) = s(2);
        
      
        % plot
        if 0
            lc  = [];
            for j = 1:7
                lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
            end

            cla
            hold on
            for j = 1:7
                plot(mbins/1000,praw(j,:,i),'.', 'MarkerEdgeColor', lc(j,:))
                plot(ftime,pou(j,:),'-', 'Color', lc(j,:), 'LineWidth', 2)
                set(gca, 'xlim', [0 1], 'ylim', [0.35 1])
            %    title(sprintf('dou - %.2f',fits(2,i)))
            end
            hold off
            
            pause
        end
    end

    save(savepath, 'fits', 'sems')
else
    load(savepath)
end


% figure
% subplot(2,2,1)
% hold on
% plot(foua(1,:),'k.')
% plot(nanrunmean(foua(1,:),5),'-r', 'LineWidth', 2)
% hold off
% title('OU-A')
% xlim([1 160])
%
% subplot(2,2,2)
% hold on
% plot(fou(2,:),'k.')
% plot(nanrunmean(fou(2,:),5),'-r', 'LineWidth', 2)
% hold off
% title('OU-tau')
% xlim([1 160])
%
% subplot(2,2,3)
% scatter(foua(1,:),fou(2,:),'.k')
% title(num2str(corr(foua(1,:)',fou(2,:)'),'%.2f'))


