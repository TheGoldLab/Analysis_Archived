function [a, ae] = getML_LIPModelroc99(Monk, recompute)
% fit LIP roc areas at 99.9% coh to a piecewise-linear model

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_LIPModelroc99_' Monk '.mat'];

if recompute
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
    bs     = 50;
    bw     = 100;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    mbins  = mean(bins,2)/1000;

    [rs rn rsd] = getML_ROCs([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, [], 0);

    a   = nans(3,length(fn));
    ae  = nans(3,2,length(fn));

    m   = rs;               % select between normalized and unnormalized spikes
    sd  = rsd./(rn.^0.5); % standard error

    for i = 1:length(fn)
        if usable(i)==1
            fprintf(sprintf('%d: %s\n', i, fn{i}))

            if strcmp(Monk, 'Cy')  % cyrus's time
                bl = 0.5;

                if      i==1    L = mbins>=0   & mbins<=0.8;
                elseif  i==2    L = mbins>=0.2 & mbins<=1.4;
                elseif  i==3    L = mbins>=1   & mbins<=1.25;
                elseif  i==4    L = mbins>=0.9 & mbins<=1.3;
                elseif  i==5    L = mbins>=0.3 & mbins<=0.75;
                elseif  i==6    L = mbins>=0   & mbins<=0.6;
                elseif  i==7    L = mbins>=0.2 & mbins<=0.8;
                elseif  i==8    L = mbins>=0.1 & mbins<=0.55;
                elseif  i==12   L = mbins>=0.1 & mbins<=0.55;
                elseif  i==13   L = mbins>=0   & mbins<=0.3;
                elseif  i==14   L = mbins>=0   & mbins<=0.6;
                elseif  i==16   L = mbins>=0   & mbins<=0.6;
                elseif  i==17   L = mbins>=0.2 & mbins<=1.4;
                elseif  i==18   L = mbins>=0.4 & mbins<=0.8;
                elseif  i==19   L = mbins>=0.1 & mbins<=0.75;
                elseif  i==21   L = mbins>=0.1 & mbins<=1.5;
                elseif  i==22   L = mbins>=0.1 & mbins<=0.5;
                elseif  i==23   L = mbins>=0.1 & mbins<=0.5;
                elseif  i==24   L = mbins>=0   & mbins<=0.5;
                elseif  i==25   L = mbins>=0   & mbins<=0.5;
                elseif  i==28   L = mbins>=0   & mbins<=0.5;
                elseif  i==30   L = mbins>=0.1 & mbins<=0.8;
                elseif  i==31   L = mbins>=0   & mbins<=0.5;
                elseif  i==32   L = mbins>=0.3 & mbins<=0.7;
                elseif  i==33   L = mbins>=0   & mbins<=0.4;
                elseif  i==34   L = mbins>=0   & mbins<=0.6;
                elseif  i==35   L = mbins>=0   & mbins<=0.5;
                elseif  i==37   L = mbins>=0.2 & mbins<=1.2;
                elseif  i==38   L = mbins>=0.1 & mbins<=0.7;
                elseif  i==40   L = mbins>=0.1 & mbins<=0.5;
                elseif  i==41   L = mbins>=0.1 & mbins<=1.1;
                elseif  i==42   L = mbins>=0   & mbins<=1;
                elseif  i==44   L = mbins>=0.1 & mbins<=1;
                elseif  i==45   L = mbins>=0   & mbins<=0.8;
                elseif  i==46   L = mbins>=0   & mbins<=0.3;
                elseif  i==47   L = mbins>=0   & mbins<=0.9;
                elseif  i==48   L = mbins>=0   & mbins<=0.5;
                elseif  i==49   L = mbins>=0   & mbins<=0.4;
                elseif  i==51   L = mbins>=0.2 & mbins<=0.8;
                elseif  i==52   L = mbins>=0   & mbins<=0.4;
                elseif  i==53   L = mbins>=0.3 & mbins<=0.9;
                elseif  i==55   L = mbins>=0.1 & mbins<=0.9;
                elseif  i==56   L = mbins>=0   & mbins<=0.7;
                elseif  i==57   L = mbins>=0   & mbins<=0.7;
                elseif  i==60   L = mbins>=0.4 & mbins<=1;
                elseif  i==61   L = mbins>=0   & mbins<=0.5;
                elseif  i==65   L = mbins>=0.55 & mbins<=0.75;
                elseif  i==66   L = mbins>=0   & mbins<=0.5;
                elseif  i==74   L = mbins>=0   & mbins<=0.5;
                elseif  i==75   L = mbins>=0   & mbins<=0.5;
                elseif  i==77   L = mbins>=0   & mbins<=1;
                elseif  i==78   L = mbins>=0.1   & mbins<=0.95;
                elseif  i==81   L = mbins>=0   & mbins<=0.65;
                elseif  i==83   L = mbins>=0   & mbins<=0.35;
                elseif  i==84   L = mbins>=0   & mbins<=0.35;
                elseif  i==85   L = mbins>=0   & mbins<=0.5;
                elseif  i==86   L = mbins>=0   & mbins<=0.4;
                elseif  i==87   L = mbins>=0   & mbins<=0.4;
                elseif  i==90   L = mbins>=0   & mbins<=0.9;
                elseif  i==91   L = mbins>=0.05& mbins<=0.4;
                elseif  i==94   L = mbins>=0   & mbins<=0.4;
                elseif  i==98   L = mbins>=0   & mbins<=0.4;
                elseif  i==99   L = mbins>=0   & mbins<=0.4;
                elseif  i==100  L = mbins>=0   & mbins<=0.3;
                elseif  i==102  L = mbins>=0   & mbins<=0.3;
                elseif  i==104  L = mbins>=0.5 & mbins<=0.9;
                elseif  i==107  L = mbins>=0.5 & mbins<=0.8;
                elseif  i==108  L = mbins>=0.1 & mbins<=0.4;
                elseif  i==109  L = mbins>=0.1 & mbins<=0.3;
                elseif  i==110  L = mbins>=0   & mbins<=0.3;
                elseif  i==112  L = mbins>=0   & mbins<=0.3;
                elseif  i==114  L = mbins>=0   & mbins<=0.3;
                elseif  i==116  L = mbins>=0   & mbins<=0.3;
                elseif  i==119  L = mbins>=0.1 & mbins<=0.55;
                elseif  i==122  L = mbins>=0   & mbins<=0.6;
                elseif  i==123  L = mbins>=0   & mbins<=0.8;
                elseif  i==138  L = mbins>=0   & mbins<=0.3;

                else
                    L = mbins>=0 & mbins<=0.6;
                end

            else
                bl = 0.5;

                if      i==1   L = mbins>=0.4  & mbins<=1;
                elseif  i==5   L = mbins>=0    & mbins<=0.8;
                elseif  i==6   L = mbins>=0    & mbins<=0.7;
                elseif  i==8   L = mbins>=0    & mbins<=0.7;
                elseif  i==12  L = mbins>=0.6  & mbins<=1.1;
                elseif  i==13  L = mbins>=0.4  & mbins<=1;
                elseif  i==15  L = mbins>=0.2  & mbins<=1.3;
                elseif  i==18  L = mbins>=0.1  & mbins<=0.8;
                elseif  i==20  L = mbins>=0.2  & mbins<=1;
                elseif  i==22  L = mbins>=0    & mbins<=1.4;
                elseif  i==23  L = mbins>=0.2  & mbins<=0.65;
                elseif  i==24  L = mbins>=0.6  & mbins<=1.4;
                elseif  i==25  L = mbins>=0    & mbins<=0.5;
                elseif  i==26  L = mbins>=0    & mbins<=1.2;
                elseif  i==27  L = mbins>=0.1  & mbins<=1.2;
                elseif  i==28  L = mbins>=0.6  & mbins<=1.2;
                elseif  i==29  L = mbins>=0.5  & mbins<=1.4;
                elseif  i==30  L = mbins>=0.4  & mbins<=0.7;
                elseif  i==31  L = mbins>=0.4  & mbins<=1.25;
                elseif  i==32  L = mbins>=0    & mbins<=1.4;
                elseif  i==33  L = mbins>=0.4  & mbins<=1.4;
                elseif  i==34  L = mbins>=0.1  & mbins<=0.4;
                elseif  i==35  L = mbins>=0.2  & mbins<=1.2;
                elseif  i==36  L = mbins>=0.8  & mbins<=1;
                elseif  i==37  L = mbins>=0.3  & mbins<=1.1;
                elseif  i==38  L = mbins>=0.35 & mbins<=1;
                elseif  i==39  L = mbins>=0.3  & mbins<=1.2;
                elseif  i==40  L = mbins>=0.3  & mbins<=1;
                elseif  i==42  L = mbins>=0.3  & mbins<=0.6;
                elseif  i==43  L = mbins>=0.3  & mbins<=1.3;
                elseif  i==44  L = mbins>=0.25  & mbins<=0.5;
                elseif  i==47  L = mbins>=0.2  & mbins<=0.7;
                elseif  i==48  L = mbins>=0    & mbins<=0.5;
                elseif  i==49  L = mbins>=0    & mbins<=0.35; bl = 0.55;
                elseif  i==50  L = mbins>=0.7  & mbins<=1;
                elseif  i==52  L = mbins>=0.4  & mbins<=1.05;
                elseif  i==53  L = mbins>=0.2  & mbins<=0.6;
                elseif  i==56  L = mbins>=0.35 & mbins<=0.65;
                elseif  i==57  L = mbins>=0.2  & mbins<=0.65;
                elseif  i==58  L = mbins>=0.3  & mbins<=0.65;
                elseif  i==59  L = mbins>=0.7  & mbins<=0.9;
                elseif  i==62  L = mbins>=0.7  & mbins<=0.8;
                elseif  i==65  L = mbins>=0    & mbins<=0.3;
                elseif  i==66  L = mbins>=0.4  & mbins<=1.3;
                elseif  i==70  L = mbins>=0.3  & mbins<=0.95;
                elseif  i==72  L = mbins>=0.3  & mbins<=0.6;
                elseif  i==74  L = mbins>=0.5  & mbins<=0.7;
                elseif  i==79  L = mbins>=0.2  & mbins<=0.9;
                elseif  i==81  L = mbins>=0.4  & mbins<=0.75;
                elseif  i==82  L = mbins>=0.95 & mbins<=1.1;
                elseif  i==85  L = mbins>=0.3  & mbins<=0.45;
                elseif  i==86  L = mbins>=0.2  & mbins<=0.4;
                elseif  i==87  L = mbins>=0.45 & mbins<=0.5;
                elseif  i==88  L = mbins>=0.2  & mbins<=0.9;
                elseif  i==89  L = mbins>=0.6  & mbins<=0.8;
                elseif  i==90  L = mbins>=0.2  & mbins<=0.6;
                elseif  i==91  L = mbins>=0.3  & mbins<=0.6;
                elseif  i==93  L = mbins>=0.5  & mbins<=0.7;
                elseif  i==95  L = mbins>=0.4  & mbins<=0.5;

                else
                    L = mbins>=0.2 & mbins<=0.8;       % Zsa Zsa's time
                end
            end


            bcon = [min(mbins(L)) 0 0; max(mbins(L)) 6 1];


            mbP   = mbins(L);
            mbP   = mbP(:);
            c     = repmat(99.9,1,sum(L));
            c     = c(:);
            r     = m(L,7,i);
            rd    = sd(L,7,i);
            r     = r(:);
            rd    = rd(:);
            rd(rd==0) = nanmean(rd(rd~=0));
            
            
            if ~all(isnan(r))
                % fit data to piecewise-linear model to get plateau
                [a(:,i) xx gof ym] = sramp_fitW2(mbP, r, rd, bl, [], bcon);
              
                % estimate se
                % Check for input etype, which indicates that we will use
                %    Monte Carlo resampling (parametric bootstrap).
                %   See Wichmann & Hill (2001)
                %       The psychometric function: II. Bootstrap-
                %           based confidence intervals and sampling
                % citype{1} is number of simulated data sets to use
                mcn = 100;
                CI  = 68;
                mcfits = zeros(size(a,1),mcn);
                for ii = 1:mcn
                    if ~mod(ii,10)
                        disp(sprintf('Bootstrap CIs, set %d', ii))
                    end
                    % compute fit on simulated data set
                    ysim = normrnd(r,rd);
                    [mcfits(:,ii) xx gof mcm] = sramp_fitW2(mbP, ysim, rd, bl, [], bcon);
                end
                ae(:,:,i) = [myprctile(mcfits',50-CI/2)' myprctile(mcfits',50+CI/2)'];
                a(:,i)
                ae(:,:,i)
                % plot and check fit
                if 0
                    clf
                    hold on
                    plot(mbins, m(:,7,i), 'o', 'MarkerSize', 6, ...
                        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
                    plot(mbP, ym, '-', 'LineWidth', 3, 'Color', 'k')
                    xlim([0 1.5]);
                    ylim([0.4 1]);
                    title(sprintf('%.2f (%.2f), %.2f (%.2f)', a(2,i), ae(2,i), a(3,i), ae(3,i)))
                    pause
                end
            end

        end
    end


    save(savepath, 'a', 'ae')
else

    load(savepath)
end

