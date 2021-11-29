function [fits, sems] = getML_ddOUn(fname, recompute)
% get coh and time dependence of neurometric performance using OU model
% INPUTS:
%   fname - 'txt' filename

if ~isempty(findstr(fname, 'LIP'))
    dtype = 'LIP';

elseif ~isempty(findstr(fname, 'MT'))
    if ~isempty(findstr(fname, 'PRe'))
        dtype = 'MTPRe';
    else
        dtype = 'MT';
    end

else
    return;
end


savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_ddOUn_' fname(1:2) dtype '.mat'];



if recompute
    if strcmp(dtype, 'MT') | strcmp(dtype, 'MTPRe')
        d1name = 'ddir';
        ecname = 'dot_dir';

        % get cumulative spike rate for MT
        bb   = 0;
        be   = 1000;
        bw   = 25;
        bs   = 25;
        bins = [zeros(size([bb+bw:bs:be]')) [bb+bw:bs:be]'];
        mbins = (bins(:,2)-bw/2)/1000;
        crtf = [0 1];
        nvf  = [0 1];
        dirf = 0;
        [rM, rSD, rN, rMn] = getML_mRate(fname, bins, crtf, nvf, dirf, 0);
        CohD               = getML_SigNoise(fname,[0 1],[0 1],[0],0);  % use to estimate the bounds for drift rate
        CohD(2,:)          = [];


        % get average response before dots on
        a      = getML_txt(fname);
        fn     = a.data{strcmp(a.name,'dat_fn')};
        usable = a.data{strcmp(a.name,'usable')};
        uid    = a.data{strcmp(a.name,'uid')};
        d1     = a.data{strcmp(a.name,d1name)};
        tkflg  = a.data{strcmp(a.name,'taskflag')};

        coh         = [0 3.2 6.4 12.8 25.6 51.2 99.9]';
        m           = rM-rM([8:14 1:7], :, :);    % select between normalized and unnormalized spikes
        %sd          = rSD./(rN.^0.5);
        sd          = rSD./(rN.^0.5)+rSD([8:14 1:7], :, :)./(rN([8:14 1:7], :, :).^0.5);    % select between normalized and unnormalized spikes

        m(8:end,:,:)  = [];
        sd(8:end,:,:) = [];

        fits = nans(2,length(fn));
        sems = nans(2,length(fn));

        global FIRA
        for i = 1:length(fn)
            if usable(i) == 1
                fprintf('%d: %s\n', i, fn{i})

                % prepare data
                r  = m(:,:,i);
                rd = sd(:,:,i);
                t  = mbins;
                t  = repmat(t',length(coh),1);
                t  = t(:);
                c  = repmat(coh/100,1,length(bins));
                c  = c(:);
                data = [c t r(:) rd(:)];
                %data = [c t r(:)];


                % estimate initial conditions
                % first estimate A term, then use 75% and 125% of that value
                % as the lower and the upper bounds
                k0        = [CohD(i);       -0.00001];
                kl        = [0.75*CohD(i);   -20     ];
                ku        = [1.25*CohD(i);   -0.00001];
                tau       = 0.1;
                [fits(:,i),f,e,o,l,g,H] = fmincon('ddOUn_err', k0, [], [], [], [], ...
                    kl, ku, [], ...
                    optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'), ...
                    data, tau);
                sems(:,i) = sqrt(diag(-((-H)^(-1))));

                % plot
                if 0
                    lc  = [];
                    for j = 1:7
                        lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
                    end

                    clf
                    %                    subplot(1,2,1)
                    hold on
                    for j = 1:length(coh)
                        plot(mbins, m(j,:,i).*mbins', '.', 'LineWidth', 2, 'MarkerEdgeColor', lc(j,:))
                        plot(mbins, ddOUn_val(fits(:,i),[repmat(coh(j)/100,length(mbins),1) mbins], tau), '-', 'LineWidth', 2, 'Color', lc(j,:))
                    end
                    hold off
                    title(sprintf('%.2f, %.2f, %.2f', fits(1,i), fits(2,i)))

                    %                     subplot(1,2,2)
                    %                     plot(coh/100, r(end-7+1:end), '.k')
                    %                     lsline

                    pause
                end
            end
        end



    elseif strcmp(dtype, 'LIP')
        d1name = 'trg_dir';
        ecname = 'trg1_dir';

        % get LIP spike rate
        bb     = 0;
        be     = 1500;
        bs     = 25;
        bw     = 50;
        bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
        mbins  = mean(bins,2)/1000;
        crtf = [1];
        nvf  = [0 1];
        dirf = 0;
        [rM, rSD, rN, rMn] = getML_mRate(fname, bins, crtf, nvf, dirf, 0);
        on                 = getML_LIPModelroc(fname(1:2), 0);
        on                 = 1000*on(3,:);
        CohD               = getML_LIPModel(fname(1:2), 0, on);    % use to estimate the bounds for drift rate
        CohD               = 1000*CohD(3,:);


        % get average response before dots on
        a      = getML_txt(fname);
        fn     = a.data{strcmp(a.name, 'dat_fn')};
        usable = a.data{strcmp(a.name, 'usable')};
        uid    = a.data{strcmp(a.name, 'uid')};
        d1     = a.data{strcmp(a.name, d1name)};
        tkflg  = a.data{strcmp(a.name, 'taskflag')};

        coh         = [0 3.2 6.4 12.8 25.6 51.2 99.9]'/100;
        Lcoh        = logical([1 1 1 1 1 1 1]);
        m           = rM-rM([8:14 1:7], :, :);    % select between normalized and unnormalized spikes
        %sd          = rSD./(rN.^0.5);
        sd          = rSD./(rN.^0.5)+rSD([8:14 1:7], :, :)./(rN([8:14 1:7], :, :).^0.5);    % select between normalized and unnormalized spikes

        m(8:end,:,:)  = [];
        sd(8:end,:,:) = [];

        fits = nans(3,length(fn));
        sems = nans(3,length(fn));

        global FIRA
        for i = 116:-1:1%length(fn):-1:1
            if usable(i) == 1
                fprintf('%d: %s\n', i, fn{i})


                if strcmp(fname(1:2), 'ZZ')
                    if      i==95    L = mbins>=0.15  & mbins<=0.35;
                    elseif  i==93    L = mbins>=0.55  & mbins<=0.9;
                    elseif  i==92    L = mbins>=0.35  & mbins<=0.75;
                    elseif  i==91    L = mbins>=0.35  & mbins<=0.6;
                    elseif  i==90    L = mbins>=0.4   & mbins<=0.55;
                    elseif  i==89    L = mbins>=0.5   & mbins<=0.8;
                    elseif  i==88    L = mbins>=0.4   & mbins<=0.6;
                    elseif  i==87    L = mbins>=0.4   & mbins<=0.625;
                    elseif  i==86    L = mbins>=0.25  & mbins<=0.4;
                    elseif  i==85    L = mbins>=0.3   & mbins<=0.45;
                    elseif  i==82    L = mbins>=0.3   & mbins<=0.5;
                    elseif  i==81    L = mbins>=0.25  & mbins<=0.375;
                    elseif  i==79    L = mbins>=0.6   & mbins<=0.75;
                    elseif  i==78    L = mbins>=0.2   & mbins<=0.8;
                    elseif  i==77    L = mbins>=0.4   & mbins<=0.6;
                    elseif  i==76    L = mbins>=0.55  & mbins<=0.9;
                    elseif  i==74    L = mbins>=0.2   & mbins<=0.45;
                    elseif  i==72    L = mbins>=0.35  & mbins<=0.55;
                    elseif  i==71    L = mbins>=0.3   & mbins<=0.8;
                    elseif  i==70    L = mbins>=0.35  & mbins<=0.8;
                    elseif  i==69    L = mbins>=0.15  & mbins<=0.7;
                    elseif  i==68    L = mbins>=0.5   & mbins<=0.7;
                    elseif  i==67    L = mbins>=0.35  & mbins<=0.65;
                    elseif  i==66    L = mbins>=0.2   & mbins<=0.7;
                    elseif  i==65    L = mbins>=0.45  & mbins<=0.65;
                    elseif  i==64    L = mbins>=0.4   & mbins<=1.5;
                    elseif  i==63    L = mbins>=0.6   & mbins<=0.95;
                    elseif  i==62    L = mbins>=0.375 & mbins<=0.55;
                    elseif  i==59    L = mbins>=0.1   & mbins<=0.35;
                    elseif  i==58    L = mbins>=0.5   & mbins<=0.8;
                    elseif  i==57    L = mbins>=0.675 & mbins<=0.9;
                    elseif  i==56    L = mbins>=0.5   & mbins<=0.675;
                    elseif  i==53    L = mbins>=0.2   & mbins<=0.6;
                    elseif  i==52    L = mbins>=0.7   & mbins<=1.1;
                    elseif  i==51    L = mbins>=0.2   & mbins<=0.8;
                    elseif  i==50    L = mbins>=0.4   & mbins<=0.6;
                    elseif  i==49    L = mbins>=0.7   & mbins<=1;
                    elseif  i==47    L = mbins>=0.35  & mbins<=0.7;
                    elseif  i==45    L = mbins>=0.6   & mbins<=1;
                    elseif  i==44    L = mbins>=0.3   & mbins<=0.5;
                    elseif  i==43    L = mbins>=0.7   & mbins<=1.2;
                    elseif  i==42    L = mbins>=0.3   & mbins<=0.8;
                    elseif  i==40    L = mbins>=0.4   & mbins<=0.95;
                    elseif  i==39    L = mbins>=0.45  & mbins<=1;
                    elseif  i==38    L = mbins>=0.4   & mbins<=0.9;
                    elseif  i==37    L = mbins>=0.4   & mbins<=1.2;
                    elseif  i==36    L = mbins>=0.5   & mbins<=1.2;
                    elseif  i==35    L = mbins>=0.4   & mbins<=1.2;
                    elseif  i==34    L = mbins>=0.7   & mbins<=1.3;
                    elseif  i==33    L = mbins>=0.1   & mbins<=1.15;
                    elseif  i==32    L = mbins>=0.1   & mbins<=1.2;
                    elseif  i==31    L = mbins>=0.1   & mbins<=1.3;
                    elseif  i==30    L = mbins>=0.2   & mbins<=0.7;
                    elseif  i==29    L = mbins>=0.5   & mbins<=1;
                    elseif  i==28    L = mbins>=0.2   & mbins<=1.2;
                    elseif  i==27    L = mbins>=0.4   & mbins<=0.9;
                    elseif  i==26    L = mbins>=0.3   & mbins<=0.7;
                    elseif  i==25    L = mbins>=0.2   & mbins<=1;
                    elseif  i==24    L = mbins>=0     & mbins<=1.3;
                    elseif  i==23    L = mbins>=0.3   & mbins<=1;
                    elseif  i==22    L = mbins>=0.5  & mbins<=1.5;
                    elseif  i==20    L = mbins>=0.3   & mbins<=1;
                    elseif  i==18    L = mbins>=0     & mbins<=0.7;
                    elseif  i==17    L = mbins>=0.2   & mbins<=0.8;
                    elseif  i==16    L = mbins>=0.75  & mbins<=1.5;
                    elseif  i==15    L = mbins>=0.2   & mbins<=1.1;
                    elseif  i==14    L = mbins>=0.2   & mbins<=1.3;
                    elseif  i==13    L = mbins>=0.5   & mbins<=1.4;
                    elseif  i==12    L = mbins>=0.1   & mbins<=1.5;
                    elseif  i==11    L = mbins>=0.2   & mbins<=1.5;
                    elseif  i==10    L = mbins>=0.1   & mbins<=1.5;
                    elseif  i==9     L = mbins>=0     & mbins<=1;
                    elseif  i==8     L = mbins>=0     & mbins<=1;
                    elseif  i==7     L = mbins>=0     & mbins<=1;
                    elseif  i==6     L = mbins>=0     & mbins<=1;
                    elseif  i==5     L = mbins>=0     & mbins<=1.5;
                    elseif  i==4     L = mbins>=0     & mbins<=1;
                    elseif  i==3     L = mbins>=0     & mbins<=0.8;
                    elseif  i==2     L = mbins>=0     & mbins<=1;
                    elseif  i==1     L = mbins>=0     & mbins<=1.5;

                    else
                        L = mbins>=0.4 & mbins<=0.6;  % Zsa Zsa's time
                    end
                end

                if strcmp(fname(1:2), 'Cy')   % pick the best time period to fit data
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
                    elseif i==97    L = mbins>=0.25    & mbins<=0.4;
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



                % prepare data
                t     = mbins(L);
                t     = repmat(t',sum(Lcoh),1);
                t     = t(:);
                c     = repmat(coh(Lcoh),1,sum(L));
                c     = c(:);
                r     = m(Lcoh,L,i)/5;
                rd    = sd(Lcoh,L,i);
                r     = r(:);
                rd    = rd(:);
                data  = [c t r(:) rd(:)];


                % estimate initial conditions
                % first estimate A term, then use 75% and 125% of that value
                % as the lower and the upper bounds
                k0        = [CohD(i)/5;     -0.00001;  1];
                kl        = [0.75*CohD(i)/5; -20;      0];
                ku        = [1.25*CohD(i)/5; -0.00001; 20];
                tau       = t(1);
                bnd       = myprctile(r,99);
                param     = [tau bnd];
                [fits(:,i),f,e,o,l,g,H] = fmincon('ddOUn_err', k0, [], [], [], [], ...
                    kl, ku, [], ...
                    optimset('LargeScale', 'off', 'Display', 'off', 'Diagnostics', 'off'), ...
                    data, param, @ddOUnb_val);

                %                 sems(:,i) = sqrt(diag(-((-H)^(-1))));

                % plot
                if 0
                    lc  = [];
                    for j = 1:7
                        lc = [lc; [0.9 0.9 0.9]-(j-1)/7];
                    end

                    clf
                    hold on
                    for j = 1:length(coh)
                        plot(mbins, m(j,:,i), '.', 'LineWidth', 2, 'MarkerEdgeColor', lc(j,:))
                        plot(mbins(L), 5*ddOUnb_val(fits(:,i),[repmat(coh(j),sum(L),1) mbins(L)], [tau,bnd]), '-', 'LineWidth', 2, 'Color', lc(j,:))
                    end
                    hold off
                    title(sprintf('%.2f, %.2f, %.2f', fits(1,i), fits(2,i)))
                    set(gca, 'xlim', [0.1 0.5])

                    pause
                end
            end
        end
    end


    save(savepath, 'fits', 'sems')
else

    load(savepath)
end
