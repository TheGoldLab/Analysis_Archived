function [m v n R] = getML_LIPfano(Monk,recompute)
% getML_LIPFano
% test if LIP fano factor is close to 1. If LIP responses is bimodal, the
% variance will be much larger than the mean.

% get the mean/variance during the last 200ms of the ramp for each coh, for
% each cell

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_LIPfano_' Monk '.mat'];


if recompute

    % load data
    a      = getML_txt([Monk 'TRain_LIP.txt']);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    ses    = a.data{strcmp(a.name,'session')};
    usable = a.data{strcmp(a.name,'usable')};
    d1     = a.data{strcmp(a.name,'trg_dir')};
    uid    = a.data{strcmp(a.name,'uid')};

    v      = nans(14,length(fn));
    m      = nans(14,length(fn));
    n      = nans(14,length(fn));

    R      = [];
    
    % get var/mean
    global FIRA
    for i = 1:length(fn)
        if usable(i)==1
            fprintf(sprintf('%d: %s\n', i, fn{i}))
            openFIRA(fn{i})

            % select for cohs
            [Lcoh, Ucoh] = selectFIRA_trialsByUniqueID('dot_coh');
            Lcoh = Lcoh(:, ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9])); % remove strange coh other than the standard cohs
            Ucoh = Ucoh(ismember(Ucoh, [0 3.2 6.4 12.8 25.6 51.2 99.9]));    % in one pre-training file there's a 9% coh condition

            % select for dir
            [Ldir, Udir] = selectFIRA_trialsByUniqueID('trg1_dir');
            Idir         = [find(round(Udir)==d1(i)) find(round(Udir)==mod(d1(i)+180,360))];

            % select for task
            Ltk = ismember(getFIRA_ecodesByName('task','id'), [3,6]);

            % select for correct/incorrect
            Lcrt =  ismember(getFIRA_ecodesByName('correct', 'id'), 1);

            % get selection array for each condition
            Lr = zeros(length(Ldir), 14);
            c = [0 3.2 6.4 12.8 25.6 51.2 99.9];
            for j = 1:7
                if ismember(c(j),Ucoh)
                    Lr(:,j)   = Ldir(:,Idir(1)) & Lcoh(:,Ucoh==c(j)) & Ltk & Lcrt;
                    Lr(:,j+7) = Ldir(:,Idir(2)) & Lcoh(:,Ucoh==c(j)) & Ltk & Lcrt;
                end
            end

            bb     = -200;
            be     = 1500;
            bs     = 50;
            bw     = 100;
            bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
            mbins  = mean(bins,2)/1000;


            % get ramp on, ramp off time, estimated by the piecewise linear
            % function
            if strcmp(Monk, 'Cy')  % cyrus's time
                if      i==1    L = mbins>=0.1 & mbins<=0.6;
                elseif  i==2    L = mbins>=0   & mbins<=1.4;
                elseif  i==3    L = mbins>=0.5 & mbins<=1.5;
                elseif  i==4    L = mbins>=0.5 & mbins<=1.5;
                elseif  i==5    L = mbins>=0   & mbins<=1;
                elseif  i==6    L = mbins>=0   & mbins<=0.6;
                elseif  i==7    L = mbins>=0   & mbins<=0.8;
                elseif  i==8    L = mbins>=0.1 & mbins<=0.55;
                elseif  i==12   L = mbins>=0.1 & mbins<=0.55;
                elseif  i==13   L = mbins>=0   & mbins<=0.5;
                elseif  i==14   L = mbins>=0   & mbins<=1.4;
                elseif  i==16   L = mbins>=0   & mbins<=0.6;
                elseif  i==17   L = mbins>=0.2 & mbins<=1.4;
                elseif  i==18   L = mbins>=0.2 & mbins<=0.9;
                elseif  i==19   L = mbins>=0.1 & mbins<=0.75;
                elseif  i==21   L = mbins>=0.1 & mbins<=1.5;
                elseif  i==22   L = mbins>=0.1 & mbins<=0.5;
                elseif  i==23   L = mbins>=0.1 & mbins<=0.5;
                elseif  i==24   L = mbins>=0   & mbins<=0.4;
                elseif  i==25   L = mbins>=0   & mbins<=0.5;
                elseif  i==28   L = mbins>=0   & mbins<=0.5;
                elseif  i==30   L = mbins>=0.1 & mbins<=0.8;
                elseif  i==31   L = mbins>=0   & mbins<=0.5;
                elseif  i==32   L = mbins>=0.3 & mbins<=0.7;
                elseif  i==33   L = mbins>=0   & mbins<=0.4;
                elseif  i==34   L = mbins>=0   & mbins<=0.6;
                elseif  i==35   L = mbins>=0.1 & mbins<=0.4;
                elseif  i==37   L = mbins>=0.2 & mbins<=1.2;
                elseif  i==38   L = mbins>=0.15& mbins<=0.75;
                elseif  i==40   L = mbins>=0.1 & mbins<=0.5;
                elseif  i==41   L = mbins>=0.1 & mbins<=1.1;
                elseif  i==42   L = mbins>=0   & mbins<=1;
                elseif  i==44   L = mbins>=0.1 & mbins<=1;
                elseif  i==45   L = mbins>=0.1 & mbins<=0.7;
                elseif  i==46   L = mbins>=0   & mbins<=0.3;
                elseif  i==47   L = mbins>=0   & mbins<=0.9;
                elseif  i==48   L = mbins>=0.1 & mbins<=0.4;
                elseif  i==49   L = mbins>=0   & mbins<=0.4;
                elseif  i==51   L = mbins>=0.5 & mbins<=1;
                elseif  i==52   L = mbins>=0   & mbins<=0.4;
                elseif  i==53   L = mbins>=0.4 & mbins<=0.6;
                elseif  i==55   L = mbins>=0.6 & mbins<=0.8;
                elseif  i==56   L = mbins>=0   & mbins<=0.7;
                elseif  i==57   L = mbins>=0   & mbins<=0.7;
                elseif  i==58   L = mbins>=0.1 & mbins<=0.3;
                elseif  i==59   L = mbins>=0.1 & mbins<=0.3;
                elseif  i==60   L = mbins>=0.5 & mbins<=0.9;
                elseif  i==61   L = mbins>=0   & mbins<=0.3;
                elseif  i==64   L = mbins>=0   & mbins<=0.3;
                elseif  i==65   L = mbins>=0.55 & mbins<=0.75;
                elseif  i==66   L = mbins>=0   & mbins<=0.3;
                elseif  i==69   L = mbins>=0.05& mbins<=0.3;
                elseif  i==70   L = mbins>=0.05& mbins<=0.3;
                elseif  i==74   L = mbins>=0.05& mbins<=0.4;
                elseif  i==75   L = mbins>=0.05& mbins<=0.25;
                elseif  i==76   L = mbins>=0.1 & mbins<=0.3;
                elseif  i==77   L = mbins>=0   & mbins<=1;
                elseif  i==78   L = mbins>=0.1 & mbins<=0.8;
                elseif  i==79   L = mbins>=0.1 & mbins<=0.35;
                elseif  i==80   L = mbins>=0.05& mbins<=0.3;
                elseif  i==81   L = mbins>=0   & mbins<=0.3;
                elseif  i==82   L = mbins>=0   & mbins<=0.3;
                elseif  i==83   L = mbins>=0   & mbins<=0.3;
                elseif  i==84   L = mbins>=0   & mbins<=0.3;
                elseif  i==85   L = mbins>=0   & mbins<=0.35;
                elseif  i==86   L = mbins>=0   & mbins<=0.35;
                elseif  i==87   L = mbins>=0   & mbins<=0.3;
                elseif  i==88   L = mbins>=0.05& mbins<=0.3;
                elseif  i==89   L = mbins>=0.05& mbins<=0.3;
                elseif  i==90   L = mbins>=0.05& mbins<=0.8;
                elseif  i==91   L = mbins>=0.05& mbins<=0.3;
                elseif  i==92   L = mbins>=0.1 & mbins<=0.35;
                elseif  i==93   L = mbins>=0.1 & mbins<=0.35;
                elseif  i==94   L = mbins>=0.1 & mbins<=0.3;
                elseif  i==95   L = mbins>=0.1 & mbins<=0.3;
                elseif  i==96   L = mbins>=0.1 & mbins<=0.3;
                elseif  i==97   L = mbins>=0.1 & mbins<=0.4;
                elseif  i==98   L = mbins>=0.05& mbins<=0.35;
                elseif  i==99   L = mbins>=0   & mbins<=0.4;
                elseif  i==100  L = mbins>=0.05& mbins<=0.3;
                elseif  i==101  L = mbins>=0.05& mbins<=0.3;
                elseif  i==102  L = mbins>=0   & mbins<=0.3;
                elseif  i==104  L = mbins>=0   & mbins<=0.35;
                elseif  i==106  L = mbins>=0   & mbins<=0.4;
                elseif  i==107  L = mbins>=0.1 & mbins<=0.7;
                elseif  i==108  L = mbins>=0.1 & mbins<=0.4;
                elseif  i==109  L = mbins>=0.1 & mbins<=0.3;
                elseif  i==110  L = mbins>=0   & mbins<=0.3;
                elseif  i==112  L = mbins>=0   & mbins<=0.3;
                elseif  i==113  L = mbins>=0   & mbins<=0.25;
                elseif  i==114  L = mbins>=0   & mbins<=0.3;
                elseif  i==116  L = mbins>=0   & mbins<=0.3;
                elseif  i==119  L = mbins>=0.1 & mbins<=0.55;
                elseif  i==120  L = mbins>=0   & mbins<=0.35;
                elseif  i==122  L = mbins>=0   & mbins<=0.3;
                elseif  i==123  L = mbins>=0.1 & mbins<=0.6;
                elseif  i==124  L = mbins>=0   & mbins<=0.4;
                elseif  i==125  L = mbins>=0   & mbins<=0.4;
                elseif  i==126  L = mbins>=0   & mbins<=0.35;
                elseif  i==127  L = mbins>=0   & mbins<=0.3;
                elseif  i==128  L = mbins>=0   & mbins<=0.3;
                elseif  i==129  L = mbins>=0   & mbins<=0.3;
                elseif  i==130  L = mbins>=0   & mbins<=0.3;
                elseif  i==131  L = mbins>=0.3 & mbins<=0.6;
                elseif  i==135  L = mbins>=0   & mbins<=0.2;
                elseif  i==136  L = mbins>=0   & mbins<=0.3;
                elseif  i==137  L = mbins>=0   & mbins<=0.3;
                elseif  i==138  L = mbins>=0   & mbins<=0.3;
                elseif  i==139  L = mbins>=0   & mbins<=0.3;
                elseif  i==141  L = mbins>=0   & mbins<=0.3;
                elseif  i==142  L = mbins>=0   & mbins<=0.3;
                elseif  i==146  L = mbins>=0.05& mbins<=0.35;
                elseif  i==147  L = mbins>=0   & mbins<=0.3;
                elseif  i==149  L = mbins>=0   & mbins<=0.2;
                else
                    L = mbins>=0 & mbins<=0.6;
                end

            else

                if      i==1   L = mbins>=0    & mbins<=1.5;
                elseif  i==8   L = mbins>=0    & mbins<=0.7;
                elseif  i==9   L = mbins>=0.3  & mbins<=0.8;
                elseif  i==15  L = mbins>=0.2  & mbins<=1.3;
                elseif  i==20  L = mbins>=0.2  & mbins<=1;
                elseif  i==27  L = mbins>=0.2  & mbins<=1.2;
                elseif  i==28  L = mbins>=0.4  & mbins<=1.2;
                elseif  i==29  L = mbins>=0.5  & mbins<=1.2;
                elseif  i==34  L = mbins>=0    & mbins<=0.4;
                elseif  i==36  L = mbins>=0.3  & mbins<=1.2;
                elseif  i==37  L = mbins>=0.3  & mbins<=1.2;
                elseif  i==38  L = mbins>=0.4  & mbins<=1.1;
                elseif  i==40  L = mbins>=0.5  & mbins<=1;
                elseif  i==42  L = mbins>=0.2  & mbins<=0.7;
                elseif  i==43  L = mbins>=0.2  & mbins<=1.5;
                elseif  i==44  L = mbins>=0.1  & mbins<=0.5;
                elseif  i==47  L = mbins>=0.4  & mbins<=0.75;
                elseif  i==49  L = mbins>=0    & mbins<=0.3;
                elseif  i==50  L = mbins>=0.5  & mbins<=1.1;
                elseif  i==51  L = mbins>=0.7  & mbins<=1.1;
                elseif  i==52  L = mbins>=0.4  & mbins<=1.1;
                elseif  i==56  L = mbins>=0.45 & mbins<=0.7;
                elseif  i==57  L = mbins>=0.1  & mbins<=0.65;
                elseif  i==58  L = mbins>=0    & mbins<=0.75;
                elseif  i==59  L = mbins>=0.6  & mbins<=1.1;
                elseif  i==62  L = mbins>=0.5  & mbins<=1;
                elseif  i==63  L = mbins>=0.5  & mbins<=1;
                elseif  i==64  L = mbins>=0    & mbins<=0.8;
                elseif  i==65  L = mbins>=0.4  & mbins<=0.8;
                elseif  i==67  L = mbins>=0.25 & mbins<=0.7;
                elseif  i==68  L = mbins>=0.25 & mbins<=0.7;
                elseif  i==69  L = mbins>=0.5  & mbins<=0.75;
                elseif  i==70  L = mbins>=0.75 & mbins<=1.25;
                elseif  i==71  L = mbins>=0.45 & mbins<=0.9;
                elseif  i==76  L = mbins>=0.55 & mbins<=0.9;
                elseif  i==77  L = mbins>=0.65 & mbins<=1.2;
                elseif  i==78  L = mbins>=0.2  & mbins<=0.6;
                elseif  i==79  L = mbins>=0.2  & mbins<=0.9;
                elseif  i==85  L = mbins>=0.2  & mbins<=0.5;
                elseif  i==86  L = mbins>=0.2  & mbins<=0.7;
                elseif  i==87  L = mbins>=0.65 & mbins<=1;
                elseif  i==89  L = mbins>=0.2  & mbins<=0.8;
                elseif  i==90  L = mbins>=0.4  & mbins<=1;
                elseif  i==91  L = mbins>=0.4  & mbins<=0.6;
                elseif  i==92  L = mbins>=0.4  & mbins<=0.6;
                elseif  i==93  L = mbins>=0.45 & mbins<=0.8;
                elseif  i==95  L = mbins>=0.2  & mbins<=0.8;

                else
                    L = mbins>=0 & mbins<=1.5;       % Zsa Zsa's time
                end
            end

            ramptime = mbins(L);
            ramptime = ramptime([1 end])'*1000;
            %ramptime = [(ramptime(2)+ramptime(1))/2 ramptime(2)+(ramptime(2)+ramptime(1))/2];
            
            % get rates
            for j = 1:14
                rb = getFIRA_ecodeTimesByName('dot_on',0);
                re = getFIRA_ecodeTimesByName('dot_off',0);

                if any(Lr(:,j))
                    r      = getFIRA_binRate(Lr(:,j), getFIRA_spikeByID(uid(i)), rb, re, rb, ramptime);
                    m(j,i) = nanmean(r);
                    n(j,i) = nansum(Lr(:,j));
                    v(j,i) = nanvar(r);
                end
            end
            
            % return spike rate for each trial, dir and coh
            r      = getFIRA_binRate(ones(length(Lcoh),1), getFIRA_spikeByID(uid(i)), rb, re, rb, ramptime);
            R = [R; i*ones(length(Lcoh),1) getFIRA_ecodesByName('dot_coh') ...
                     Ldir(:,1) r];
        end
    end

    save(savepath, 'm', 'n', 'v', 'R')
else
    load(savepath)
end

%
%
% figure
% lc  = [];
% for i = 1:7
%     lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
% end
%
% hold on
% for i=1:7
%     x = m([i i+7],:);
%     y = v([i i+7],:);
%     y(x<0) = [];
%     x(x<0) = [];
%     x = x(:);
%     y = y(:);
%
%     % get hypothetical variance curve
%     xlo = x(x<=myprctile(x,20));
%     xhi = x(x>=myprctile(x,80));
%
%     N    = 10000;
%     xsim = nans(N,1);
%     ysim = nans(N,1);
%     len  = length(xlo)-1;
%     M    = 20;
%     for j = 1:10000
%         p       = round(M*rand);   % portion drawn from LO group
%         iLo     = round(len*rand(1,p))+1;
%         iHi     = round(len*rand(1,M-p))+1;
%         xsim(j) = mean([xlo(iLo);xhi(iHi)]);
%         ysim(j) = var([xlo(iLo);xhi(iHi)]);
%     end
%     [xsim I] = sort(xsim);
%     ysim     = ysim(I);
%
%     if i<8
%         plot(x,y,'o', 'MarkerFaceColor', lc(i,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 0.5)
%         Lgd = ~isnan(x) & ~isnan(y) & y~=0 & x~=0;
%         plot(xsim, ysim, 'Color', lc(i,:), 'LineWidth', 1)
%     else
%         plot(x,y,'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc(i-7,:), 'MarkerSize', 0.5)
%         Lgd = ~isnan(x) & ~isnan(y) & y~=0 & x~=0;
%         plot(xsim, ysim, 'Color', lc(i-7,:), 'LineWidth', 1)
%     end
% end
% hold off
% set(gca, 'xscale', 'log', 'yscale', 'log', 'xlim', [0.1 2000], 'ylim', [0.1 2000])
%


%
%
%
% figure
% lc  = [];
% for i = 1:7
%     lc = [lc; [0.9 0.9 0.9]-(i-1)/7];
% end
%
% hold on
% for i=1:14
%     x = m(i,:);
%     y = v(i,:);
%     y(x<0) = [];
%     x(x<0) = [];
%     x = x(:);
%     y = y(:);
%
%
%     if i<8
%         plot(x,y,'o', 'MarkerFaceColor', lc(i,:), 'MarkerEdgeColor', 'none')
%     else
%         plot(x,y,'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', lc(i-7,:))
%     end
% end
% hold off
% set(gca, 'xscale', 'log', 'yscale', 'log')
%
% x = m(:,:);
% y = v(:,:);
% y(x<0) = [];
% x(x<0) = [];
% x = x(:);
% y = y(:);
% Lgd = ~isnan(x) & ~isnan(y) & y~=0 & x~=0;
% b = regress(y(Lgd),x(Lgd));


