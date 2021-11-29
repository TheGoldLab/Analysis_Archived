function [R, m, v] = getML_LIPSlopeDist(Monk, recompute)
% get distribution of the slope of LIP spike rate during the 'ramp'. The begin and end
% time of the ramp is computed separately for each coherence by fitting a
% piecewise linear function to the average rate data within a session.
% This time is then used for computing rate for each trial.
%
% We use this analysis to check whether there's a bimodal distribution in
% the LIP response
%
%
% OUTPUT:
%   R = [cell #, session #, coh, pref_dir, correct, slope];
%   m = mean for each coh
%   v = variance for each coh

%%
savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_LIPSlopeDist_' Monk '.mat'];

if recompute
    %% get data
    a      = getML_txt([Monk 'TRain_LIP.txt']);
    fn     = a.data{strcmp(a.name,'dat_fn')};
    ses    = a.data{strcmp(a.name,'session')};
    usable = a.data{strcmp(a.name,'usable')};
    uid    = a.data{strcmp(a.name,'uid')};
    d1     = a.data{strcmp(a.name,'trg_dir')};

    crtf   = [1];
    nvf    = [0 1];
    dirf   = 0;

    bb     = -200;
    be     = 1500;
    bs     = 25;
    bw     = 100;
    bins   = [[bb:bs:be-bw]' [bb+bw:bs:be]'];
    mbins  = mean(bins,2)/1000;

    [rM, rSD, rN, rMn] = getML_mRate([Monk 'TRain_LIP.txt'], bins, crtf, nvf, dirf, 0);
    m      = rM-rM([8:14 1:7], :, :);    % select between normalized and unnormalized spikes
    sd     = rSD./(rN.^0.5);
    n      = (rN(1:7,:,:)+rN(8:14,:,:))./2;

    coh    = [0 3.2 6.4 12.8 25.6 51.2 99.9]';

    R      = [];
    m      = nans(length(fn), 2*length(coh));  % mean slope
    v      = nans(length(fn), 2*length(coh));  % variance of slope
    slope  = cell(length(fn), 1);
    
    %% get rate
    global FIRA
    for i = 1:length(fn)
        if usable(i)==1
            openFIRA(fn{i})
            fprintf(sprintf('%d: %s\n', i, fn{i}))

         
            % ramp on ramp off time estimated by the piecewise-linear
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


            t   = mbins(L);
            on  = t(1)*1000;
            off = t(end)*1000; 
            
            % plot
            if 0
                cla
                hold on
                plot(mbins, m(ind,:,i), 'ko', 'MarkerFaceColor', 'k')
                plot(mbins(L),  y, 'k-', 'LineWidth', 2)
                hold off
                uinput = input('','s');
                if strcmp(uinput,'q')
                    break;
                end
            end

            if on>off
                on  = off;
                off = on+200;
            end
            

            % get linear slope for individual trial between ramp on and off
            slp        = nans(length(FIRA.ecodes.data(:,1)),1);
            crt        = getFIRA_ecodesByName('correct');
            spikei     = getFIRA_spikeByID(uid(i));
            dotontime  = getFIRA_ecodesByName('dot_on');
            
            % estimate rate with alpha function
            N    = 500;
            tau  = 25;
            afun = alphafun(1:N+1,'one',1,tau)';
            
            % x axis data (time) for fit
            x = [round(on)/1000:0.001:round(off)/1000]';
            A = [ones(size(x)) x];            
            
            for j = 1:length(FIRA.ecodes.data(:,1))
                if crt(j)>=0
                    sp = FIRA.spikes.data{j,spikei}-dotontime(j);
                    sp = sp(sp>=on & sp<=off);

                    if ~isempty(sp)
                        % prepare data to get slope
                        y            = zeros(round(off)-round(on)+1,1);
                        y(round(sp)-round(on)+1) = 1;
                        y            = conv2(y,[zeros(N,1); afun], 'same')*1000/100;
                        

                        % get slope
                        slp_all = A\y;
                        slp(j)  = slp_all(2);
                        
                        % plot
                        if 0
                            cla
                            hold on
                            plot(x, y, 'k')
                            plot(x, slp_all(1)+x*slp_all(2), 'r')
                            hold off
                            uinput = input('','s');
                            if strcmp(uinput,'q')
                                break;
                            end
                        end
                    else
                        slp(j) = 0;
                    end
                end
            end


            TrialInfo = [repmat(i,length(FIRA.ecodes.data),1) repmat(ses(i),length(FIRA.ecodes.data),1) ...
                getFIRA_ecodesByName('dot_coh') round(getFIRA_ecodesByName('trg1_dir'))==d1(i) ...
                getFIRA_ecodesByName('correct')];
            R  = [R; TrialInfo slp];
            
            
            for j = 1:length(coh)
                Lcoh     = getFIRA_ecodesByName('dot_coh')==coh(j);
                Ldir     = round(getFIRA_ecodesByName('trg1_dir'))==d1(i);
                Lcrt     = getFIRA_ecodesByName('correct')==1;
                m(i,j)   = nanmean(slp(Ldir&Lcrt&Lcoh));
                v(i,j)   = nanvar(slp(Ldir&Lcrt&Lcoh));

                m(i,j+7) = nanmean(slp(Ldir&Lcrt&Lcoh));
                v(i,j+7) = nanvar(slp(Ldir&Lcrt&Lcoh));
            end
         end
    end
   
    save(savepath, 'R', 'm', 'v', 'slope')
else
    load(savepath)
end







