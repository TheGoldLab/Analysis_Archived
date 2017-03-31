function [stpddat_, sudat_] = getLCP_stpds(siteData, ptype, utype, num_shuffled)
%   get spike-triggered pupil diameter/slope
%
% INPUT:
%  siteData celery:
%  - siteData{1}: trialsxcols matrix, cols are:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correcdt
%   4 ... beep on time, wrt to fix start time
%   5 ... trial begin time, wrt to fix start time
%   6 ... trial end time, wrt to fix start time
%   7 ... trial wrt time (cpu clock)
%   8 ... LFP index corresponding to fix start time (coded above)
%
%  - siteData{2}: Analog:
%   dim1: trial
%   dim2: sample
%   dim3: 1 = x, 2 = y, 3 = pupil, 4 = std pupil, 5 = 
%   [remember first sample is considered time=0 (i.e., wrt fix start time)]
%     = eyedat(Lgood,:,:);
%
%  - siteData{3}: spikes, re-coded wrt fix start time
%
%  - siteData{4}: LFP
%
%  ptype:           'pupil' (default) or 'slope'
%  utype:           'single' (default), 'multi', or 'both'
%  num_shuffled:   shuffle trials, for stats
%
% OUTPUTS:
%  stpddat_ is:pupil measure
%   1: per unit (multi first)
%   2: NUM_BINS
%   3: mean, sem, shuffled mean
%
% OLD:
%  pdat_ is: min/max
%   1: unit
%   2: min, max
%   3: index (time), p-val Mann-Whitney, p-val shuffled min/max

% Window size for spike-triggered PD
WIN_START  = -1000;
WIN_END    =  2000;
wi         = WIN_START:WIN_END;
WIN_SIZE   = length(wi);
%BSI        = -WIN_START-25+(0:51); % indices for baseline subtraction
% BNUM       = 1000; % for bootstrap
% bdat       = nans(BNUM,2); % min/max bootstraps
% sdat       = nans(WIN_SIZE,4);

%% find no beep trials -- spike-triggered average
LnoBeep    = isnan(siteData{1}(:,4)); % only non-beep trials
FnoBeep    = find(LnoBeep);
num_noBeep = length(FnoBeep);

%% get pupil data ...
if nargin < 2 || isempty(ptype) || strcmp(ptype, 'pupil')
    % default -- pupil diameter
    pdat = siteData{2}(FnoBeep,:,4); % trial/sample
else
    % pupil slope
    pdat = siteData{2}(FnoBeep,:,5); % trial/sample
end    
num_samps = size(pdat,2);

%% get unit type
if nargin < 3 || isempty(utype) || strcmp(utype, 'single')    
    num_units = size(siteData{3},2)-1; % first is always mulit-unit (or empty)
    if num_units < 1
        stpddat_ = [];
        sudat_ = [];
        % pdat_    = [];
        return
    end
    ui = 2:num_units+1;
elseif strcmp(utype, 'multi')
    num_units = 1;
    ui = 1;
else
    num_units = size(siteData{3},2); % first is always mulit-unit (or empty)
    ui = 1:num_units;
end

%% possibly get shuffled trials
if nargin < 4 || isempty(num_shuffled) || num_shuffled <= 0
    num_shuffled = 0;
else
    Fsh = nans(num_noBeep, num_shuffled);
    for ss = 1:num_shuffled
        Fsh(:,ss) = randperm(num_noBeep); % shuffled trials
    end
end

%% Loopy units & collect data
stpddat_ = nans(num_units, WIN_SIZE, 4); % 3:mean/sem/shuffled mean/p
sudat_   = (1:num_units)';
%pdat_    = nans(num_units, 3, 2); % 3: time/MW-p/Shuffled-p; 3: min/max, 
Lug      = true(num_units,1);
for uu = 1:num_units
    
    % get/count spikes
    sps = siteData{3}(LnoBeep,ui(uu));
    
    % count spikes
    tsdat = zeros(num_noBeep, 1);
    for tt = 1:num_noBeep % each trial
        if ~isempty(sps{tt})
            Lg = sps{tt}+WIN_END>=0 & sps{tt}+WIN_START<=num_samps;
            sps{tt} = sps{tt}(Lg);
            tsdat(tt) = sum(Lg);
        end
    end

    % if enough total spikes...
    if sum(tsdat) < 200
        Lug(uu) = false;
    else
        
        % make matrix for this session, this unit
        tmpdat = nans(sum(tsdat), WIN_SIZE);
        ind1   = 1;
        if num_shuffled > 0
            nts   = min(floor(200000/sum(tsdat)), num_shuffled);
            shdat = nans(sum(tsdat).*nts, WIN_SIZE);
            ind2  = 1;
        else
            nts = 0;
        end
        disp([uu ui(uu) sum(tsdat) nts])

        for tt = find(tsdat>0)' % loop through all good trials
            for pp = 1:length(sps{tt}) % each spike
                spi = floor(sps{tt}(pp)) + wi;
                Lti = spi>=0 & spi<num_samps;
                if any(Lti)
                    % save real (baseline subtracted)
                    tmpdat(ind1,Lti) = pdat(tt,spi(Lti)+1);
                    
                    if num_shuffled > 0
                        % save shuffled (baseline subtracted)
                        shdat(ind2+(0:nts-1)',Lti) = ...
                            pdat(Fsh(tt,1:nts)',spi(Lti)+1);
                        ind2 = ind2 + nts;
                    end
                end
                ind1 = ind1 + 1;
            end
        end
        
        % save the mean/sem (after peri-spike mean subtraction)
        %tmpdat           = tmpdat - repmat(nanmean(tmpdat(:,BSI),2), 1, WIN_SIZE);
        stpddat_(uu,:,1) = nanmean(tmpdat);
        stpddat_(uu,:,2) = nanse(tmpdat);
        
        % if shuffled use for p values and such
        if num_shuffled > 0
            for ww = 1:WIN_SIZE
                Lg1 = isfinite(tmpdat(:,ww));
                Lg2 = isfinite(shdat(:,ww));
                stpddat_(uu,ww,3) = mean(shdat(Lg2,ww));
                stpddat_(uu,ww,4) = ranksum(tmpdat(Lg1,ww), shdat(Lg2,ww));
            end
        end
    end
end
sudat_   = sudat_(Lug);
stpddat_ = stpddat_(Lug,:,:);

%             % get shuffled min/max
%             nsh = size(shdat,1); % total number of shuffled events
%             npb = min(size(tmpdat,1), 5000); % number to get per bootstrap       
%             for xx = 1:BNUM
%                 td        = shdat(randperm(nsh,npb),:);
%                 Lg        = ~isfinite(td);
%                 td(Lg)    = 0;
%                 sdat(:,3) = sum(td);
%                 sdat(:,4) = sum(Lg);
%                 sdat(:,3) = sdat(:,3)./sdat(:,4) - ...
%                     (sdat(:,1)-sdat(:,3))./(sdat(:,2)-sdat(:,4));
%                 bdat(xx,:) = [min(sdat(:,3)) max(sdat(:,3))];
%             end
% 
%             % find peak, p-value at peak
%             vals = stpddat_(uu,:,1) - stpddat_(uu,:,3);
%             
%             % min
%             Vmn  = nanmin(vals);
%             Fmn  = find(vals==Vmn,1);
%             Lg1  = isfinite(tmpdat(:,Fmn));
%             Lg2  = isfinite(shdat(:,Fmn));
%             pdat_(uu,:,1) = [ ...
%                 wi(Fmn) ...
%                 ranksum(tmpdat(Lg1,Fmn), shdat(Lg2,Fmn)) ...
%                 sum(bdat(:,1)<=Vmn)./BNUM];
% 
%             % max
%             Vmx  = nanmax(vals);
%             Fmx  = find(vals==Vmx,1);
%             Lg1  = isfinite(tmpdat(:,Fmx));
%             Lg2  = isfinite(shdat(:,Fmx));
%             pdat_(uu,:,2) = [ ...
%                 wi(Fmx) ...
%                 ranksum(tmpdat(Lg1,Fmx), shdat(Lg2,Fmx)) ...
%                 sum(bdat(:,2)>=Vmx)./BNUM];
%             
%         else
%             vals = stpddat_(uu,:,1);
%             
%             % min
%             Fmn  = find(vals==nanmin(vals),1);
%             Lg   = isfinite(tmpdat(:,Fmx));
%             pdat_(uu,:,1) = [wi(Fmn) signtest(tmpdat(Lg,Fmn)) nan];
%             
%             % max
%             Fmx  = find(vals==nanmax(vals),1);
%             Lg   = isfinite(tmpdat(:,Fmx));
%             pdat_(uu,:,2) = [wi(Fmx) signtest(tmpdat(Lg,Fmx)) nan];
%         end        
%pdat_    = pdat_(Lug,:,:);

% 
%         subplot(num_units, 1, uu); cla reset; hold on;
%         plot(wi([1 end]), [0 0], 'k:');
%         plot([0 0], [-.4 .4], 'k:');
%         co = {'r' 'g' 'b'};
%         for ii = 1:3
%             stdif = stpddat_(:,uu,1,ii)-stpddat_(:,uu,2,ii);
%             if ii == 2
%                 stdif = stdif.*100;
%             end
%             plot(wi, stdif, '-', 'Color', co{ii});
%             Lp = stpddat_(:,uu,3,ii)<0.05;
%             plot(wi(Lp), stdif(Lp), 'k.');
%         end
%         ylim(max(abs(stdif)).*[-1.1 1.1]);
%     end
% end
% 
% % 
% %             bdat_(bb,uu,3,ii) =
% %         end
% %             stpddat_   = nans(nbins, num_units, 3); % 3:mean/sem; 4:raw/shuffled
% % 
% %             % save the mean/sem/std
% %             stpddat_(:,uu,:,ii) = [ ...
% %                 nanmean(tmpdat(:,:,ii))', ...
%                 nanse(tmpdat(:,:,ii))'];
%         end
% sdat_      = nans(5, num_units);
% bdat_      = nans(nbins, num_units, 3, 2);

%  sdat_ is:
%   1: peak time/peak pval/mean dif/sem dif/std dif
%   2: per unit

            % stats
%             smpd = nanrunmean(stpddat_(:,uu,1,ii), 50); % smoothed mean pd
%             peak_time = find(smpd(-WIN_START+1:end-200)==...
%                 max(smpd(-WIN_START+1:end-200)),1);
%             Lg = isfinite(tmpdat(:,-WIN_START+peak_time));
%             [H,P] = ttest(tmpdat(Lg,-WIN_START+peak_time));
%             sdat_(:,uu) = [peak_time, P, ...
%                 mean(tmpdat(Lg,-WIN_START+peak_time)), ...
%                 nanse(tmpdat(Lg,-WIN_START+peak_time)), ...
%                 std(tmpdat(Lg,-WIN_START+peak_time))];
            
%         subplot(num_units, 1, uu); cla reset; hold on;
%         hold on
%         plot(tax, bs(:,1), 'r-')
%         plot(tax, bs(:,2), 'g-')
%         plot(stpddat_(:,uu,1), 'b-')
%         ylim([-0.02 0.02])
% 
%     end
% end

                %                     if false
                %                         tax = 1:num_samps;
                %                         subplot(2,1,1); cla reset; hold on;
                %                         plot(pdat(FnoBeep(tt),:,1));
                %                         plot(sps{tt}([pp pp]), [-2 2], 'r-')
                %                         plot([sps{tt}(pp)+1000 sps{tt}(pp)+1000], [-2 2], 'r--')
                %                         if any(Lti)
                %                             plot(tax(spi(Lti)+1), pdat(FnoBeep(tt),spi(Lti)+1), 'g-');
                %                             axis([-3000 4500 -3 3]);
                %                             subplot(2,1,2); cla reset; hold on;
                %                             plot(wi(Lti), tmpdat(ind2,Lti,1), 'r-');
                %                             plot([0 0], [-3 3], 'k:');
                %                             plot([1000 1000], [-3 3], 'k:');
                %                             axis([-100 1100 -3 3]);
                %                         else
                %                             axis([-3000 4500 -3 3]);
                %                             subplot(2,1,2); cla reset; hold on;
                %                             axis([-100 1100 -3 3]);
                %                         end
                %                         pause(0.1); %r=input('next')
                %                     end


