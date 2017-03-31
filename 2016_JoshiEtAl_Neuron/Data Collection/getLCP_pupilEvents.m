function pdat_ = getLCP_pupilEvents(ed)
%
% pupil-triggered spikes/LFP
%
% ed is:
%   dim1: trial
%   dim2: sample
%   dim3: 1 = x, 2 = y, 3 = pupil (raw z-score), 4 = standardized z-score,
%           5 = moving slope
%   [remember first sample is considered time=0 (i.e., wrt fix start time)]
%     = ed(Lgood,:,:);
%
% OUTPUT:
% pdat_ is mxn matrix of pupil event data
%   rows are events, columns are:
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. end time of event (wrt fix start time)
%   4. magnitude at start of event (raw z-score)
%   5. magnitude at end of event (raw z-score)
%   6. magnitude at start of event (corrected z-score)
%   7. magnitude at end of event (corrected z-score)
%   8. time of subsequent max slope
%   9. magnitude of subsequent max slope (corrected z/sample)

pdat_ = []; %cell(ntr, 2); % 2 columns are rising, falling
for tt = 1:size(ed,1)
    
    % find PD peaks - zero crossings of first derivative
    lasti = find(isfinite(ed(tt,:,4)),1,'last');
    pd    = ed(tt,1:lasti,3);
    pdc   = ed(tt,1:lasti,4);
    pds   = ed(tt,1:lasti,5);
    d1pds = sign(pds)';
    Fpeak = find(d1pds==0 | (d1pds.*[d1pds(2:end);nan])==-1);

    %plot(pd, 'k-');
    %plot(pds.*100, 'g-');

    if length(Fpeak)>1
        Ldiff = diff(Fpeak)<75;
        if any(Ldiff)
            for ff = find(Ldiff)'
                Fpeak(ff:ff+1) = nan;
            end
            Fpeak = Fpeak(isfinite(Fpeak));
        end
        
        if length(Fpeak)>1
                        
            for ff = 1:length(Fpeak)-1
                
                % identify subsequent slope
                evt   = pds(Fpeak(ff):Fpeak(ff+1));
                sMaxI = Fpeak(ff) + find(abs(evt)==max(abs(evt)),1);
                pdat_ = cat(1, pdat_, ...
                    [tt, Fpeak(ff)-1, Fpeak(ff+1)-1, ...
                    nanmean(pd(max(1,Fpeak(ff)-76):min(lasti,Fpeak(ff)+75))), ...
                    nanmean(pd(max(1,Fpeak(ff+1)-76):min(lasti,Fpeak(ff+1)+75))), ...
                    pdc(Fpeak(ff)), pdc(Fpeak(ff+1)), ...
                    sMaxI, pds(sMaxI)]);
                
                %if pds(sMaxI) < 0
                %  co = 'r';
                %else
                %  co = 'b';
                %end
                %plot(Fpeak(ff), pd(Fpeak(ff)), 'o', 'Color', co);
                %plot(Fpeak(ff+1), pd(Fpeak(ff+1)), '+', 'Color', co);
                %plot(sMaxI, pd(sMaxI), '*', 'Color', co);
            end
        end
        %         r = input('next')
        %         cla reset; hold on;
    end
end
