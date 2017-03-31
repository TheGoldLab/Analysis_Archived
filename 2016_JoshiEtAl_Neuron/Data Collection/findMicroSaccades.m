function sacs_ = findMicroSaccades(Xin, Yin, offset, show_fig)
% function sacs_ = findMicroSaccades(Xin, Yin, offset, show_fig)
%
% Finds the first "num_sacs" saccades in the given X, Y eye movement
%   traces. "A" for acceleration (the other findSaccades uses only
%   velocity information).
%
% Arguments:
%  Xin         ... raw analog horizontal data
%  Yin         ... raw analog vertical data
%  store_rate  ... samples/sec
%  show_fig     ... optional argument to print debugging figure
%
% Returns an nx8 vectors, rows are saccades found, columns are:
%  1 .. onset time
%  2 .. duration (ms)
%  3 .. maximum speed
%  4 .. distance vect (the length of the saccade calculated as the linear vector)

% SACCADE CRITERIA
MIN_VEL = 0.015;    % deg/ms
MIN_DUR = 6;       % ms
MIN_SEP = 10;      % ms ... minimum separation time between events
maxi    = length(Xin)-50; % max index to identify start of event

if nargin < 3 || isempty(offset)
    offset = 0;
end

% Gaussian kernels
global gk1 gk2
if isempty(gk1)
    gk1 = normpdf(-35:35,0,15);
    gk1 = gk1./sum(gk1);
    gk2 = normpdf(-25:25,0,8);
    gk2 = gk2./sum(gk2);
end

%% get more smoothed version to find saccades
XinS1 = conv(Xin,gk1,'same');
YinS1 = conv(Yin,gk1,'same');              

%% compute velocity, use to identify potential saccades
vel1  = sqrt(diff(XinS1).^2+diff(YinS1).^2);
Fvel  = find(vel1>MIN_VEL);
if isempty(Fvel)
    nsac  = 0;
    sacs_ = [];
else
    % potential saccades are continuous runs of thresholded velocity
    FvelD  = diff(Fvel);
    psacs  = [1;find(FvelD>1)+1];
    npsacs = length(psacs);
    pd     = nans(npsacs, 2);
    for ss = 1:npsacs
        if ss < npsacs
            endi = Fvel(psacs(ss+1)-1)+1;
        else
            endi = Fvel(end)+1;
        end
        if ss > 1 && Fvel(psacs(ss))-pd(ss-1,2)<MIN_SEP
            pd(ss-1,2) = endi;
        else
            pd(ss,:) = [Fvel(psacs(ss)) endi];
        end
    end
    pd(~isfinite(pd(:,1))|pd(:,1)<50|pd(:,1)>maxi|...
        diff(pd,[],2)<MIN_DUR,:) = [];
    
    %% compute stuff to return .. use less smoothed versions
    XinS2 = conv(Xin,gk2,'same');
    YinS2 = conv(Yin,gk2,'same');
    vel2  = sqrt(diff(XinS2).^2+diff(YinS2).^2);
    nsac  = size(pd,1);
    sacs_ = nans(nsac, 4);
    for ss=1:nsac
        
        sacs_(ss,:) = [ ...
            pd(ss,1)+offset, ...                 % onset
            pd(ss,2)-pd(ss,1), ...               % duration
            max(vel2(pd(ss,1):pd(ss,2)-1,1)) ... % max velocity
            sqrt((XinS2(pd(ss,2))-XinS2(pd(ss,1))).^2 + ...   % distance
            (YinS2(pd(ss,2))-YinS2(pd(ss,1))).^2)];
    end
end

%% show stuff
if nargin > 3 && show_fig
    
    xax = offset+(1:length(Xin))';
    
    subplot(2,1,1); cla reset; hold on;
    plot(xax, Xin);
    plot(xax, XinS1,'b--');
    plot(xax, Yin,'r');
    plot(xax, YinS1,'r--');
    
    if nsac>0
        plot(xax, XinS2,'k');
        plot(xax, YinS2,'k');
        
        for ss = 1:nsac
            plot(sacs_(ss,1).*[1 1], [-1 1], 'g-');
            plot((sacs_(ss,1)+sacs_(ss,2)).*[1 1], [-1 1], 'k--');
        end
    end
    ylim([-1.5 1.5]);
    
    subplot(2,1,2); cla reset; hold on;
    plot(xax(1:end-1), vel1,'k--');
    if nsac>0
        plot(xax(1:end-1), vel2,'k:');
    end
    plot(xax([1 end-1]), MIN_VEL.*[1 1], 'r--');
    ylim([0 0.03]);
end