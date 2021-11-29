function [fits, sems, th] = getML_psyPerformanceDD(fn, recompute)
% use     [fits(:,i),sems(:,i)] = ctPsych_fit(@ddExp3, d(:,1:3), d(:,4));
% to get behavioral threshold

if nargin < 1
    return;
elseif nargin < 2
    recompute = 1;
end


if recompute
    a      = getML_txt(fn);
    fname  = a.data{strcmp(a.name,'dat_fn')};
    
    fits   = nans(3,length(fname));
    sems   = nans(3,length(fname));
    th     = nans(length(fname),1);     % threshold

    % Cyrus's first session of training used only 1500ms trial, so it can't
    % be used to fit the time dependent function (b/c viewing time is not a
    % variable).
    if ~isempty(findstr('Cy', fn))  
        ses = 2:length(fname);
    else
        ses = 1:length(fname);
    end

    % get fits
    global FIRA
    for i = ses
            fprintf('%d: %s\n', i, fname{i})
            openFIRA(fname{i})
            
            %% get data for fit
            % data(1) = coh  (0 .. 1)
            % data(2) = time (fractional seconds)
            % data(3) = dot dir: left (-1) / right (1)
            % data(4) = pct or correct (1) / error (0)
            % data(5) = (optinal) n
            d = [getFIRA_ecodesByName('dot_coh')./100 ...
                (getFIRA_ecodesByName('dot_off')-getFIRA_ecodesByName('dot_on'))./1000 ...
                sign(cos(pi/180.*getFIRA_ecodesByName('dot_dir'))) ...
                getFIRA_ecodesByName('correct')];

            % select trials
            Lgood = d(:,4)>=0 & d(:,2)<=1;
            d     = d(Lgood,:);

            % get fits
            [fits(:,i),sems(:,i)] = ctPsych_fit(@ddExp3, d(:,1:3), d(:,4));

            % get threshold at 1 sec
            c   = [0:0.0001:2]';
            d   = [c ones(size(c)) ones(size(c))];
            val = ddExp3([fits(1:end-1,i); 0],d);   % force lapse to be 0 and compute threshold as if there's no lapse
            [xxx, I] = min(abs(val-0.8161)); % find value closest to 81.61%
            th(i) = c(I);
    end

    
    % save
    s = which('getML_psyPerformanceDD.m');
    save([s(1:end-2) fn(1:end-4) '.mat'], 'fits', 'sems', 'th')
else
    % load
    s = which('getML_psyPerformanceDD.m');
    load([s(1:end-2) fn(1:end-4) '.mat'])
end





