function [fits, sems] = getML_psyPerformancecg(fn, recompute)
% use  cGauss to get behavioral threshold

if nargin < 1
    return;
elseif nargin < 2
    recompute = 1;
end

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_psyPerformancecg_' fn(1:2) '.mat'];

if recompute
    a      = getML_txt(fn);
    fname  = a.data{strcmp(a.name,'dat_fn')};
    
    spec   = [1 1 0 0]; 
    fits   = nans(sum(spec),length(fname));
    sems   = nans(sum(spec),length(fname));
   
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
        [fits(:,i),sems(:,i)] = ctPsych_fit(@cGauss, d(:,1:3), d(:,4), [], spec);
        
        th(i) = fits(1,i);
    end

    
    % save
    save(savepath, 'fits', 'sems', 'th')
else

    % load
    load(savepath)
end





