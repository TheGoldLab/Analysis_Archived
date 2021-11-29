function [w, r, roc, wi, roci] = getML_tuningProperties(fn,recompute)
% get tuning properties for MT and LIP cells
% returns:
%   w   - tuning width of the cell estimated using the [vector average von mises] method
%   r   - the vector strength from vector average
%   roc - the predictive index on sac dir/dot dir


savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_tuningProperties_' fn(1:end-4) '.mat'];
if recompute
    if ~isempty(findstr('MT', fn))
        tdname = 'ddir';
        fdname = 'dot_dir';
        teb    = 'dot_on';
        tee    = 'dot_off';

    elseif ~isempty(findstr('LIP', fn))
        tdname = 'trg_dir';
        fdname = 'trg1_dir';
        teb    = 'trg_off';
        tee    = 'fp_off';

    else
        return;
    end

    a = getML_txt(fn);
    if findstr(fn, 'TRain')
        fname  = a.data{find(strcmp(a.name, 'tune_fn'))};
        task   = a.data{find(strcmp(a.name, 'task_tune'))}; % passive fixation 0 or 3, memory saccade either 2 or 5
    else
        fname  = a.data{find(strcmp(a.name, 'dat_fn'))};
        task   = [];
    end
    
    usable = a.data{find(strcmp(a.name, 'usable'))};
    d1     = a.data{find(strcmp(a.name, tdname))};
    uid    = a.data{find(strcmp(a.name, 'uid'))};
    
    th   = nans(length(fname),1);
    w    = nans(length(fname),2);
    wi   = nans(length(fname),2);  
    r    = nans(length(fname),1);
    roc  = nans(length(fname),1);
    roci = nans(length(fname),1);
    

    % get data
    global FIRA
    for i = 1:length(fname)
        if usable(i) == 1
            openFIRA(fname{i});
            fprintf('%d: %s\n', i, fname{i})

            if ~isnan(fname{i}) & ~isempty(getFIRA_spikeByID(uid(i))) & ~(strcmp(fn, 'ZZTRain_MT.txt') & i==43) %troublesome session
                % remove trials other than 99% for MT pre-training data
                if findstr(fn, 'PRe')
                    Lcoh = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('dot_coh'))>99;
                    FIRA.ecodes.data = FIRA.ecodes.data(Lcoh,:);
                end
                
                % get rate (select only correct trials and for task)
                trials          = find(~isnan(FIRA.ecodes.data(:,1)));            % remove nans
                [Ld, Ud]        = selectFIRA_trialsByUniqueID(fdname, trials);    % select for motion/target dir
                Ld(:,isnan(Ud)) = [];
                Ud(isnan(Ud))   = [];
                Lc              = getFIRA_ecodesByName('correct', [], trials)==1; % select for correct trials
                
                if ~isempty(task)
                    Ltsk        = getFIRA_ecodesByName('task', [], trials)==task(i);    % select for task (passive fix or mem sac)
                else
                    Ltsk        = logical(ones(size(Lc)));
                end
                
                bt              = getFIRA_ecodeTimesByName(teb, 0);
                be              = getFIRA_ecodeTimesByName(tee, 0);

                rate            = getFIRA_rate(trials, getFIRA_spikeByID(uid(i)), bt, be);

                % get mean rate for each dir
                y  = nans(length(Ud),1);
                ys = nans(length(Ud),1);
                n  = nans(length(Ud),1);

                for j = 1:length(Ud)
                    y(j)  = nanmean(rate(Ld(:,j)&Lc&Ltsk));
                    ys(j) = nanse(rate(Ld(:,j)&Lc&Ltsk));
                    n(j)  = nansum(Ld(:,j)&Lc&Ltsk);
                end

               
                % get tuning
                L = ~isnan(y) & ~isnan(ys) & n>0;
                [th(i), w(i,1), stat]            = getTuning(Ud(L)*pi/180, y(L), ys(L), n(L), 'VectorAvg');
                r(i)                             = stat(1);
                weight                           = 1./(n(L).^0.5); % weight by the number of trials but not by ys
                weight                           = nanmean(y)*weight./sum(weight); % make weight big enough so that it'll work for likelihood fit
                [thvm, w(i,2), ym, k, ki]        = vonMisesFitW(Ud(L)*pi/180, y(L), weight);
                [tmp wi(i,1)] = vonMises(ki(:,2),Ud(L)*pi/180);
                [tmp wi(i,2)] = vonMises(ki(:,1),Ud(L)*pi/180);
                w(i,:)  = w(i,:)*180/pi;
                wi(i,:) = wi(i,:)*180/pi;
                
                
                % plot and check von mises fit
                if 0
                    cla
                    polar(Ud, y, 'ko')
                    hold on
                    polar(Ud(L), ym, 'k')
                    hold off
                    title(sprintf('%.1f, %.1f, %.1f', thvm*180/pi, w(i,2), k(2)))
                    pause
                end
                
                

                % get roc areas
                % pick trials +/- 25 deg th(i) as pref,
                % and +/- 25 deg of 180+th(i) as null
                ecd   = getFIRA_ecodesByName(fdname);
                dd    = 20;

                % pref
                ang = mod(d1(i)*180/pi,360);
                if ang>dd | ang<360-dd
                    LP = ecd>=ang-dd & ecd<=ang+dd;
                elseif ang>=0 & ang<=dd
                    LP = (ecd>=ang-dd & ecd<=ang+dd) | ecd>=360-dd+ang;
                elseif ang<0 & ang>=360-dd
                    LP = (ecd>=ang-dd & ecd<=ang+dd) | ecd<=ang-360+dd;
                end

                ang = mod(th(i)*180/pi+180,360);
                if ang>dd | ang<360-dd
                    LN = ecd>=ang-dd & ecd<=ang+dd;
                elseif ang>=0 & ang<=dd
                    LN = (ecd>=ang-dd & ecd<=ang+dd) | ecd>=360-dd+ang;
                elseif ang<0 & ang>=360-dd
                    LN = (ecd>=ang-dd & ecd<=ang+dd) | ecd<=ang-360+dd;
                end


                roc1 = rocN(rate(LP), rate(LN), 1000);

                % old way:
                % compute standard error using the formular in Hanley
                % McNeil 1982, The meaning and Use of the Area under a
                % Receiver Operating Characteristic (ROC) Curve. Radiology
                % 143: 29-36, April 1982.
                A  = roc1;
                n1 = sum(LP);
                n2 = sum(LN);
                Q1 = A/(2-A);
                Q2 = 2*A^2/(1+A);

                ri1 = ((A*(1-A) + (n1-1)*(Q1-A^2) + (n2-1)*(Q2-A^2))/(n1*n2))^0.5;

                
                % pref
                ang = d1(i);
                LP  = ecd==ang;
                ang = mod(d1(i)+180,360);
                LN  = ecd==ang;
                
                if any(LP) & any(LN)
                    roc2 = rocN(rate(LP), rate(LN), 1000);
                else
                    roc2 = nan;
                end
                
                % old way:
                % compute standard error using the formular in Hanley
                % McNeil 1982, The meaning and Use of the Area under a
                % Receiver Operating Characteristic (ROC) Curve. Radiology
                % 143: 29-36, April 1982.
                A  = roc2;
                n1 = sum(LP);
                n2 = sum(LN);
                Q1 = A/(2-A);
                Q2 = 2*A^2/(1+A);

                ri2 = ((A*(1-A) + (n1-1)*(Q1-A^2) + (n2-1)*(Q2-A^2))/(n1*n2))^0.5;

                
                if roc1<0.5
                    roc1 = 1-roc1;
                end
                
                if roc2<0.5
                    roc2 = 1-roc2;
                end
                
                
                [roc(i) ind] = nanmax([roc1 roc2]);
                if ind==1
                    roci(i) = ri1;
                else
                    roci(i) = ri2;
                end
                
            end
        end
    end

    save(savepath, 'r', 'w', 'roc', 'wi', 'roci')
else
    load(savepath)
    
end



















