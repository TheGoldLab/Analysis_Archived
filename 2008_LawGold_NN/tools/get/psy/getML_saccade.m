function [sac sacd]= getML_saccade(fn,recompute)
% get saccade metrics
% returns (columns):
%   1 - Drift during fixation along the direction of motion
%   2 - sd of eye positions along motion direction during fixation
%   3 - sd of eye positions in the orthogonal direction to the motion axis during fixation
%   4 - Saccade velocity (avg)
%   5 - Saccade velocity (max)
%   6 - Accuracy of saccade response
%   7 - Saccade latency
%   8 - Saccade duration
%   9 - Time to attain fixation

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_saccade_' fn(1:2) '.mat'];
if recompute
    a     = getML_txt(fn);
    fname = a.data{find(strcmp(a.name, 'dat_fn'))};

    sac   = nans(length(fname),6);
    sacd  = nans(length(fname),6);

    % get data
    global FIRA
    for i = 1:length(fname)
        openFIRA(fname{i});
        fprintf('%d: %s\n', i, fname{i})    
        
        bt  = getFIRA_ecodesByName('dot_on');
        et  = getFIRA_ecodesByName('dot_off');
        dd  = getFIRA_ecodesByName('dot_dir');
        Igd = find(getFIRA_ecodesByName('correct')==1);
              
        dt  = nans(length(bt),1);    % drift
        sdd = nans(length(bt),1);    % sd along motion dir
        sdo = nans(length(bt),1);    % sd perpendicular to motion
        for j = Igd'
            % get indices of horizontal, vertical eye position data
            eyeX = strcmp('horiz_eye', FIRA.analog.name);
            eyeY = strcmp('vert_eye', FIRA.analog.name);

            % check that we have appropriate analog data
            if any(eyeX) && ~isempty(FIRA.analog.data(j, eyeX).values) && ...
                    any(eyeY) && ~isempty(FIRA.analog.data(j, eyeY).values) && ...
                        FIRA.ecodes.data(j, strcmp('fp_off', FIRA.ecodes.name))>0
            
                % get eye positions between dots on and dots off
                bt_samp = ceil((bt(j) - FIRA.analog.data(j, 1).start_time) * FIRA.analog.store_rate(1)/1000);
                et_samp = ceil((et(j) - FIRA.analog.data(j, 1).start_time) * FIRA.analog.store_rate(1)/1000);
                
                if bt_samp<et_samp & et_samp<=length(FIRA.analog.data(j,1).values)
                    hori = FIRA.analog.data(j,eyeX).values(bt_samp:et_samp);
                    vert = FIRA.analog.data(j,eyeY).values(bt_samp:et_samp);
                    
                    % when vert eye position moves -0.5 deg more than the trim mean,
                    % consider it as a bad trial. This is an eye-tracker
                    % problem and happened to Cyrus's data set a lot.
                    if any(vert-mytrimmean(vert,90)<-0.5);
                        break;
                    end
                    
                    % compute sd of eye position during motion viewing
                    % along the direction of motion and along the
                    % orthogonal direction of motion
                    sdd(j) = nanstd(cos(pi/180*dd(j))*hori-sin(pi/180*dd(j))*vert);
                    sdo(j) = nanstd(sin(pi/180*dd(j))*hori+cos(pi/180*dd(j))*vert);                  
                    
                    % compute drift during fixation
%                     dt(j) = dot([hori(end)-hori(1) vert(end)-vert(1)], [cos(dd(j)*pi/180) sin(dd(j)*pi/180)]);
%                     dt(j) = dt(j)./((et(j)-bt(j))/1000);
                    dt(j) = sqrt(sdd(j)^2+sdo(j)^2);

                    if 0
                        subplot(1,3,1)
                        cla
                        hold on
                        plot(hori,vert,'k-')
                        plot(mean(hori),mean(vert),'r.','MarkerSize', 10)
                        plot(mytrimmean(hori,90),mytrimmean(vert,90),'b.','MarkerSize', 10)
                        hold off
                        %scatter(hori,vert)
                        xlim([-1 1])
                        ylim([-1 1])
                        subplot(1,3,2)
                        cla
                        plot(hori)
                        ylim([-1 1])
                        subplot(1,3,3)
                        cla
                        plot(vert)
                        ylim([-1 1])
                        pause
                    end
                else
                    dt(j) = nan;
                    sd(j) = nan;                
                end
    
            else
                dt(j) = nan;
                sd(j) = nan;             
            end
        end
        sac(i,1)  = nanmean(dt);
        sac(i,2)  = nanmean(sdd);
        sac(i,3)  = nanmean(sdo);
        sacd(i,1) = nanse(dt);
        sacd(i,2) = nanse(sdd);       
        sacd(i,3) = nanse(sdo);       
        
        % get saccade metrics
        itd = getFIRA_ecodeColumnByName('trg1_dir');
        id  = getFIRA_ecodeColumnByName('sac_dur');
        iva = getFIRA_ecodeColumnByName('sac_vavg');
        ivm = getFIRA_ecodeColumnByName('sac_vmax');
        itx = getFIRA_ecodeColumnByName('t1_x');
        ity = getFIRA_ecodeColumnByName('t1_y');
        irx = getFIRA_ecodeColumnByName('sac_endrx');
        iry = getFIRA_ecodeColumnByName('sac_endry');
        ifx = getFIRA_ecodeColumnByName('fix_time');
        
        Lc  = getFIRA_ecodesByName('correct')==1;
        
        sac(i,4)  = nanmean(FIRA.ecodes.data(Lc,iva));
        sac(i,5)  = nanmean(FIRA.ecodes.data(Lc,ivm));
        sacd(i,4) = nanse(FIRA.ecodes.data(Lc,iva));
        sacd(i,5) = nanse(FIRA.ecodes.data(Lc,ivm));

        
        %%% get saccade accuracy along the direction of motion
        dx  = FIRA.ecodes.data(:,itx)-FIRA.ecodes.data(:,irx);
        dy  = FIRA.ecodes.data(:,ity)-FIRA.ecodes.data(:,iry);
        td  = FIRA.ecodes.data(:,itd);
        Lgd = abs(dy)<5;    % remove trials with >5 deg vertical deviation.
                                                        % Most of these trials are blinking trials
        Lx  = FIRA.ecodes.data(:,itx)<0;
        dev = cos(pi/180*td).*dx-sin(pi/180*td).*dy;
        dev(Lx) = -1*dev(Lx);
        sac(i,6)  = nanmean(dev(Lc&Lgd));
        sacd(i,6) = nanse(dev(Lc&Lgd));
        
        % get other saccade parameters
        Lc   = getFIRA_ecodesByName('correct')==1;
        ivl  = getFIRA_ecodeColumnByName('sac_lat');
        Linf = isinf(FIRA.ecodes.data(:,ivl));
        sac(i,7)  = mynanmedian(FIRA.ecodes.data(Lc&~Linf,ivl));
        sacd(i,7) = nanse(FIRA.ecodes.data(Lc&~Linf,ivl));
        sac(i,8)  = nanmean(FIRA.ecodes.data(Lc,id));
        sacd(i,8) = nanse(FIRA.ecodes.data(Lc,id));   
        sac(i,9)  = nanmean(FIRA.ecodes.data(Lc,ifx));
        sacd(i,9) = nanse(FIRA.ecodes.data(Lc,ifx));   
        
    end
    
        
    save(savepath, 'sac', 'sacd')
else
    load(savepath)
    
end
