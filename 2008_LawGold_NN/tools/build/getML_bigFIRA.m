function data = getML_bigFIRA(utxt_fn, type)
% make big FIRA, binned by session
%  utxt_fn : filename of utxt file
%  type    : 'MT', 'LIP', or 'psy', determine whether 'dot_dir' or
%            'trg1_dir' is changed
%            If type is 'ML', MT and LIP data is put together in one signal
%            file.

if ~strcmp(type, 'ML')
    global FIRA

    % get unit info
    utxt   = getML_txt(utxt_fn);
    fname  = utxt.data{:, strcmp(utxt.name,'dat_fn')};
    ses    = utxt.data{:, strcmp(utxt.name,'session')};

    date   = utxt.data{strcmp(utxt.name,'date')};
    date   = datenum(date);
    date   = date - date(1);

    if strcmp(type, 'MT')
        d1       = utxt.data{:, strcmp(utxt.name,'ddir')};
        ec_dname = 'dot_dir';
        usable = utxt.data{:, strcmp(utxt.name,'usable')};
        uid    = utxt.data{:, strcmp(utxt.name,'uid')};
    elseif strcmp(type, 'LIP')
        d1       = utxt.data{:, strcmp(utxt.name,'trg_dir')};
        ec_dname = 'trg1_dir';
        usable = utxt.data{:, strcmp(utxt.name,'usable')};
        uid    = utxt.data{:, strcmp(utxt.name,'uid')};
    elseif strcmp(type,'psy')
        d1       = [];
        ec_dname = 'dot_dir';
    else
        error(sprintf('%s is not an appropriate data type.\n', type))
        return;
    end

    
    % 1   ... session
    % 2   ... session (in "days since first session")
    % 3   ... correct (-2=broken fix, -1=no choice, 0 = error, 1 = correct)
    % 4   ... dot dir (angle, in degrees)
    % 5   ... coherence (pct)
    % 6   ... viewing time (ms)
    % 7   ... saccade latency
    % 8   ... saccade velocity (avg)
    % 9   ... saccade velocity (max)
    % 10  ... saccade accuracy (dist. from target)
    % 11  ... trial begin time
    % 12  ... time to attain fixation
    % 13  ... # rewards (beeps)



    data = [];
    % get big data matrix
    for i = 1:length(fname)
        openFIRA(fname{i});
        fprintf([int2str(i) ': ' fname{i} '\n'])

        
        % get column indexes
        icoh    = getFIRA_ecodeColumnByName('dot_coh');
        idon    = getFIRA_ecodeColumnByName('dot_on');
        idoff   = getFIRA_ecodeColumnByName('dot_off');
        icrt    = getFIRA_ecodeColumnByName('correct');
        idir    = getFIRA_ecodeColumnByName(ec_dname);
        ilat    = getFIRA_ecodeColumnByName('sac_lat');
        iva     = getFIRA_ecodeColumnByName('sac_vavg');
        ivm     = getFIRA_ecodeColumnByName('sac_vmax');
        irew    = getFIRA_ecodeColumnByName('num_rewards');
        itbe    = getFIRA_ecodeColumnByName('trial_begin');
        ifix    = getFIRA_ecodeColumnByName('fix_time');
        itx     = getFIRA_ecodeColumnByName('t1_x');
        ity     = getFIRA_ecodeColumnByName('t1_y');
        it2x    = getFIRA_ecodeColumnByName('t2_x');
        it2y    = getFIRA_ecodeColumnByName('t2_y');
        irx     = getFIRA_ecodeColumnByName('sac_endrx');
        iry     = getFIRA_ecodeColumnByName('sac_endry');
        itd     = getFIRA_ecodeColumnByName('trg1_dir');
        
        
        % get saccade accuracy along the direction of motion
        dev       = nans(length(FIRA.ecodes.data), 1);
        Lcrt      = FIRA.ecodes.data(:,icrt)==1;
        dx        = FIRA.ecodes.data(:,itx)-FIRA.ecodes.data(:,irx);
        dy        = FIRA.ecodes.data(:,ity)-FIRA.ecodes.data(:,iry);
        dev(Lcrt) = sqrt(dx(Lcrt).^2+dy(Lcrt).^2);
        Lcrt      = FIRA.ecodes.data(:,icrt)==0;
        dx        = FIRA.ecodes.data(:,it2x)-FIRA.ecodes.data(:,irx);
        dy        = FIRA.ecodes.data(:,it2y)-FIRA.ecodes.data(:,iry);
        dev(Lcrt) = sqrt(dx(Lcrt).^2+dy(Lcrt).^2);
        
        
        data = [data;...
            repmat(i,size(FIRA.ecodes.data(:,1))) ...
            repmat(date(i),size(FIRA.ecodes.data(:,1))) ...
            FIRA.ecodes.data(:,icrt) ...
            FIRA.ecodes.data(:,idir) ...
            FIRA.ecodes.data(:,icoh) ...
            (FIRA.ecodes.data(:,idoff)-FIRA.ecodes.data(:,idon)) ...
            FIRA.ecodes.data(:,ilat) ...
            FIRA.ecodes.data(:,iva)  ...
            FIRA.ecodes.data(:,ivm)  ...
            dev                      ...
            FIRA.ecodes.data(:,itbe) ...
            FIRA.ecodes.data(:,ifix) ...
            FIRA.ecodes.data(:,irew)];
    end

    
% if type is 'ML', put MT and LIP data together (added 1/11/07)
elseif strcmp(type,'ML')     
    % get behavioral sessions info
    monk = utxt_fn(1:2);
    
    a   = getML_txt([monk 'TRain_psy.txt']);
    fn  = a.data{:, strcmp(a.name,'dat_fn')};
    ses = a.data{:, strcmp(a.name,'session')};

    date = a.data{strcmp(a.name,'date')};
    date = datenum(date);
    date = date - date(1);

    % get MT units info
    mt    = getML_txt([monk 'TRain_MT.txt']);
    fn_m  = mt.data{:, strcmp(mt.name,'dat_fn')};
    ses_m = mt.data{:, strcmp(mt.name,'session')};
    d1_m  = mt.data{:, strcmp(mt.name,'ddir')};
    gd_m  = mt.data{:, strcmp(mt.name,'usable')};
    uid_m = mt.data{:, strcmp(mt.name,'uid')};
    tr_m  = mt.data{:, strcmp(mt.name,'train')};
    
    fn_m(gd_m~=1|tr_m~=1)  = [];   % remove bad cells
    ses_m(gd_m~=1|tr_m~=1) = [];
    d1_m(gd_m~=1|tr_m~=1)  = [];
    uid_m(gd_m~=1|tr_m~=1) = [];
   
    % get LIP units info
    lip    = getML_txt([monk 'TRain_LIP.txt']);
    fn_l  = lip.data{:, strcmp(lip.name,'dat_fn')};
    ses_l = lip.data{:, strcmp(lip.name,'session')};
    d1_l  = lip.data{:, strcmp(lip.name,'trg_dir')};
    gd_l  = lip.data{:, strcmp(lip.name,'usable')};
    uid_l = lip.data{:, strcmp(lip.name,'uid')};
    
    fn_l(gd_l~=1)  = [];   % remove bad cells
    ses_l(gd_l~=1) = [];
    d1_l(gd_l~=1)  = [];
    uid_l(gd_l~=1) = [];
   
    
    global FIRA
    data = [];
    for i = 1:length(fn)
        openFIRA(fn{i});
        fprintf([int2str(i) ': ' fn{i} '\n'])

        % make big matrix
        % 1) for each session create a big matrix with all the behavioral data
        % 2) then find out how many MT/LIP cells in that session
        % 3) fill in data from column 10-27 accordingly

        %% The columns are as follow:
        % 1   ... session
        % 2   ... session (in "days since first session")
        % 3   ... correct (-2=broken fix, -1=no choice, 0 = error, 1 = correct)
        % 4   ... dot dir (angle, in degrees)
        % 5   ... coherence (pct)
        % 6   ... viewing time (ms)
        % 7   ... saccade latency
        % 8   ... saccade velocity (avg)
        % 9   ... saccade velocity (max)
        % 10  ... saccade accuracy (dist. from target)
        % 11  ... trial begin time
        % 12  ... time to attain fixation
        % 13  ... # rewards (beeps)
        % 14  ... MT1 spike rate (sp/s) of 200 ms preceding dots on
        % 15  ... MT1 spike rate during dots
        % 16  ... MT1 spike rate 200 ms preceding fp off
        % 17  ... MT2 spike rate (sp/s) of 200 ms preceding dots on
        % 18  ... MT2 spike rate during dots
        % 19  ... MT2 spike rate 200 ms preceding fp off
        % 20  ... LIP1 spike rate (sp/s) of 200 ms preceding dots on
        % 21  ... LIP1 spike rate during dots
        % 22  ... LIP1 spike rate 200 ms preceding fp off
        % 23  ... LIP2 spike rate (sp/s) of 200 ms preceding dots on
        % 24  ... LIP2 spike rate during dots
        % 25  ... LIP2 spike rate 200 ms preceding fp off
        % 26  ... LIP3 spike rate (sp/s) of 200 ms preceding dots on
        % 27  ... LIP3 spike rate during dots
        % 28  ... LIP3 spike rate 200 ms preceding fp off
        % 29  ... LIP4 spike rate (sp/s) of 200 ms preceding dots on
        % 30  ... LIP4 spike rate during dots
        % 31  ... LIP4 spike rate 200 ms preceding fp off
        %%

         % get column indexes
        icoh    = getFIRA_ecodeColumnByName('dot_coh');
        idon    = getFIRA_ecodeColumnByName('dot_on');
        idoff   = getFIRA_ecodeColumnByName('dot_off');
        icrt    = getFIRA_ecodeColumnByName('correct');
        idir    = getFIRA_ecodeColumnByName(ec_dname);
        ilat    = getFIRA_ecodeColumnByName('sac_lat');
        iva     = getFIRA_ecodeColumnByName('sac_vavg');
        ivm     = getFIRA_ecodeColumnByName('sac_vmax');
        irew    = getFIRA_ecodeColumnByName('num_rewards');
        itbe    = getFIRA_ecodeColumnByName('trial_begin');
        ifix    = getFIRA_ecodeColumnByName('fix_time');
        itx     = getFIRA_ecodeColumnByName('t1_x');
        ity     = getFIRA_ecodeColumnByName('t1_y');
        irx     = getFIRA_ecodeColumnByName('sac_endrx');
        iry     = getFIRA_ecodeColumnByName('sac_endry');
        itdir   = getFIRA_ecodeColumnByName('trg1_dir');
        itb     = getFIRA_ecodeColumnByName('trial_begin');
        ite     = getFIRA_ecodeColumnByName('trial_end');
        itw     = getFIRA_ecodeColumnByName('trial_wrt');
        ifix    = getFIRA_ecodeColumnByName('fix_time');

        
        % get saccade accuracy along the direction of motion
        dx      = FIRA.ecodes.data(:,itx)-FIRA.ecodes.data(:,irx);
        dy      = FIRA.ecodes.data(:,ity)-FIRA.ecodes.data(:,iry);
        Lx      = FIRA.ecodes.data(:,itx)<0;
        td      = FIRA.ecodes.data(:,itd);
        dev     = cos(pi/180*td).*dx-sin(pi/180*td).*dy;
        dev(Lx) = -1*dev(Lx);
        
        PSY = [repmat(i,size(FIRA.ecodes.data(:,1))) ...
                repmat(date(i),size(FIRA.ecodes.data(:,1))) ...
                FIRA.ecodes.data(:,icrt) ...
                FIRA.ecodes.data(:,idir) ...
                FIRA.ecodes.data(:,icoh) ...
                (FIRA.ecodes.data(:,idoff)-FIRA.ecodes.data(:,idon)) ...
                FIRA.ecodes.data(:,ilat) ...
                FIRA.ecodes.data(:,iva)  ...
                FIRA.ecodes.data(:,ivm)  ...
                dev                      ...
                FIRA.ecodes.data(:,itbe) ...
                FIRA.ecodes.data(:,ifix) ...
                FIRA.ecodes.data(:,irew)];
        
        % get columns to find begin trial (for MT and LIP cells)
        var  = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('var'));
        base = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('base'));
        crt  = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'));
        
        
        
        % get MT neural data
        MTs = nans(length(tw),6);
        I   = find(ses_m==ses(i));
        if ~isempty(I)
            for j = 1:length(I)
                % open FIRA file
                openFIRA(fn_m{I(j)})
                %fprintf([int2str(i) ': ' fn_m{I(j)} '\n'])

                
                trials = [1:size(FIRA.ecodes.data,1)]';
                % spike rate (sp/s) of 200 ms preceding dots on
                r1 = getFIRA_rate(trials,...
                    getFIRA_spikeByID(uid_m(I(j))),...
                    getFIRA_ecodesByName('dot_on')-200,...
                    getFIRA_ecodesByName('dot_on'));
                % spike rate during dots
                r2 = getFIRA_rate(trials,...
                    getFIRA_spikeByID(uid_m(I(j))),...
                    getFIRA_ecodeTimesByName('dot_on',  0),...
                    getFIRA_ecodeTimesByName('dot_off', 0));
                % spike rate 200 ms preceding fp off
                r3 = getFIRA_rate(trials,...
                    getFIRA_spikeByID(uid_m(I(j))),...
                    getFIRA_ecodeTimesByName('fp_off', -200),...
                    getFIRA_ecodeTimesByName('fp_off', 0));
                
                
                % figure out where to put data (note that some cells emerge
                % from the background in the middle of the session, so need
                % to figure out when it appeared from the absolute time
                % (i.e. trial_wrt)
                var_m = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('var'));
                crt_m = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'));
                base_m = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('base'));
                      
                Ibe  = find(var==var_m(1) & crt==crt_m(1) & base==base_m(1));                
                count = 0;
                while length(Ibe)>1
                    count = count+1;
                    Ibe  = find(var==var_m(count+1) & crt==crt_m(count+1) & base==base_m(count+1));
                  end
                Ibe = Ibe-count;
               
                Iend = Ibe+length(trials)-1;
               
                MTs(Ibe:Iend,(j-1)*3+1:(j-1)*3+3) = [r1 r2 r3];
            end
        end

        
        % get LIP neural data
        LIPs = nans(length(tw),12);
        I    = find(ses_l==ses(i));
        if ~isempty(I)
            for j = 1:length(I)
                % open FIRA file
                openFIRA(fn_l{I(j)})
                %fprintf([int2str(i) ': ' fn_l{I(j)} '\n'])

                trials = [1:size(FIRA.ecodes.data,1)]';
                % spike rate (sp/s) of 200 ms preceding dots on
                r1 = getFIRA_rate(trials,...
                    getFIRA_spikeByID(uid_l(I(j))),...
                    getFIRA_ecodesByName('dot_on')-200,...
                    getFIRA_ecodesByName('dot_on'));
                % spike rate during dots
                r2 = getFIRA_rate(trials,...
                    getFIRA_spikeByID(uid_l(I(j))),...
                    getFIRA_ecodeTimesByName('dot_on',  0),...
                    getFIRA_ecodeTimesByName('dot_off', 0));
                % spike rate 200 ms preceding fp off
                r3 = getFIRA_rate(trials,...
                    getFIRA_spikeByID(uid_l(I(j))),...
                    getFIRA_ecodeTimesByName('fp_off', -200),...
                    getFIRA_ecodeTimesByName('fp_off', 0));
                
                
                % figure out where to put data (note that some cells emerge
                % from the background in the middle of the session, so need
                % to figure out when it appeared from the absolute time
                % (i.e. trial_wrt)
                var_l  = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('var'));
                crt_l  = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('correct'));
                base_l = FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('base'));
                
                Ibe  = find(var==var_l(1) & crt==crt_l(1) & base==base_l(1));
                count = 0;
                while length(Ibe)>1
                    count = count+1;
                    Ibe   = find(var==var_l(count+1) & crt==crt_l(count+1) & base==base_l(count+1));
                end
                Ibe = Ibe-count;
                
                Iend = Ibe+length(trials)-1;
               
                LIPs(Ibe:Iend,(j-1)*3+1:(j-1)*3+3) = [r1 r2 r3];               
            end
        end
        data = [data; PSY MTs LIPs];
    end
    
else
    data = [];
end
