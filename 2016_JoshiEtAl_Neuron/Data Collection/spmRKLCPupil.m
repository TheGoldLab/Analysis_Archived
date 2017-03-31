function spmRKLCPupil(func)

% function spmRKcmd(func):  THIS has TARG acquired ECODE, which
% spmRKcmdnew.m does not
% File describing how to interpret data for FIRA
% Use input argument "func" to switch functions:
%   'init'      ... fills FIRA.spm
%   'trial'     ... parse current trial (in FIRA.trial)
%   'cleanup'   ... whatever you need to do at the end
%
% Returns:
%   nada, but can change FIRA
%

% created by jig 10/21/04
% modified from Sharath's spm724rg.m
% LD 2-21-2007
% modified from Long's spmLDdots.m
% RK 3-19-2007
try
global FIRA EC
    declareEC_ecodes_RKLCPupil;
%     declareEC_ecodes_RKcmd;


    %%if called with empty FIRA.spm, fill it and return
    
    if strcmp(func, 'init')

        % useful ecode markers for tagged ecodes
        cb = 7000;
        cm = 6000;
        cx = 8000;

        % FIRA.spm is a struct with fields corresponding to "dataTypes"
        % that have values that are argument lists to the given dataType
        % constructor method (empty means use defaults).
        % (or, if you prefer, an nx2 cell array, with the first column
        % corresponding to the fieldnames, the second to the arglist.)
        %
        % Thus, to create FIRA with dataType "ecodes" including a list of
        % computed values called "one" and "two", use:
        % spm.ecodes = {'compute', {'one' 'value'; 'two' 'value}}
        %
        % buildFIRA_init then converts these descriptors into actual dataType
        % objects used to parse data. Note that the arglists to buildFIRA_init
        % can override the choices of dataTypes here (so you can use the same
        % base spm file to create different FIRAs)
        %
        % see FIRA/build/dataTypes for possible dataTypes
        %     FIRA/build/buildFIRA_init.m for details of how the spm struct is
        %       used
        FIRA.spm = struct( ...
            'trial',  {{ ...
            'startcd', 1005,                             ...
            'startcd_offset', -1,                    ... % offset to start code
            'endcd_offset', -1,                      ... % offset to end code
            'anycd',   [EC.CORRECTCD EC.WRONGCD EC.NOCHCD EC.BRFIXCD], ...
            'allcd',   [] ...
            'textra',  0, ...
            }}, ...
            'ecodes', {{  ...
            'extract', {   ...  % values to extract & save in FIRA{2}: <name> <type> <method> <args>
            'fp_on'     'time'      @getEC_time         {EC.FPONCD      1}; ...
            'eyefix'    'time'      @getEC_time         {EC.EYINWD      1}; ...
            'tgt_on'    'time'      @getEC_time         {EC.TRGC1CD     1}; ...
            'tgt_off'   'time'      @getEC_time         {EC.TARGOFFCD   1}; ...
            'fp_off'    'time'      @getEC_time         {EC.FPOFFCD     1}; ...
            'fdbkon'    'time'      @getEC_time         {EC.FDBKONCD    1}; ...
            'targ_aq'   'time'      @getEC_time         {EC.TRGACQUIRECD 1}; ...
            'beep_on'   'time'      @getEC_time         {EC.BEEPON       1}; ...
            'stim_time' 'time'      @getEC_time         {EC.ELESTM       1}; ...
            'OLcorrect' 'time'      @getEC_time         {EC.CORRECTCD    1}; ...
            'OLerror'   'time'      @getEC_time         {EC.WRONGCD      1}; ...
            'OLncerr'   'time'      @getEC_time         {EC.NOCHCD       1}; ...
            'OLbrfix'   'time'      @getEC_time         {EC.BRFIXCD      1}; ...         
            'fp_x'      'value'     @getEC_tagged       {EC.I_FIXXCD    cb    cm    cx  0.1};   ...
            'fp_y'      'value'     @getEC_tagged       {EC.I_FIXYCD    cb    cm    cx  0.1};   ...
            't1_x'      'value'     @getEC_tagged       {EC.I_TRG1XCD   cb    cm    cx  0.1};   ...
            't1_y'      'value'     @getEC_tagged       {EC.I_TRG1YCD   cb    cm    cx  0.1};   ...
            't2_x'      'value'     @getEC_tagged       {EC.I_TRG2XCD   cb    cm    cx  0.1};   ...
            't2_y'      'value'     @getEC_tagged       {EC.I_TRG2YCD   cb    cm    cx  0.1};   ...
            'taskid'    'id'        @getEC_tagged       {EC.I_TASKIDCD  cb    cm    cx  1};     ...
            'trialid'   'id'        @getEC_tagged       {EC.I_TRIALIDCD cb    cm    cx  1};     ...
            'seed_base' 'value'     @getEC_tagged       {EC.I_DTVARCD   cb    cm    cx  1};     ...
            }, ...
            'compute', { ...                % names/types of fields to compute & save later in FIRA{2}
            'cor_object',   'id'; ...       % whether target is correct or fixation point
            'choice'        'id'; ...       % which object was chosen:  0=fixation point, 1=target 1, -1=no choice
            'correct'       'id'; ...       % offline score: 1=correct, 0=error, -1=nc, -2=brfix
            'OLscore'       'value'; ...    % online score: 1=correct, 0=error, -1=nc, -2=brfix
            'scorematch'    'value' ; ...
            'fix_time'      'value' ; ...   % time to attain fixation from beginning of trial (first fp_on)
            'sac_endx'      'value' ; ...
            'sac_endy'      'value' ; ...
            'sac_lat'       'value' ; ...
            'sac_lat2',     'value' ; ...
            'sac_lat3',     'value' ; ...
            'sac_dur'       'value' ; ...
            'sac_vmax'      'value' ; ...
            'sac_amp1'       'value' ; ...
            'sac_amp2'       'value' ; ...
            'analog_flag'       'id'; ...
            }, ...
            'tmp', { ...  % values to extract & use temporarily: <name> <type> <method> <args>
            }}}, ...
            'spikes',  [], ...
            'analog', [] );

        % parse the trial ... check whether anything is given in FIRA.trial
    elseif strcmp(func, 'trial')

        % get this trial index
        tti = FIRA.raw.trial.good_count;
        taskid = getFIRA_ec(tti, {'taskid'});

        % compute the saccade-related values
        [sacs_, bf_] = getFIRA_saccades(tti, 3, true , 'horiz_eye', 'vert_eye', 'fp_off', 'fp_x', 'fp_y', 1000);
        
        if ~isempty(sacs_);   %%there are saccade parameters

            % set the latency
            setFIRA_ec(tti, 'sac_lat', sacs_(1,1));
                      
            % set the x-coord of the saccade end point
            setFIRA_ec(tti, 'sac_endx', sacs_(1,5));

            % set the y-coord of the saccade end point
            setFIRA_ec(tti, 'sac_endy', sacs_(1,6));

            % set the saccade maximum velocity
            setFIRA_ec(tti, 'sac_vmax', sacs_(1,3));

            % set the saccade amplitude (raw distance)
            setFIRA_ec(tti, 'sac_amp1', sacs_(1,7));  

            % set the saccade amplitude (linear vector)
            setFIRA_ec(tti, 'sac_amp2', sacs_(1,8));
                        
        end


        %set the correct object
        setFIRA_ec(tti, {'cor_object'}, 0); %%correct object is the fixation point

        %set the choice & offline score
        
        if ~isempty('fp_off')       %%%if it got to the fp_off point in the task then it's correct
            setFIRA_ec(tti, {'choice'}, 0);  %%%stayed with the fixation point
            setFIRA_ec(tti, {'correct'}, 1); %%%it's correct
        elseif ~isempty(bf_) ; %%%there's a broken fixation 
                setFIRA_ec(tti, {'choice'}, -1); 
                setFIRA_ec(tti, {'correct'}, -2); 
        end

        %%%GET extra SACCADE TIMES
        [nSacs,m] = size(sacs_);
        if nSacs>1 && ~isinf( sacs_(2,1) ) && ~isnan( sacs_(2,1) )
            setFIRA_ec(tti, 'sac_lat2', sacs_(2,1));

        end

        if nSacs>2 && ~isinf( sacs_(3,1) ) && ~isnan( sacs_(3,1) )
            setFIRA_ec(tti, 'sac_lat3', sacs_(3,1));
   
        end
        
        currentscore = getFIRA_ec(tti,{'correct'});
        if ~isnan(getFIRA_ec(tti, {'OLcorrect'}))
            OLscore = 1;
        elseif ~isnan(getFIRA_ec(tti, {'OLerror'}))
            OLscore = 0;
        elseif ~isnan(getFIRA_ec(tti, {'OLncerr'}))
            OLscore = -1;
        elseif ~isnan(getFIRA_ec(tti, {'OLbrfix'}))
            OLscore = -2;
        else
            OLscore = -3;   %will this ever happen?
        end
        setFIRA_ec(tti, {'OLscore'}, OLscore);
        if currentscore == OLscore
            setFIRA_ec(tti, {'scorematch'}, 1);
        else
            setFIRA_ec(tti, {'scorematch'}, 0);
        end

  

        %% CHECK to see if the trial is a continuous recording!
        %% Some data seem continuous before toggle was set, these data
        %% are probably correct SSD or VGS trials which went 'til end
        %% anyway
        if FIRA.analog.data(tti,3).length < getFIRA_ec(tti, {'trial_end'})
           setFIRA_ec(tti, {'analog_flag'}, 0)
        elseif FIRA.analog.data(tti,3).length >= getFIRA_ec(tti, {'trial_end'})
            setFIRA_ec(tti, {'analog_flag'}, 1)
        end
       

        % %cleanup
    else

    end
    
catch
    evalin('base', 'e=lasterror')
end