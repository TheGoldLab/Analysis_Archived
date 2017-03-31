% checkLCP_files
%
% utility for going through raw data files to check for single & multi-unit
%   spiking data

base_dir =  fullfile(dirnames('local'), 'Data', 'Projects', '2013_LCPupil', ...
    'Data', 'Raw');
monks    = {'Cicero', 'Oz'};
sites    = {'IC', 'LC', 'SC'};

DBdir    = fullfile(base_dir, 'DB');

%% For each monkey
for mm = 1:length(monks)
    
    % For each site
    for ss = 1:length(sites)
        
        % change to monkey/site directory, get lists of files for
        % single-/multi- unit data, then return here
        ms_dir = fullfile(base_dir, monks{mm}, sites{ss});
        cd(ms_dir);
        [su_list, mu_list] = feval(str2func(['fileList_' monks{mm}(1:2) '_' sites{ss}]));
        
        % Input directory of nex files
        inDir  = fullfile(base_dir, monks{mm}, sites{ss}, 'nex');

        disp(sprintf('%s, %s', monks{mm}, sites{ss}))
        disp('SINGLE')
        for uu = 1:size(su_list,1)           
            if isempty(dir(fullfile(inDir, [su_list{uu,1} '.nex']))) && ...
                    ~isempty(dir(fullfile(DBdir, [su_list{uu,1} '.nex'])))
                movefile(fullfile(DBdir, [su_list{uu,1} '.nex']), ...
                    fullfile(inDir, [su_list{uu,1} '.nex']));
            end
        end
        disp('MULTI')
        for uu = 1:size(mu_list,1)           
            if isempty(dir(fullfile(inDir, [mu_list{uu,1} '.nex']))) && ...
                    ~isempty(dir(fullfile(DBdir, [mu_list{uu,1} '.nex'])))
                movefile(fullfile(DBdir, [mu_list{uu,1} '.nex']), ...
                    fullfile(inDir, [mu_list{uu,1} '.nex']));
            end
        end
    end
end
