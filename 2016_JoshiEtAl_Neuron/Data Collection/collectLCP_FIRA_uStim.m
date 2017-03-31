function collectLCP_FIRA_uStim(base_dir, remake_FIRA)
%% function collectLCP_FIRA_uStim(base_dir, remake_FIRA)
%
%   Make FIRA files for Sidd's LC-pupil data set
%   8/14/14

% get the appropriate data directory
if nargin < 1 || isempty(base_dir)
    base_dir =  fullfile(dirnames('local'), 'Data', 'Projects', '2013_LCPupil', ...
        'Data', 'uStim');
end
this_dir = pwd;
sites = { ...
    'Cicero' 'LC' {'short'}; ...
    'Cicero' 'IC' {'short'}; ...
    'Cicero' 'SC' {'short'}; ...
    'Oz'     'LC' {'long' 'short'}; ...
    'Oz'     'IC' {'long' 'short'}; ...
    'Oz'     'SC' {'short'}};
num_sites = size(sites,1);

%% CHANGE THIS TO REBUILD ALL FIRA FILES FROM NEX
if nargin < 2 || isempty(remake_FIRA)
    remake_FIRA = false;
end

%%
global FIRA

%% For each site
for ss = 1:num_sites
    
    % for each stim type
    for tt = 1:size(sites{ss,3},2)
        
        % change to monkey/site directory, get lists of files for
        % single-/multi- unit data, then return here
        ms_dir = fullfile(base_dir, sites{ss,1}, sites{ss,2}, sites{ss,3}{tt});
        files = dir(fullfile(ms_dir, 'nex', '*.nex'));
        
        % loop through file list
        for ff = 1:size(files,1)
            disp([ss tt ff])
            
            % need to rebuild FIRA from the raw data files
            if remake_FIRA
                
                bNex(fullfile(ms_dir, 'nex', files(ff).name), 'spmRKLCPupil', ...
                    fullfile(ms_dir, 'mat', files(ff).name(1:end-4)), [], ...
                    'all', 0, 1, [], []); % [4,14:16]
                
            else
                
                % just get the mat file
                load(fullfile(ms_dir, 'mat', files(ff).name(1:end-4)));
                FIRA = data;
                clear data;
            end
            
            % Make the clean file
            cleanLCP_FIRA(fullfile(ms_dir, 'clean', files(ff).name(1:end-4)));
        end
    end
end

% be nice
cd(this_dir)