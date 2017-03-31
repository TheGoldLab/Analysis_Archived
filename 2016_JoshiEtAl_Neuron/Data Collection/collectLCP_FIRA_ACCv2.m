function collectLCP_FIRA_ACCv2(base_dir, remake_FIRA)
% function collectLCP_FIRA_ACCv2(base_dir, remake_FIRA)
%
%   Make FIRA files for Sidd's Sprout ACC data (Neuron v2)
%
%   10/14/15

%% get current directory and directory containing the raw data
this_dir = pwd;
if nargin < 1 || isempty(base_dir)
    base_dir =  fullfile(dirnames('local'), 'Data', 'Projects', '2013_LCPupil', ...
        'Data', 'Recording');
end
% need just this directory of raw data
base_dir = fullfile(base_dir, 'Sprout', 'ACC');

%% CHANGE THIS TO REBUILD ALL FIRA FILES FROM NEX
if nargin < 2 || isempty(remake_FIRA)
    remake_FIRA = false;
end

%%
global FIRA
FIRA = [];

% change to monkey/site directory, get lists of files for
% single-/multi- unit data, then return here
cd(base_dir);
file_list = feval('Sprout_2015_Sidd_FileNames');

% loop through file list
for uu = 1:size(file_list,1)
    
    % need to rebuild FIRA from the raw data files
    if remake_FIRA
        
        % multi-unit
        if ~isempty(file_list{uu,2})
            bNex(fullfile(base_dir, 'nexv2', ['Sprout_' file_list{uu,1}], file_list{uu,2}), ...
                'spmRKLCPupil', [], 'all', 'all', 0, 1, [], []); % [4,14:16]
            FIRA_multi = FIRA;
            filename   = file_list{uu,2};
        else
            FIRA_multi = [];
        end
        
        % single unit
        if ~isempty(file_list{uu,3})
            bNex(fullfile(base_dir, 'nexv2', ['Sprout_' file_list{uu,1}], file_list{uu,3}), ...
                'spmRKLCPupil', [], 'all', 'all', 0, 1, [], []); % [4,14:16]
            
            Funits = find(ismember(FIRA.spikes.unit, file_list{uu,4}));
            FIRA.spikes.channel = FIRA.spikes.channel(Funits);
            FIRA.spikes.unit    = FIRA.spikes.unit(Funits);
            FIRA.spikes.id      = FIRA.spikes.id(Funits);
            FIRA.spikes.data    = FIRA.spikes.data(:, Funits);
            filename = file_list{uu,3};
        else
            FIRA = [];
        end
        
        % pssibly combine
        if ~isempty(FIRA_multi) && ~isempty(FIRA)
            FIRA.spikes.channel = cat(2, -1, FIRA.spikes.channel);
            FIRA.spikes.unit    = cat(2, -1, FIRA.spikes.unit);
            FIRA.spikes.id      = cat(2, -1, FIRA.spikes.id);
            FIRA.spikes.data    = cat(2, FIRA_multi.spikes.data(:,1), ...
                FIRA.spikes.data);
            
        elseif ~isempty(FIRA_multi)
            FIRA = FIRA_multi;
        end
        
        % save it!
        saveFIRA(fullfile(base_dir, 'mat', [filename '_OUT']));
        
    else
        
        % just get the mat file
        if ~isempty(file_list{uu,3})
            filename = file_list{uu,3};
        else
            filename = file_list{uu,2};
        end
        load(fullfile(base_dir, 'mat', [filename '_OUT']));
        FIRA = data;
        clear data;
    end
    
    % Make the clean file
    cleanLCP_FIRA(fullfile(base_dir, 'clean', filename));
end

cd(this_dir);
