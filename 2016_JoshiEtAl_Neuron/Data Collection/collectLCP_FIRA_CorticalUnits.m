function collectLCP_FIRA_CorticalUnits(base_dir, remake_FIRA)
%% collectLCP_FIRA_CorticalUnits(base_dir, remake_FIRA)
%
%   Make FIRA files for Yin's MT/ACC/PCC pupil data
%   9/4/14

if nargin < 1 || isempty(base_dir)
    base_dir =  fullfile(dirnames('local'), 'Data', 'Projects', '2013_LCPupil', ...
        'Data', 'Recording');
end
this_dir = pwd;
sets = { ...
    'Cheetah'   'MT'; ...
    'Cheetah'   'PCC'; ...
    'Sprout'    'ACC'; ...
    'Sprout'    'PCC'; ...
    'Atticus'   'ACC'};

%% CHANGE THIS TO REBUILD ALL FIRA FILES FROM NEX
if nargin < 2 || isempty(remake_FIRA)
    remake_FIRA = false;
end

%%
global FIRA
FIRA = [];

%% For each data set
for ss = 2:size(sets,1)
    
    % change to monkey/site directory, get lists of files for
    % single-/multi- unit data, then return here
    ms_dir = fullfile(base_dir, sets{ss,1}, sets{ss,2});
    cd(ms_dir);
    file_list = feval(str2func(['fileList_' sets{ss,1} '_' sets{ss,2}]));
    
    % loop through file list
    for uu = 1:size(file_list,1)
        
        % need to rebuild FIRA from the raw data files
        if remake_FIRA
            
            bNex(fullfile(ms_dir, 'nex', file_list{uu,1}), ...
                'spmRKLCPupil', [], 'all', 'all', 0, 1, [], []); % [4,14:16]
            
            % shuffle units so multi is first, single are next
            Funits = cat(2, ...
                find(ismember(FIRA.spikes.id, file_list{uu,3})), ...
                find(ismember(FIRA.spikes.id, file_list{uu,2})));
            FIRA.spikes.channel = FIRA.spikes.channel(Funits);
            FIRA.spikes.unit    = FIRA.spikes.unit(Funits);
            FIRA.spikes.id      = FIRA.spikes.id(Funits);
            FIRA.spikes.data    = FIRA.spikes.data(:, Funits);
            
            % Conditionally add 'blank' multi unit
            if isempty(file_list{uu,3})
                FIRA.spikes.channel = cat(2, -1, FIRA.spikes.channel);
                FIRA.spikes.unit    = cat(2, -1, FIRA.spikes.unit);
                FIRA.spikes.id      = cat(2, -1, FIRA.spikes.id);
                FIRA.spikes.data    = cat(2, cell(size(FIRA.spikes.data,1),1), ...
                    FIRA.spikes.data);
            end
            
            % use only trials within prescribed "good" range
            if ~isempty(file_list{uu,4})
                Lgood = false(size(FIRA.ecodes.data,1), 1);
                wrts  = FIRA.ecodes.data(:,4)./1000; % in sec
                for bb = 1:size(file_list{uu,4},1)
                    Lgood(wrts>=file_list{uu,4}(bb,1) & ...
                        wrts<=file_list{uu,4}(bb,2)) = true;
                end
                FIRA.ecodes.data = FIRA.ecodes.data(Lgood,:);
                FIRA.spikes.data = FIRA.spikes.data(Lgood,:);
                FIRA.analog.data = FIRA.analog.data(Lgood,:);
                FIRA.dio         = FIRA.dio(Lgood,:);
            end
            
            % save it!
            saveFIRA(fullfile(ms_dir, 'mat', [file_list{uu,1} '_OUT']));
            
        else
            
            % just get the mat file
            load(fullfile(ms_dir, 'mat', [file_list{uu,1} '_OUT']));
            FIRA = data;
            clear data;
        end
        
        % Make the clean file
        cleanLCP_FIRA(fullfile(ms_dir, 'clean', file_list{uu,1}));
    end
end

% be nice
cd(this_dir);