function collectLCP_FIRA_SubcorticalUnits(base_dir, remake_FIRA)
%% collectLCP_FIRA_SubcorticalUnits(base_dir)
%
%   Make FIRA files for Sidd's LC-pupil data set
%   8/14/14

if nargin < 1 || isempty(base_dir)
    base_dir =  fullfile(dirnames('local'), 'Data', 'Projects', '2013_LCPupil', ...
        'Data', 'Recording');
end
this_dir = pwd;
monks    = {'Cicero', 'Oz'};
sites    = {'IC', 'LC', 'SC', 'SubC'};

%% CHANGE THIS TO REBUILD ALL FIRA FILES FROM NEX
if nargin < 2 || isempty(remake_FIRA)
    remake_FIRA = false;
end

%%
global FIRA

%% For each monkey
for mm = 1:length(monks)
    
    % For each site
    for ss = 1:length(sites)
        
        % change to monkey/site directory, get lists of files for
        % single-/multi- unit data, then return here
        ms_dir = fullfile(base_dir, monks{mm}, sites{ss});
        cd(ms_dir);
        [su_list, mu_list] = feval(str2func(['fileList_' monks{mm} '_' sites{ss}]));
        
        % need to rebuild FIRA from the raw data files
        if remake_FIRA
            % Input directory of nex files
            inDir  = fullfile(ms_dir, 'nex');
           
            % Output directory for mat files
            outDir = fullfile(ms_dir, 'mat');
            
            % make FIRA for single-unit data
            for uu = 1:size(su_list,1)
                bNex(fullfile(inDir, su_list{uu,1}), 'spmRKLCPupil', ...
                    fullfile(outDir, [su_list{uu,1} '_OUT']), ...
                    'all', 'all', 0, 1, [], []); % [4,14:16]
            end
            
            % make FIRA for multi-unit data
            for uu = 1:size(mu_list,1)
                bNex(fullfile(inDir, mu_list{uu,1}), 'spmRKLCPupil', ...
                    fullfile(outDir, [mu_list{uu,1} '_OUT']), ...
                    'all', 'all', 0, 1, [], []); % [4,14:16]
            end
        end
        
        %% Combine su/mu data into "clean" files       
        % First check if all sus have mus ... if not, add to mu_list
        for uu = 1:size(su_list,1)
            if ~any(strncmp(su_list{uu,1}, mu_list(:,1), strfind(su_list{uu,1}, '-')))
                mu_list = cat(1, mu_list, {su_list{uu,1} []});
            end
        end

        % now loop through the files
        num_files = size(mu_list, 1);
        for ff = 1:num_files
            
            % load mu data file .. weird stuff to avoid load error
            in_dir = fullfile(pwd, 'mat');
            D      = dir(fullfile(in_dir, [mu_list{ff,1} '*.mat']));
            A      = load(fullfile(in_dir, D.name));
            if isfield(A, 'FIRA')
                FIRA = A.FIRA;
            elseif isfield(A, 'data')
                FIRA = A.data;
            end
            clear A;
            
            % possibly clear spikes if no multi-unit data given
            if isempty(mu_list{ff,2})
                FIRA.spikes.data = cell(size(FIRA.spikes.data));
            else
                FIRA.spikes.data = FIRA.spikes.data(:, mu_list{ff,2});
            end
            
            % load su data file -- just append su spike data
            basename = mu_list{ff,1}(1:strfind(mu_list{ff,1},'-')-1);
            for ii = find(strncmp(basename, su_list(:,1), length(basename)))'
                D = dir(fullfile(in_dir, [su_list{ii,1} '*.mat']));
                A = load(fullfile(in_dir, D.name));
                if isfield(A, 'data')
                    A.FIRA = A.data;
                end
                if size(A.FIRA.ecodes.data, 1) == size(FIRA.ecodes.data, 1)
                    FIRA.spikes.data = cat(2,FIRA.spikes.data(:,mu_list{ff,2}), ...
                        A.FIRA.spikes.data(:,su_list{ii,2}));
                else
                    disp('DANGER! DANGER! DANGER!')
                end
                clear A;
            end
            disp(sprintf('file %d, %d units', ff, size(FIRA.spikes.data,2)))

            % call cleanup to do the heavy lifting on FIRA
            cleanLCP_FIRA_tmp(fullfile(ms_dir, 'clean', basename));
        end
    end
end

cd(this_dir);



%  for uu = 1:size(su_list,1)
%                 if isempty(dir(fullfile(inDir, [su_list{uu,1} '.nex'])))
%                     disp(sprintf('%s, %s, SU NO : %s', monks{mm}, sites{ss}, [su_list{uu,1} '.nex']))
%                 else
%                     disp(sprintf('%s, %s, SU YES: %s', monks{mm}, sites{ss}, [su_list{uu,1} '.nex']))
%                 end
%             end
%             
%             for uu = 1:size(mu_list,1)
%                 if isempty(dir(fullfile(inDir, [mu_list{uu,1} '.nex'])))
%                     disp(sprintf('%s, %s, MU NO : %s', monks{mm}, sites{ss}, [mu_list{uu,1} '.nex']))
%                 else
%                     disp(sprintf('%s, %s, MU YES: %s', monks{mm}, sites{ss}, [mu_list{uu,1} '.nex']))
%                 end
%             end
% 
%         end
%     end
% end
% 
