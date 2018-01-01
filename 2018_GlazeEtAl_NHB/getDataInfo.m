function [fileList_, dataDir_, rawDataDir_, sessionsPerSubject_] = ...
   getDataInfo(recompute)
% function fileList_ = getDataInfo(recompute)
%
% Finds raw data files with at least two sessions of task with feedback
%
% returns (output args 2+ are optional):
%  fileList_            ... list of data filenames
%  dataDir_             ... directory name of parsed data
%  rawDataDir_          ... directory name of raw data
%  sessionsPerSubject_  ... number of sessions per subject

%% SET DATA DIRECTORIES HERE
dataDir    = fullfile('data');
rawDataDir = fullfile(dataDir, 'raw');
filename   = 'dataInfo';

%% possibly just load pre-computed list
if (nargin < 1 || ~recompute) && ...
      exist(fullfile(dataDir, filename), 'file')
   
   % load from file
   load(fullfile(dataDir, filename), 'fileList_');
   
else
   
   %% OTHERWISE CHECK DATA FILES
   % get list in order they were used to generate model
   %  simulations
   load(fullfile(dataDir, 'subj_info.mat'), 'subjids');
   
   % Set up data arrays
   num_files          = size(subjids,1);
   fileList_          = cell(num_files, 1);
   Lgood              = false(num_files, 1);
   sessionsPerSubject = nans(num_files, 1);
   
   % loop through the files
   for ff = 1:num_files
      
      % concatenate path and filename
      datafile = ['data_' subjids{ff} '.mat'];
      fullname = fullfile(rawDataDir, datafile);
      
      % check for data file
      if exist(fullname, 'file')
         
         % load it
         load(fullname);
         
         % check for good session
         sessionsPerSubject(ff) = session_ind - 1;
         if sessionsPerSubject(ff)>1 && mean(data1.signaled)>0
            fileList_{ff} = datafile;
            Lgood(ff)     = true;
         end
      end
   end
   
   % save the selected list
   fileList_ = fileList_(Lgood);
   save(fullfile(dataDir, filename), 'fileList_', 'Lgood', 'sessionsPerSubject');
end

if nargout > 1
   dataDir_ = dataDir;
   
   if nargout > 2
      rawDataDir_ = rawDataDir;
      
      if nargout > 3
         sessionsPerSubject_ = sessionsPerSubject(Lgood);
      end
   end
end
