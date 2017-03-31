%% collectLCP_FIRA
%
% Script for collecting data from raw .nex files for Joshi et al 2016.
% Note that each routine fills a folder (should already be created) called
%   "clean" located at the same level as each "raw" folder containing the 
%   raw data files
%

%% This is where all the files should be located
base_dir =  fullfile(dirnames('local'), 'Data', 'Projects', '2013_LCPupil', 'Data');

%% Flag determining whether to remake "FIRA" from the raw data files
%   see the assocaited spm*.m file for more details
remake_FIRA = false;

%% Collect spiking data from LC, subC, IC, SC
unit_dir = fullfile(base_dir, 'Recording');
collectLCP_FIRA_SubcorticalUnits(unit_dir, remake_FIRA);

%% Collect spiking data from ACC/PCC
unit_dir = fullfile(base_dir, 'Recording');
collectLCP_FIRA_CorticalUnits(unit_dir, remake_FIRA);
collectLCP_FIRA_ACCv2(unit_dir, remake_FIRA)

%% Collect ustim data
ustim_dir =  fullfile(base_dir, 'uStim');
collectLCP_FIRA_uStim(unit_dir, remake_FIRA);
