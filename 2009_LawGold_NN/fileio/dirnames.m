function [home_dir_, lab_dir_, current_dir_, tmat_dir_] = dirnames
% [home_dir_, lab_dir_, current_dir_, tmat_dir_] = dirnames

        % home/lab directories
        home_dir_ = fullfile(filesep, 'Users', 'jigold', 'Desktop');
        lab_dir_  = fullfile(filesep, 'Users', 'jigold', 'Desktop');

        % current project directory
        current_dir_ = fullfile(filesep, 'Users', 'jigold', 'Desktop');
        
        % current temporary data storage directory
        tmat_dir_ = fullfile(filesep, 'Users', 'jigold', 'Desktop');
        
return
