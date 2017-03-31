function [dir_, fnames_] = getLCP_cleanDataDir(monkey, site)

dir_ = fullfile(dirnames('local'), 'Data', 'Projects', '2013_LCPupil', ...
    'Data', 'Recording', monkey, site, 'clean');

if nargout > 1
    D = dir(fullfile(dir_, '*.mat'));
    fnames_ = {D.name};
end