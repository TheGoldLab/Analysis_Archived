
Monk = 'Cyrus';
Mk   = 'Cy';



%% build cell pairs
%
dataPath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/nex'];
savePath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/mat'];
a = getML_txt(['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/' Mk 'TRain_cellpairs.txt']);

fname = a.data{find(strcmp(a.name, 'dat_fn'))};
tflag = a.data{find(strcmp(a.name, 'taskflag'))};
for i =1:length(a.data{1})
    if ~exist(fname{i})   % only do it if file does not exist
        fprintf('Building %sTRain_psy %d:%s\n', Mk, i, fname{i})
        fn = fname{i};
        tf = tflag(i);
        bNex([dataPath '/' fn(1:end-3) 'nex'], '724j', [savePath '/' fn],...
            'all', [15 16], [], true, tf);
    else
        fprintf('%d: %s has already been made.\n', i, fname{i})
    end
end


%% build asso psy
%
dataPath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/nex'];
savePath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/mat'];
a = getML_txt(['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/' Mk 'TRain_asso.txt']);

fname = a.data{find(strcmp(a.name, 'dat_fn'))};
tflag = a.data{find(strcmp(a.name, 'taskflag'))};
for i = 1:length(a.data{1})
    if ~exist(fname{i})   % only do it if file does not exist
        fprintf('Building %sTRain_psy %d:%s\n', Mk, i, fname{i})
        fn = fname{i};
        tf = tflag(i);
        bNex([dataPath '/' fn(1:end-3) 'nex'], '724j', [savePath '/' fn],...
            'all', [15 16], [], true, tf);
    else
        fprintf('%d: %s has already been made.\n', i, fname{i})
    end
end


%% build asso neuro
%
dataPath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/nex'];
savePath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/mat'];
a = getML_txt(['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/' Mk 'TRain_asso.txt']);

fname = a.data{find(strcmp(a.name, 'neuro_fn'))};
tflag = a.data{find(strcmp(a.name, 'taskflag'))};
for i = 1:length(a.data{1})
    if ~exist(fname{i}) & ~isempty(fname{i}) % only do it if file does not exist
        fprintf('Building %sTRain_psy %d:%s\n', Mk, i, fname{i})
        fn = fname{i};
        tf = tflag(i);
        bNex([dataPath '/' fn(1:end-3) 'nex'], '724j', [savePath '/' fn],...
            'all', [15 16], [], true, tf);
    else
        fprintf('%d: %s has already been made.\n', i, fname{i})
    end
end




%% build psy
%
dataPath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/nex'];
savePath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/mat'];
a = getML_txt(['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/' Mk 'TRain_psy.txt']);

fname = a.data{find(strcmp(a.name, 'dat_fn'))};
tflag = a.data{find(strcmp(a.name, 'taskflag'))};
for i = 1:length(a.data{1})
    if ~exist(fname{i})   % only do it if file does not exist
        fprintf('Building %sTRain_psy %d:%s\n', Mk, i, fname{i})
        fn = fname{i};
        tf = tflag(i);
        bNex([dataPath '/' fn(1:end-3) 'nex'], '724j', [savePath '/' fn],...
            'all', [15 16], [], true, tf);
    else
        fprintf('%d: %s has already been made.\n', i, fname{i})
    end
end



%% build LIP training
%
dataPath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/nex'];
savePath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/mat'];
a = getML_txt(['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/' Mk 'TRain_LIP.txt']);

fname = a.data{find(strcmp(a.name, 'dat_fn'))};
usable= a.data{find(strcmp(a.name, 'usable'))};
tflag = a.data{find(strcmp(a.name, 'taskflag'))};
for i = 1:length(a.data{1}) 
    if ~exist(fname{i})   % only do it if file does not exist
        if usable(i) & ~isempty(fname{i})
            fprintf('Building %sTRain_LIP %d:%s\n', Mk, i, fname{i})
            fn = fname{i};
            tf = tflag(i);
            bNex([dataPath '/' fn(1:end-3) 'nex'], '724j', [savePath '/' fn], 'all', [15 16], [], true, tf);
        end
    else
        fprintf('%d: %s has already been made.\n', i, fname{i})
    end
end


%% build MT training
%
dataPath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/nex'];
savePath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/mat'];
a = getML_txt(['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/' Mk 'TRain_MT.txt']);

fname = a.data{find(strcmp(a.name, 'dat_fn'))};
usable= a.data{find(strcmp(a.name, 'usable'))};
tflag = a.data{find(strcmp(a.name, 'taskflag'))};
for i = 1:length(a.data{1})
    if ~exist(fname{i})   % only do it if file does not exist
        if usable(i)
            fprintf('Building %sTRain_MT %d:%s\n', Mk, i, fname{i})
            fn = fname{i};
            tf = tflag(i);
            bNex([dataPath '/' fn(1:end-3) 'nex'], '724j', [savePath '/' fn], 'all', [15 16], [], true, tf);
        end
    else
        fprintf('%d: %s has already been made.\n', i, fname{i})
    end
end


%% build MT pre-training
%
dataPath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/PRe_train/nex'];
savePath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/PRe_train/mat'];
a = getML_txt(['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/PRe_train/' Mk 'PRe_MT.txt']);

fname = a.data{find(strcmp(a.name, 'dat_fn'))};
usable= a.data{find(strcmp(a.name, 'usable'))};
tflag = a.data{find(strcmp(a.name, 'taskflag'))};
for i =1:length(a.data{1})
    if ~exist(fname{i})   % only do it if file does not exist
        if usable(i)==1
            fprintf('Building %sPRe_MT %d:%s\n', Mk, i, fname{i})
            fn = fname{i};
            tf = tflag(i);
            bNex([dataPath '/' fn(1:end-3) 'nex'], '724j', [savePath '/' fn], 'all', [15 16], [], true, tf);
        end
    else
        fprintf('%d: %s has already been made.\n', i, fname{i})
    end
end



%% build LIP tuning
%
dataPath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/nex'];
savePath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/mat'];
a = getML_txt(['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/' Mk 'TRain_LIP.txt']);

fname = a.data{find(strcmp(a.name, 'tune_fn'))};
usable= a.data{find(strcmp(a.name, 'usable'))};
tflag = a.data{find(strcmp(a.name, 'tuneflag'))};
for i = 1:length(a.data{1})
    if ~exist(fname{i}) % only do it if file does not exist
        if usable(i) & ~isempty(fname{i})
            fprintf('Building %sTRain_LIP %d:%s\n', Mk, i, fname{i})
            fn = fname{i};
            tf = tflag(i);
            bNex([dataPath '/' fn(1:end-3) 'nex'], '724j', [savePath '/' fn], 'all', [15 16], [], true, tf);
        end
    else
        fprintf('%d: %s has already been made.\n', i, fname{i})
    end
end





%% build MT tuning
%
dataPath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/nex'];
savePath = ['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/mat'];
a = getML_txt(['/work/mirror_lab/Data/Physiology/MTLIP/' Monk '/TRain/' Mk 'TRain_MT.txt']);

fname = a.data{find(strcmp(a.name, 'tune_fn'))};
usable= a.data{find(strcmp(a.name, 'usable'))};
tflag = a.data{find(strcmp(a.name, 'tuneflag'))};
for i = 1:length(a.data{1})
    if ~exist(fname{i}) & ~isempty(fname{i})% only do it if file does not exist
        if usable(i)
            fprintf('Building %sTRain_MT %d:%s\n', Mk, i, fname{i})
            fn = fname{i};
            tf = tflag(i);
            bNex([dataPath '/' fn(1:end-3) 'nex'], '724j', [savePath '/' fn], 'all', [15 16], [], true, tf);
        end
    else
        fprintf('%d: %s has already been made.\n', i, fname{i})
    end
end