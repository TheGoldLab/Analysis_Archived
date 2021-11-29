%% make big FIRA file with 5 columns
%  [days since the first session, coh, dot_dur, dot_dir, correct]
%%


% load session info
global utxt
global FIRA
utxt = getML_txt('ZZTRain_psy.txt');

fname  = utxt.data{strcmp(utxt.name,'dat_fn')};
date   = utxt.data{strcmp(utxt.name,'date')};
date   = datenum(date);
date   = date - date(1);

data = [];


% get performance
for i = 1:length(fname)
    openFIRA(fname{i});
    fprintf([int2str(i) ': ' fname{i} '\n'])

    % get column indexes
    icoh    = getFIRA_ecodeColumnByName('dot_coh');
    idon    = getFIRA_ecodeColumnByName('dot_on');
    idoff   = getFIRA_ecodeColumnByName('dot_off');
    icrt    = getFIRA_ecodeColumnByName('correct');
    iddir   = getFIRA_ecodeColumnByName('dot_dir');

    data = [data;...
        repmat(date(i),size(FIRA.ecodes.data(:,1))) ...
        FIRA.ecodes.data(:,icrt) ...
        FIRA.ecodes.data(:,iddir) ...
        FIRA.ecodes.data(:,icoh) ...
        (FIRA.ecodes.data(:,idoff)-FIRA.ecodes.data(:,idon))];
end

ZZBigEC = data;