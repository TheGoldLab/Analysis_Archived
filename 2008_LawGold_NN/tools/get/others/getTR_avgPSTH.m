function [ym yse yn y n] = getTR_avgPSTH(fn, fn_type, Bses, binsize, normalized, recompute, savefn)

% check input
if nargin<7    
    % default save path
    savefn = which('getTR_avgPSTH.m');
    savefn = [savefn(1:end-2) '.mat'];
end

if nargin<6
    recompute = 1;
end

if nargin<5
    normalized = 0;
end

if nargin<4
    binsize = 100;
end

if nargin<3
    return;
end



% spreadsheet and FIRA column names for direction
if strcmp(fn_type, 'MT')
    utxtcol = 'ddir';
    firacol = 'dot_dir';
elseif strcmp(fn_type, 'LIP')
    utxtcol = 'trg_dir';
    firacol = 'trg1_dir';
elseif strcmp(fn_type, 'PReMT')
    utxtcol = 'ddir';
    firacol = 'dot_dir';
end



% define global variables
global FIRA
global utxt
utxt = getML_txt(fn);

fname  = utxt.data{strcmp(utxt.name,'dat_fn')};
usable = utxt.data{strcmp(utxt.name,'usable')};
uid    = utxt.data{strcmp(utxt.name,'uid')};
d1  = utxt.data{strcmp(utxt.name, utxtcol)};
if ~strcmp(fn_type, 'PReMT')
    ses    = utxt.data{strcmp(utxt.name,'session')};
else
    ses    = [1:length(fname)]';
    Bses   = [1 length(fname)+1];
end

ec0_     = 'dot_on';
ec1_     = 'dot_off';
os0_     = -200;
os1_     = 0;


maxdd = 1500;
% if recompute is set true or y and n doesn't exist, compute
if recompute | ~exist(savefn)
    % get psth for each cell, assume that max dot_dur is 1500ms
    y = nans(ceil((maxdd-os0_)/binsize),14,length(utxt.data{1}));
    n = nans(ceil((maxdd-os0_)/binsize),14,length(utxt.data{1}));
    for i = 1:length(utxt.data{1})
        if usable(i)==1
            % load data
            openFIRA(fname{i});
            fprintf('%s\n',fname{i})

            % get selection arrays
            % for trials.
            trials = find(~isnan(FIRA.ecodes.data(:,1)));
            % for stimulus conditions (trg1_bdir, dot_coh and task)
            [Lcoh, Ucoh]   = selectFIRA_trialsByUniqueID('dot_coh');
            Ltsk           = ismember(FIRA.ecodes.data(:,getFIRA_ecodeColumnByName('task')), [3,6]);
            [Lcrt, Ucrt]   = selectFIRA_trialsByUniqueID('correct');
            Lcrt           = Lcrt(:,Ucrt==1);   % correct trials only
            [Ltrgd, Utrgd] = selectFIRA_trialsByUniqueID(firacol);

            ind_d = [find(round(Utrgd)==d1(i)) find(round(Utrgd)==mod(d1(i)+180,360))];
            L     = zeros(length(Ltrgd), 2*length(Ucoh));
            for k = 1:2
                for m = 1:length(Ucoh)
                    L(:,(k-1)*length(Ucoh)+m) = Ltrgd(:,ind_d(k)) & Lcoh(:,m) & Lcrt & Ltsk;
                end
            end

            % get psth
            [a b c d] = getFIRA_rasterPSTH( trials,...
                L,...
                getFIRA_spikeByID(uid(i)),...
                getFIRA_ecodeTimesByName(ec0_, os0_, trials),...
                getFIRA_ecodeTimesByName(ec1_, os1_, trials),...
                getFIRA_ecodeTimesByName(ec0_, 0, trials),...
                binsize, 0);
            
            % normalize responses
            if normalized
                % subtract rate with the baseline response 200ms before dots on
                bins = os0_:binsize:maxdd-binsize;
                bl   = mean(mean(c(bins<0,:)));
                c    = c-bl;
                % normalize rate with the max of preferred 99.9% and 51.2% responses
                maxr = max(max(c(bins(:,2)<=1000,6:7)));  % assume that col 6 and 7 are col for 51.2 and 99.9%
                                            % and only use values under
                                            % 1000ms b/c long vt data are noisy
                c    = c./maxr;
            end
            
            if length(Ucoh)==7
                y(:,:,i) = c(1:size(y,1),1:size(y,2));
                n(:,:,i) = d(1:size(n,1),1:size(n,2));
            else    % if any of the coh is missing
                Lcc = ismember([0 3.2 6.4 12.8 25.6 51.2 99.9], Ucoh);
                Lcc = [Lcc Lcc]';   % make it selection array for both preferred and null
                y(:,Lcc,i) = c(1:size(y,1),:);
                n(:,Lcc,i) = d(1:size(n,1),:);
            end 
        end
    end
    save(savefn, 'y', 'n') 
else
    load(savefn)
end

% compute average response
ym  = nans(ceil((maxdd-os0_)/binsize),14,size(Bses,1));
yse = nans(ceil((maxdd-os0_)/binsize),14,size(Bses,1));
yn  = nans(1,size(Bses,1));

for i=1:size(Bses,1)
    Lses       = ses>=Bses(i,1) & ses<Bses(i,2);
    ym(:,:,i)  = nanmean(y(:,:,Lses),3);
    yn(i)      = sum(Lses & ~isnan(shiftdim(y(1,1,:),2)));
    yse(:,:,i) = nanstd(y(:,:,Lses), ones(1,sum(Lses)), 3)./sum(~isnan(y(:,:,Lses)),3);     
end