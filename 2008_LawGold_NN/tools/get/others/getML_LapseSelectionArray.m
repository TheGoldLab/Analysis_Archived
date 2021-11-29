function L = getML_LapseSelectionArray(monk, type, thres)
% get a selection array for sessions with performance at 99.9% coh
% higher than 90% correct

if nargin<2
    type  = 'psy';
end

if nargin<3
    thres = 0.9;
end



%fits = getML_psyPerformanceCT([monk 'TRain_psy.txt'], 0);
%L    = fits(4,:)<=thres;
p = getML_performanceAsso([monk 'TRain_psy.txt'], 0)';
L = p>=thres;

if strcmp(lower(type), 'psy')
    return;
else
    % if type is 'MT', or 'LIP', get selection array for neural data
    %  remove neural data in sessions with behavioral lapse rate smaller
    %  than or equal to thres
    a     = getML_txt([monk, 'TRain_psy.txt']);
    ses_b = a.data{strcmp(a.name,'session')};

    a     = getML_txt([monk, 'TRain_' type '.txt']);
    ses_n = a.data{strcmp(a.name,'session')};

    L = ismember(ses_n, ses_b(L));
end


%% old way, use set threshold with time constant
% % get time constant
% [p n pci]    = getML_performanceAsso([monk 'TRain_asso.txt'], 0);
% [pt nt ptci] = getML_performanceAsso([monk 'TRain_psy.txt'], 0);
% N            = length(p);
% 
% lx   = [-N:1:-1 1:size(pt,2)];
% ly   = 1-[p pt];
% lyci = 1-[pci ptci];
% lys  = mean(abs([lyci(1,:)-ly; lyci(2,:)-ly]));
% lys(lys==0) = 0.001;   % make it a small number
% L    = ~isnan(lx) & ~isnan(ly) & ~isnan(lys);
% b = exp_fit2W(lx(L), ly(L), lys(L), [0.5; 5], [0.4 0.6;0.01 100]);
% 
% 
% a     = getML_txt([monk, 'TRain_psy.txt']);
% ses_b = a.data{strcmp(a.name,'session')};
% L     = ses_b>=round(1.6*b(2))-N;
% 
% 
% if strcmp(lower(type), 'psy')
%     return;
% else
%     % if type is 'MT', or 'LIP', get selection array for neural data
%     %  remove neural data in sessions with behavioral lapse rate smaller
%     %  than or equal to thres
% 
%     a     = getML_txt([monk, 'TRain_' type '.txt']);
%     ses_n = a.data{strcmp(a.name,'session')};
%     L = ismember(ses_n, ses_b(L));
% end