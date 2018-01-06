function fcs_ = getBIAS_fcs(monks, Larrays)
%
% 
% get filtered choices

if nargin < 1 || isempty(monks)
    monks = getBIAS_monks;
elseif ischar(monks)
    monks = {monks};
end

monklist = getBIAS_monks;
wkdat    = FS_loadProjectFile('2008_Bias', 'figBIAS_wkFits');

%wkdat    = loadBIAS_figData('figBIAS_wkFits');
fcs_     = cell(length(monks), 1);
for mm = 1:length(monks)
    if nargin > 1 && ~isempty(Larrays)
        fcs_{mm} = wkdat{strcmp(monks{mm}, monklist),3}(Larrays(:,mm),1);
    else
        fcs_{mm} = wkdat{strcmp(monks{mm}, monklist),3}(:,1);
    end    
end

if length(fcs_) == 1
    fcs_ = fcs_{1};
end