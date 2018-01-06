function [axs_,fig_] = getBIAS_axes(num, wid, hts, cols, psh, psw, fs, al)
% function [axs_,fig_] = getBIAS_axes(num, wid, hts, cols, psh, psw, fs, al)
%
% panel separation height,width

% Figure number
if nargin < 1 || isempty(num)
    num = 1;
end

% total width
if nargin < 2 || isempty(wid)
    wid = 6.0;
end

% heights of each column
if nargin < 3 || isempty(hts)
    hts = 1.5;
end

% number of columns in each row
if nargin < 4 || isempty(cols)
    cols = {1};
end

% panel separation height
if nargin < 5 || isempty(psh)
    psh = 0.5;
end

% panel separation width
if nargin < 6 || isempty(psw)
    psw = 0.5;
end

% font size
if nargin < 7 || isempty(fs)
    fs = 12;
end

% author list
if nargin < 8 || isempty(al)
    al = 'Gold et al';
end

% get figure
fig_ = figure;
wysifig(fig_);

% add label
ax = axes('Units', 'inches', 'position', [7 10 .2 .2]);
set(ax,'Visible','off');
h=[text(0,1,al), text(0,0,sprintf('Figure %d',num))];
set(h,'FontSize',14);

% useful variables
if length(hts) > 1
    nrow = length(hts);
elseif length(cols) > 1
    nrow = length(cols);
    hts  = hts*ones(1,nrow);
else
    nrow = 1;
end

lt   = (8.5 - wid)/2;
bot  = 9.5 - sum(hts) - psh*(nrow-1);
bts  = cumsum([bot flipdim(hts,2)+psh]);
bts  = bts(end-1:-1:1);

axi   = 1;
axs_  = [];

% loop through the rows
for rr = 1:nrow

    % if cols{rr} is a scalar, make that many equally
    %   spaced columns
    if isscalar(cols{rr}) && cols{rr} >= 1
        cols{rr} = 1/cols{rr}.*ones(1,cols{rr});
    end

    % compute width, left position per panel (column)
    ncol = length(cols{rr});
    if ncol >= 1
        wids = cols{rr}*(wid-psw*(ncol-1));
        lts  = lt+cumsum([0 psw+wids(1:end-1)]);
    end

    % loop through the columns
    for cc = 1:length(cols{rr})
        axs_(axi) = axes('Units', 'inches', 'position', ...
            [lts(cc), bts(rr), wids(cc), hts(rr)]);
        axi = axi+1;
    end
end

set(axs_, 'FontSize', fs);