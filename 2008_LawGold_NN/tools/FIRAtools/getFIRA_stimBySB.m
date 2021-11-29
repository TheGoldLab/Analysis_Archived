function stim = getFIRA_stimBySB(aperture, direction, speed, coherence, duration, sb)
%function [stim, stim_rot, stim_cart] = getFIRA_stimBySB(aperture, direction, speed, coherence, duration, sb)
%
% This function takes parameters of the dots stimulus and return location
% of the dots at different time, computed using the seed base.
%
% USAGE:
%   INPUTS:
%       aperture  - dots aperture (deg * 10)
%       direction - dot direction (angle)
%       speed     - dot speed     (deg/sec^2 * 10)
%       coherence - dot coherence (0-1000)
%       duration  - dot duration  (second)
%       sb        - seed base
%
%   OUTPUTS:
%       stim      - a 2-d matrix of the actual stimulus. Each row represent
%                   the location of the dots in pixel(projected along the direction
%                   axis) at each screen frame. The default screen reflesh
%                   rate is 85.
%       
%
% modify from dots8.m to get stimulus from seed base

% created by jcl on 11/15/06

%%%                                       %%%
%       Check inputs and set parameters     %
%%%                                       %%%

% default parameters
MON_REFRESH			= 85;	% frames per second
MAX_DOTS_PER_FRAME 	= 150;	% by trial and error.  Depends on graphics card
mon_horizontal_cm  	= 40;	% horizontal dimension (cm) of SCREEN_RECT
view_dist_cm 		= 60;	% IMPORTANT to set this (cm) VIEW DISTANCE
hori_pix            = 1024;
PIX_PER_DEG         = hori_pix * (1 ./ (2 * atan2(mon_horizontal_cm / 2, view_dist_cm))) * pi/180;

% stimulus parameters
rseed   = [sb, floor(coherence/10)+(cos(direction*pi/180)>=0)];
AP_D    = aperture;
d_ppd 	= AP_D/10 * PIX_PER_DEG;	% size of aperture
AP_RECT = [0 0 d_ppd d_ppd];
ap_off  = [AP_RECT(1) AP_RECT(2)];  % aperture upper/left corner
coh   	= coherence/1000;           % DtCoh is on 0..999.  We want 0..1
ndots 	= max(3,min([MAX_DOTS_PER_FRAME, ceil(16.7 * AP_D * AP_D * 0.01 / MON_REFRESH)]));




%%%                                       %%%
%              Get dots position            %
%%%                                       %%%

% get dots position
rand('state', prod(rseed));
dxdy 	= (speed/10) * (10/AP_D) * (3/MON_REFRESH)...
* [cos(pi*direction/180.0) sin(pi*direction/180.0)];
dxdy 	= dxdy(ones(ndots,1),:);

% ARRAYS, INDICES for loop
nloops	= 3;	% number of loops
rs		= rand( ndots, 3, nloops); % array of dot positions raw; third col is history
as 		= zeros(ndots, 4, nloops); % array of dot positions within aperture
%as(:,3,:) = DOT_WHITE_INDEX; 	 % if we want to specify dot CLUT indices
as(:,4,:) = 3; 				 % if we want to specify dot sizes

% divide dots into nl sets...
loopi     = 1; 	% loops through the three sets of dots
numframe  = floor(abs(duration)*MON_REFRESH);   % use floor function b/c if remaining duration is shorter than a frame, no dot will be shown on the next one
cont      = numframe;       
went_flag = 0;
wf        = 0;

% tmp array
pts      = zeros(ndots, 2);
past_pts = nans(ndots, 2, cont);	% keep track of last set of dots (to erase)

% THE MAIN LOOP
while cont
	
	% compute new locations
	L = rand(ndots,1) < coh;		    % coherence is specified on 0..1
	
	% NEW jig 12/27/04 -- limit lifetime
	num = sum(L);
	if(num)
		[y,I] 		= sort(rs(:,3,loopi));
		L	  		= isnan(L);
		L(I(1:num)) = logical(1);
	end
	% update lifetimes
	rs(L ,3,loopi) = rs(L,3,loopi)+1;
	rs(~L,3,loopi) = 0;
	
	% get current raw pts
	pts = rs(:,1:2,loopi);
	
	% offest selected dots
	pts(L,:) = pts(L,:) + dxdy(L,:);
	
	% get new random locations for the rest
	if any(~L)
		pts(~L,:) = rand(sum(~L),2);	
    end
 
	% wrap around -- jig changed 12/27/04 to real wrapping

	Llo = pts < 0;
	Lhi = pts > 1;
	Lany = any(Llo|Lhi,2);
	if any(Lany)
		pts(Llo) = 1.0-0.1*rand(sum(Llo(:)),1);
		pts(Lhi) =     0.1*rand(sum(Lhi(:)),1);
		pts(flipdim(Llo|Lhi,2)) = rand(sum(sum(Llo|Lhi)),1);
	end

	% re-save the raw pts
	rs(:,1:2,loopi) = pts;
	
	% old-style (random) wrap
	%	L = sum((rs(:,1:2,loopi) > 1 | rs(:,1:2,loopi) < 0)')' ~= 0;  
	%	if any(L)
	%		rs(L,1:2,loopi) = rand(sum(L),2);	
	%	end
	
	% convert from raw (0..1) to pixels
	as(:,1:2,loopi) = round(d_ppd*pts); % pix/ApUnit
	
	% save these pts as "past"
	past_pts(:,:,numframe-cont+1) = rs(:,1:2,loopi);
	
	% update the loop pointer
	loopi = loopi+1;
	if loopi > nloops
		loopi = 1;
	end

    % check for end of loop
    cont = cont - 1;
end




%%%                                       %%%
%            Process dots position          %
%%%                                       %%%

% remove dots outside of the unit circle (x-0.5)^2+(y-0.5)^2>0.5^2,
% this is equivalent to overlaying a circular mask on the stimulus
x = past_pts(:,1,:);
x = x(:);
y = past_pts(:,2,:);
y = y(:);
z = repmat([1:size(past_pts,3)], size(past_pts,1), 1);
z = z(:);

L = (x-0.5).^2+(y-0.5).^2>0.5^2;
x(L) = [];
y(L) = [];
z(L) = [];
stim_cart = [z,x,y];

% rotate stimulus space to the direction axis,
% so now x-axis will be the direction axis (Pref-Null axis)
% and y-axis will be orthogonal to the direction axis
if cos(pi/180*direction)>=0     % get unit vector for direction axis
    angle = direction;
    uvx  = [cos(pi/180*angle) sin(pi/180*angle)];
    uvy  = [-sin(pi/180*angle) cos(pi/180*angle)];
else
    angle = mod(direction+180,360);
    uvx  = [cos(pi/180*angle) sin(pi/180*angle)];
    uvy  = [-sin(pi/180*angle) cos(pi/180*angle)];
end
xp = (x-0.5)*uvx(1)+(y-0.5)*uvx(2)+0.5;  % projection of [x,y] along the direction axis
yp = (x-0.5)*uvy(1)+(y-0.5)*uvy(2)+0.5;  % projection of [x,y] along the orthogonal axis
stim_rot = [z,xp,yp];

% get spatio-temporal stimulus matrix, where y-axis is displacement along
% the direction axis, x-axis is time
% It can be used to compute the motion energy using fft2
xpix      = round(xp*round(d_ppd));     % location of dots in pixel, ranging from 0 to round(d_ppd)
ares      = round(d_ppd)+1;             % resolution of xp
tres      = numframe;           % resolution of time, or the frame rate
stim      = zeros(tres,ares);           % z,xp

% add 3-pixel wide dots
dsize = 3;
for i = 1:dsize
    Lbd       = (xpix+i)>ares;    % ignore part of that dot that is beyond the border
    I         = (xpix(~Lbd)+i-1)*tres+z(~Lbd);
    stim(I)   = 1;
end


