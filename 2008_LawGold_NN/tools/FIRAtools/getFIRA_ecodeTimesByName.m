function times = getFIRA_ecodeTimesByName(ec, offset, trials)
% function varargout = gsGUI_ecodeTimesByName(key, varargin)
%
%   assumes there is a GUI with a menu and an edit box
%   'set' makes the menu entries the "time" ecode fields from FIRA
%   'get' gets an array of times for each of the given trials
%       of the time of the event chosen in the menu plus the offset in
%       the text box
%
% Usage:
%              gsGUI_ecodeTimesByName('setf', menu_h, default_index, ...
%                                       edit_h, edit_value);
%      times = gsGUI_ecodeTimesByName('getf', menu_h, edit_h, ...
%                                       trials);

% modified from gsGUI_ecodeTimesByName.m (created by jig) by ctl on 8/10/05



% get the times
if nargin==0 | nargin>3
    return
elseif nargin<3
    times = getFIRA_ecodesByName(ec, 'time');
else
    times = getFIRA_ecodesByName(ec, 'time', trials);
end

% add offset, if given
if ~isempty(offset)
    times = times+offset;
end