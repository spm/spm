function spm_demo_proceed(tag,str)
% prompt for OK and activate correct figure
% FORMAT spm_demo_proceed(tag,str)
%
% tag - graphics tag
% str - string for dialogue box
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_demo_proceed.m 4136 2010-12-09 22:22:28Z guillaume $

% get figure
%--------------------------------------------------------------------------
try, tag; catch, tag = 'MFM';         end
try, str; catch, str = [tag ' demo']; end

% get figure
%--------------------------------------------------------------------------
drawnow
uiwait(warndlg(str,'Proceed with demonstration?'));
spm_figure('GetWin',tag);
clf