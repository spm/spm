function spm_demo_proceed(tag,str)
% prompt for OK and activate correct figure
% FORMAT spm_demo_proceed(tag,str)
%
% tag - graphics tag
% str - string for dialogue box
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_demo_proceed.m 1380 2008-04-11 18:55:18Z karl $

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