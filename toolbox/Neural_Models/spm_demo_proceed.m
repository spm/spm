function spm_demo_proceed(tag,str)
% prompt for OK and activate correct figure
% FORMAT spm_demo_proceed(tag,str)
%
% tag - graphics tag
% str - string for dialogue box
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

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
