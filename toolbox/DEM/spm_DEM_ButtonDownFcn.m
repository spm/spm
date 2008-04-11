function spm_DEM_ButtonDownFcn
% ButtonDownFcn to play a sound on button press
% FORMAT spm_DEM_ButtonDownFcn
%
% Require gcbo to have UserData.Y and Userdata.FS
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_ButtonDownFcn.m 1380 2008-04-11 18:55:18Z karl $

% default
%--------------------------------------------------------------------------
S = get(gcbo,'Userdata');
sound(S{1},S{2});

