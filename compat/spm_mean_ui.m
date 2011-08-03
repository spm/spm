function spm_mean_ui
% Prompt for a series of images and averages them
%__________________________________________________________________________
%
% This function is deprecated. Use spm_mean instead.
%__________________________________________________________________________
% Copyright (C) 1998-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner, Andrew Holmes
% $Id: spm_mean_ui.m 4419 2011-08-03 18:42:35Z guillaume $

persistent runonce
if isempty(runonce)
    warning('spm_mean_ui is deprecated. Use spm_mean instead.');
    runonce = 1;
end

spm_mean;
