function new = spm_eeg_inv_transform_points(M, old)
% new = spm_eeg_inv_transform_points(M, old)
% Applies homogenous transformation to cortical mesh. 
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_transform_points.m 2389 2008-10-23 11:15:23Z vladimir $

old(:,4) = 1;
new = old * M';
new = new(:,1:3);