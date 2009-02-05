function new = spm_eeg_inv_transform_points(M, old)
% Applies homogenous transformation to a set of 3D points
% FORMAT new = spm_eeg_inv_transform_points(M, old)
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_transform_points.m 2696 2009-02-05 20:29:48Z guillaume $

old(:,4) = 1;
new = old * M';
new = new(:,1:3);