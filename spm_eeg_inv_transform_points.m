function new = spm_eeg_inv_transform_points(M, old)
% Applies homogenous transformation to a set of 3D points
% FORMAT new = spm_eeg_inv_transform_points(M, old)
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_transform_points.m 2863 2009-03-11 20:25:33Z guillaume $

old(:,4) = 1;
new = old * M';
new = new(:,1:3);
