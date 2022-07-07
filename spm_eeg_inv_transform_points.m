function new = spm_eeg_inv_transform_points(M, old)
% Apply homogeneous transformation to a set of 3D points
% FORMAT new = spm_eeg_inv_transform_points(M, old)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


old(:,4) = 1;
new = old * M';
new = new(:,1:3);
