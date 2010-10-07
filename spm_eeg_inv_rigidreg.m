function [M1] = spm_eeg_inv_rigidreg(data1, data2)
% Computes homogenious transformation matrix based on two sets
% of points from two coordinate systems
%
% FORMAT [M1] = spm_eeg_inv_rigidreg(data1, data2)
% Input:
% data1      - locations of the first set of points corresponding to the
%            3D surface to register onto 
% data2      - locations of the second set of points corresponding to the
%            second 3D surface to be registered 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_rigidreg.m 4082 2010-10-07 16:11:48Z guillaume $
 
M       = spm_detrend(data1');
S       = spm_detrend(data2');
[U A V] = svd(S'*M);
R1      = V*U';
if det(R1) < 0
    B      = eye(3);
    B(3,3) = det(V*U');
    R1     = V*B*U';
end
t1 = mean(data1,2) - R1*mean(data2,2);
M1 = [R1 t1; 0 0 0 1];
