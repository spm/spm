function [pt_vx] = spm_eeg_inv_mm2vx(pt_mm,M)
% Apply a transformation (mm -> vox) to a set of points
% FORMAT [pt_vx] = spm_eeg_inv_mm2vx(pt_mm,M);
% Transforms point(s) 'pt_mm' mm coordinates into voxel coordinates 'pt_vx'
% according to the transformation matrix 'M'.
% 'M' is the 4x4 affine transformation matrix : from vx to mm
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_mm2vx.m 2696 2009-02-05 20:29:48Z guillaume $

if size(pt_mm,1)==3
    Npt   = size(pt_mm,2);
elseif size(pt_mm,2)==3;
    pt_mm = pt_mm';
    Npt   = size(pt_mm,2);
else
    error('Wrong format')
end
    
pt_vx = M\[pt_mm; ones(1,Npt)];
pt_vx = pt_vx(1:3,:);

