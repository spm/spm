function [pt_mm] = spm_eeg_inv_vx2mm(pt_vx,M);

%=======================================================================
% FORMAT [pt_mm] = spm_eeg_inv_vx2mm(pt_vx,M);
% Transforms point(s) 'pt_vx' voxel coordinates into mm coordinates 'pt_mm'
% according to the transformation matrix 'M'.
% 'M' is the 4x4 affine transformation matrix : from vx to mm
%=======================================================================

if size(pt_vx,1)==3
    Npt = size(pt_vx,2);
elseif size(pt_vx,2)==3;
    pt_vx = pt_vx';
    Npt = size(pt_vx,2);
else
    error('Wrong vector format')
end
	
pt_mm = M*[pt_vx;ones(1,Npt)];
pt_mm = pt_mm(1:3,:);
