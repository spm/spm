function M1 = spm_eeg_inv_headcoordinates(nas, lpa, rpa)
% Returns the homogeneous coordinate transformation matrix
% that converts the specified fiducials in any coordinate system (e.g. MRI)
% into the rotated and translated headccordinate system.
%
% FORMAT M1 = spm_eeg_inv_headcoordinates(nas, lpa, rpa)
%
% The headcoordinate system in CTF is defined as follows:
% the origin is exactly between lpa and rpa
% the X-axis goes towards nas
% the Y-axis goes approximately towards lpa, orthogonal to X and in the plane spanned by the fiducials
% the Z-axis goes approximately towards the vertex, orthogonal to X and Y
%__________________________________________________________________________

% Robert Oostenveld
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% ensure that they are row vectors
lpa = lpa(:)';
rpa = rpa(:)';
nas = nas(:)';

% compute the origin and direction of the coordinate axes in MRI coordinates

% follow CTF convention
origin = [lpa+rpa]/2;
dirx = nas-origin;
dirx = dirx/norm(dirx);
dirz = cross(dirx,lpa-rpa);
dirz = dirz/norm(dirz);
diry = cross(dirz,dirx);

% compute the rotation matrix
rot = eye(4);
rot(1:3,1:3) = inv(eye(3) / [dirx; diry; dirz]);
% compute the translation matrix
tra = eye(4);
tra(1:4,4)   = [-origin(:); 1];
% compute the full homogeneous transformation matrix from these two
M1 = rot * tra;
