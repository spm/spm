function surf = spm_eeg_inv_TesBin(n,ctr_vol,P,info);

%=======================================================================
% FORMAT ts = spm_eeg_inv_TesBin(n,ctr_vol,P,info);
%
% Generate a mesh covering a binary volume
%   1. Generate a simple spherical mesh
%   2. The spherical mesh is projected radially on the bin volume
%
% Then, "elastic" mesh correction is performed to correct for overlong edges.
%
% Input : 
% n          - number of vertices on each surface n = Npts^2*5/4+2
% ctr_vol    - centre of bin volume for the radial projection (mm)
% P          - filename of bin volume
% info       - information string
%
% Output :
% ts         - tesselation structure
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips & Jeremie Mattout
% $Id: spm_eeg_inv_TesBin.m 716 2007-01-16 21:13:50Z karl $

if nargin > 4
    error('Wrong input arguments for ''TesBin''.') ;
end
    
% Load volume information
%--------------------------------------------------------------------------
Vv       = spm_vol(P);
VOX      = sqrt(sum(Vv.mat(1:3,1:3).^2));
    
% A few definitions
%--------------------------------------------------------------------------
ho       = 1 ;            % Trilinear interpolation shoud be enough
n_div    = 1000;          % # of sample in the radial direction
trsh_vol = .9;            % Thershold for surface detection
rad      = round(4/3*max(Vv.dim(1:3).*VOX/2)); % Radius (mm) of original sphere
dr       = rad/n_div;     % size of sampling step (mm)
d_li     = 0:dr:rad;      % Sampling location on radius (mm)
nd_li    = length(d_li); unit = ones(1,nd_li);

% Create a tessalated sphere
%--------------------------------------------------------------------------
tsph = spm_eeg_inv_TesSph(rad,n);             % tesselated sphere of radius rad,
vert = [tsph.vert ; ones(1,tsph.nr(1))];      % centered around [0 0 0]
vert = spm_matrix([ctr_vol 0 pi/2 0])*vert;   % Rotate the sphere by 90deg around y axis and centre sphere around the "centre" of the brain,     
vert = vert(1:3,:) ;                          % All this in mm. Leave them as a 3xNvert matrix.

srf_vert = zeros(3,tsph.nr(1)) ; % vertices at the surface of brain, in vx !
spm_progress_bar('Init',tsph.nr(1),	['Generate ',info ],'Vertices projected');

for i = 1:tsph.nr(1)
	or       = ctr_vol'-vert(:,i) ; or = or/norm(or) ; % direction from the point toward the centre
	line     = vert(:,i)*unit + or*d_li;
    line_vx  = spm_eeg_inv_mm2vx(line,Vv.mat);
    val_line = spm_sample_vol(Vv,line_vx(1,:),line_vx(2,:),line_vx(3,:),ho);
	srf_vert(:,i) = line_vx(:,min(find(val_line>trsh_vol))); % first point to intercept the surface
   	spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


% tesllation structure
%--------------------------------------------------------------------------
srf_vert_mm      = spm_eeg_inv_vx2mm(srf_vert,Vv.mat);
surf.XYZmm       = srf_vert_mm(1:3,:);
surf.tri         = tsph.tri;
surf.nr          = tsph.nr;
surf.M           = Vv.mat;
surf.info.str    = info;
surf.info.bin_fn = P;
surf.Centre      = ctr_vol;

return
