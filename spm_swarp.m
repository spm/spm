function that = spm_swarp(this,def,M)
% Warp surface.
% FORMAT that = spm_swarp(this,def)
% this - a gifti object
% def  - a deformation (nifti object or filename)
% that - the warped gifti object
%
% FORMAT that = spm_swarp(this,def,M)
% this - a gifti object
% def  - a deformation field (nx*ny*nz*1*3)
% M    - mapping from voxels to world, for deformation field
% that - the warped gifti object
%
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_swarp.m 2705 2009-02-06 15:37:19Z john $

if ~isa(this,'gifti'), this = gifti(this); end

if nargin<2, that = this; return; end

if ischar(def) || isa(def,'nifti'),
    if ischar(def),
        def = nifti(def);
    end
    y   = def(1).dat(:,:,:,:,:);
    M   = def(1).mat;
else
    y = def;
    if nargin<3, M = eye(4); end
end

%vrt = this.vertices;
vrt = subsref(this,substruct('.','vertices'));

iM  = inv(M);
vrt = iM(1:3,1:4)*[vrt'; ones(1,size(vrt,1))];
xyz = {double(vrt(1,:)'),double(vrt(2,:)'),double(vrt(3,:)')};
vrt = [spm_bsplins(y(:,:,:,1,1),xyz{:},[1 1 1 0 0 0]),...
       spm_bsplins(y(:,:,:,1,2),xyz{:},[1 1 1 0 0 0]),...
       spm_bsplins(y(:,:,:,1,3),xyz{:},[1 1 1 0 0 0])];

that = subsasgn(this,substruct('.','vertices'),vrt);

