function spm_write_flow(p,vox,bb)
% Writes out flow-fields for use in diffeomorphic registration
% FORMAT u = spm_write_flow(p,vox,bb)
%     p   - parameters from spm_preproc.m
%     vox - voxel sizes
%     bb  - bounding box
%
%_______________________________________________________________________
% Copyright (C) 2006 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$

if nargin<3, bb  = NaN*ones(2,3); end;
if nargin<2, vox = NaN*ones(1,3); end;
if nargin<1, p   = spm_select(Inf,'.*seg_sn.mat','Select parameter files'); end;


if iscell(p),
    % Cellstring
    for i=1:numel(p),
        write_flowfield(load(deblank(p{i})),vox,bb);
    end;
elseif ischar(p),
    % Matrix of strings
    for i=1:size(p,1),
        write_flowfield(load(deblank(p(i,:))),vox,bb);
    end;
elseif isstruct(p),
    for i=1:numel(p),
        write_flowfield(p(i),vox,bb);
    end;
else
    error('Unable to do anything with this"');
end;

function write_flowfield(sn,vox,bb)
[Dsp,mat] = sn2disp(sn,vox,bb);
dm        = size(Dsp{1});
dat       = file_array;
[pth,nam,ext] = fileparts(sn.VF.fname);
pth = '.';
dat.fname  = fullfile(pth,['u_' nam '.nii']);
dat.dim    = [dm 1 3];
dat.dtype  = 'float32-le';
dat.scl_slope = 1;
dat.scl_inter = 0;
NU         = nifti;
NU.dat     = dat;
NU.mat     = mat;
NU.mat0    = mat;
NU.descrip = 'Velocity Field';
create(NU);
NU.dat(:,:,:,1,1) = Dsp{1};
NU.dat(:,:,:,1,2) = Dsp{2};
NU.dat(:,:,:,1,3) = Dsp{3};
return;
%=======================================================================

%=======================================================================
function [Dsp,mat] = sn2disp(sn,vox,bb)
% Convert a SPM _sn.mat file into a displacement field, and return it.

[bb0,vox0] = bbvox_from_V(sn.VG(1));

if any(~finite(vox)), vox = vox0; end;
if any(~finite(bb)),  bb  = bb0;  end;
bb  = sort(bb);
vox = abs(vox);

% Adjust bounding box slightly - so it rounds to closest voxel.
bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
bb(:,3) = round(bb(:,3)/vox(3))*vox(3);

M   = sn.VG(1).mat;
vxg = sqrt(sum(M(1:3,1:3).^2));
ogn = M\[0 0 0 1]';
ogn = ogn(1:3)';

% Convert range into range of voxels within template image
x   = (bb(1,1):vox(1):bb(2,1))/vxg(1) + ogn(1);
y   = (bb(1,2):vox(2):bb(2,2))/vxg(2) + ogn(2);
z   = (bb(1,3):vox(3):bb(2,3))/vxg(3) + ogn(3);

og  = -vxg.*ogn;
of  = -vox.*(round(-bb(1,:)./vox)+1);
M1  = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
M2  = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
mat = sn.VG(1).mat*inv(M1)*M2;

st = size(sn.Tr);

Dsp = single(0);
Dsp(numel(x),numel(y),numel(z)) = 0;
Dsp = {Dsp; Dsp; Dsp};

if (prod(st) ~= 0),
    basX = spm_dctmtx(sn.VG(1).dim(1),st(1),x-1);
    basY = spm_dctmtx(sn.VG(1).dim(2),st(2),y-1);
    basZ = spm_dctmtx(sn.VG(1).dim(3),st(3),z-1);
    for j=1:length(z)
        tx = reshape( reshape(sn.Tr(:,:,:,1),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
        ty = reshape( reshape(sn.Tr(:,:,:,2),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
        tz = reshape( reshape(sn.Tr(:,:,:,3),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );

        X1 = -basX*tx*basY';
        Y1 = -basX*ty*basY';
        Z1 = -basX*tz*basY';

        Dsp{1}(:,:,j) = single(X1);
        Dsp{2}(:,:,j) = single(Y1);
        Dsp{3}(:,:,j) = single(Z1);
    end;
end;


%_______________________________________________________________________
function [bb,vx] = bbvox_from_V(V)
% Return the default bounding box for an image volume

vx = sqrt(sum(V.mat(1:3,1:3).^2));
o  = V.mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V.dim(1:3)-o)];
return;
%_______________________________________________________________________

