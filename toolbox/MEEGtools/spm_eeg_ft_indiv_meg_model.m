function  D = spm_eeg_ft_indiv_meg_model(varargin)
% Function for making an individually fitted MEG head model based on
% subject's sMRI.
%
% FORMAT  D = spm_eeg_ft_indiv_meg_model(D, val, Msize, sMRI, hs_source, vm_source)
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL

% Robert Oostenveld, Vladimir Litvak,  
% $Id: spm_eeg_ft_indiv_meg_model.m 2419 2008-10-30 19:40:32Z vladimir $

[Finter] = spm('FnUIsetup','FT based MEG head model',0);

[D,val] = spm_eeg_inv_check(varargin{:});

if val == 0
    val = 1;
end

if ~isfield(D, 'inv')
    D.inv = {struct([])};
end

if ~isfield(D.inv{val}, 'modality')
    modality = spm_eeg_modality_ui(D, 1);
    if ~strcmp(modality, 'MEG')
        error('Only MEG is supported');
    end
    D.inv{val}(1).modality = modality;
end

if nargin>2
    Msize = varargin{3};
else
    Msize = spm_input('Mesh size (vertices)', '+1','3000|4000|5000|7200', [1 2 3 4]);
end

if nargin>3
    sMRI = varargin{4};
else
    sMRI =  spm_select(1,'image', 'Select the subject''s structural image');
end

% Segment structural if necessary
%======================================================================
[spmvol, spmfid, mesh] = spm_eeg_inv_meshing(sMRI, Msize, 'MEG');
%%
[p f] = fileparts(strtok(mesh.sMRI, ',')); x = '.nii';

seg = [];
for c = 1:5
    mri = ft_read_fcdc_mri(fullfile(p, ['c' num2str(c) f x]));
    %mri.anatomy = flipdim(mri.anatomy,  1);
    if c == 1
        anatomy =  mri.anatomy;
    else
        anatomy = anatomy+mri.anatomy;
    end
    switch c
        case 1
            seg.gray  = mri.anatomy;
            seg.transform = mri.transform;
            seg.dim = mri.dim;
        case 2
            seg.white  = mri.anatomy;
        case 3
            seg.csf = mri.anatomy;
    end
end

mri.anatomy = anatomy;

cfg = [];
cfg.spheremesh = 1000;
ftvol = ft_prepare_singleshell(cfg, seg);

ftvol = forwinv_convert_units(ftvol, 'mm');
%%
cfg = [];
cfg.downsample = 2;
mri = ft_volumedownsample(cfg, mri);

seg0 = mri.anatomy;

seg1 = mri.anatomy > 0.9;

bol3 = strel_bol(3);
bol5 = strel_bol(5);
bol7 = strel_bol(7);

seg2 = convn(seg1, bol3, 'same');

seg3 = seg2>40;

iml = bwlabeln(seg3);

seg4 = (iml==mode(iml(iml(:)>0)));

seg4(:,:,(end-0)) = 1;
seg4(:,:,(end-1)) = 1;
seg4(:,:,(end-2)) = 1;

seg5 = imfill(seg4, 'holes');
seg5(:,:,(end-0)) = 0;
seg5(:,:,(end-1)) = 0;
seg5(:,:,(end-2)) = 0;

seg6 = imerode(seg5, bol5);
iml = bwlabeln(seg6);
seg7 = imdilate(iml==1, bol5);

% make the surface
seg = seg7;
npnt = 1000;
[mrix, mriy, mriz] = ndgrid(1:size(seg,1), 1:size(seg,2), 1:size(seg,3));

ori(1) = mean(mrix(seg(:)));
ori(2) = mean(mriy(seg(:)));
ori(3) = mean(mriz(seg(:)));
[pnt, tri] = triangulate_seg(seg, npnt, ori);
% apply the coordinate transformation from voxel to head coordinates
pnt(:,4) = 1;
pnt = (mri.transform * (pnt'))';
pnt = pnt(:,1:3);

ftfid = [];
ftfid.pnt = pnt;
ftfid.tri = tri;
ftfid.fid = spmfid.fid;
%%
D.inv{val}.mesh = mesh;

if nargin>4
    hs_source = varargin{5};
else
    hs_source = spm_input('Which head surface?', '+1','SPM|FT', [1 0], 0);
end

if hs_source
    D.inv{val}.datareg.fid_mri = spmfid;
else
    D.inv{val}.datareg.fid_mri = ftfid;
end

if nargin>4
    vm_source = varargin{5};
else
    vm_source = spm_input('Which volume model?', '+1','SPM|FT', [1 0], 1);
end

if vm_source
    D.inv{val}.forward.vol = spmvol;
else
    D.inv{val}.forward.vol = ftvol;
end

% check meshes and display
%--------------------------------------------------------------------------
spm_eeg_inv_checkmeshes(D);

save(D);

%%
function se = strel_bol(r)

% STREL_BOL constructs a 3D sphere with the specified radius
% that can be used as structural element in 3D image processing
%
% See STREL, IMERODE, IMDILATE (image processing toolbox)

dim = [2*r+1, 2*r+1, 2*r+1];
se  = zeros(dim);
for i=1:dim(1)
  for j=1:dim(2)
    for k=1:dim(3)
      x = i-1-r;
      y = j-1-r;
      z = k-1-r;
      if norm([x y z])<=r
        se(i,j,k) = 1;
      end
    end
  end
end


function [pnt, tri] = triangulate_seg(seg, npnt, ori);

% TRIANGULATE_SEG constructs a triangulation of the outer surface of a
% segmented volume. It starts at the center of the volume and projects the
% vertices of an evenly triangulated sphere onto the outer surface. The
% resulting surface is star-shaped from the center of the volume.
%
% Use as
%   [pnt, tri] = triangulate_seg(seg, npnt)
%
% See also KSPHERE

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: triangulate_seg.m,v $
% Revision 1.3  2006/07/26 11:05:58  roboos
% use find('last') for matlab 7 and higher, and regular find for older matlab versions
%
% Revision 1.2  2006/04/03 10:39:24  roboos
% added origin of the projection towards the surface as third (optional) input argument
%
% Revision 1.1  2005/11/03 11:12:32  roboos
% new implementation, using a projection of a ksphere triangulation from the center of the volume
%

seg = (seg~=0);
dim = size(seg);
len = ceil(sqrt(sum(dim.^2))/2);

if nargin<3
  ori(1) = dim(1)/2;
  ori(2) = dim(2)/2;
  ori(3) = dim(3)/2;
end

% start with a unit sphere with evenly distributed vertices
[pnt, tri] = ksphere(npnt);

for i=1:npnt
  % construct a sampled line from the center of the volume outward into the direction of the vertex
  lin = (0:0.5:len)' * pnt(i,:);
  lin(:,1) = lin(:,1) + ori(1);
  lin(:,2) = lin(:,2) + ori(2);
  lin(:,3) = lin(:,3) + ori(3);
  % round the sampled line towards the nearest voxel indices, which allows
  % a quick nearest-neighbour interpolation/lookup
  lin = round(lin);
  % exclude indices that do not ly within the volume
  sel = lin(:,1)<1 | lin(:,1)>dim(1) | ...
        lin(:,2)<1 | lin(:,2)>dim(2) | ...
        lin(:,3)<1 | lin(:,3)>dim(3);
  lin = lin(~sel,:);
  sel = sub2ind(dim, lin(:,1), lin(:,2), lin(:,3));
  % interpolate the segmented volume along the sampled line
  int = seg(sel);
  % find the last sample along the line that is part of the segmentation
  try
    % for matlab 7 and higher
    sel = find(int, 1, 'last');
  catch
    % for older matlab versions
    sel = find(int);
    sel = sel(end);
  end
  pnt(i,:) = lin(sel,:);
end

% undo the shift of the origin from where the projection is done
% pnt(:,1) = pnt(:,1) - ori(1);
% pnt(:,2) = pnt(:,2) - ori(2);
% pnt(:,3) = pnt(:,3) - ori(3);

% fast unconditional re-implementation of the standard Matlab function
function [s] = sub2ind(dim, i, j, k)
s = i + (j-1)*dim(1) + (k-1)*dim(1)*dim(2);

function [pnt, tri] = ksphere(N);

% KSPHERE returns a triangulated sphere with K vertices that are
% approximately evenly distributed over the sphere. 
%
% Use as
%   [pnt, tri] = ksphere(K);
%
% This implements the recommended algorithm in spherical coordinates
% (theta, phi) according to "Distributing many points on a sphere"
% by E.B. Saff and A.B.J. Kuijlaars, Mathematical Intelligencer 19.1
% (1997) 5--11
%
% See also http://www.math.niu.edu/~rusin/known-math/97/spherefaq

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: ksphere.m,v $
% Revision 1.1  2005/11/01 09:56:00  roboos
% new implementation
%

for k=1:N
  h = -1 + 2*(k-1)/(N-1);
  theta(k) = acos(h);
  if k==1 || k==N
    phi(k) = 0;
  else
    phi(k) = mod((phi(k-1) + 3.6/sqrt(N*(1-h^2))), (2*pi));
  end
end
az = phi;
el = theta - pi/2;
[x, y, z] = sph2cart(az', el', 1);
pnt = [x, y, z];
tri = convhulln(pnt);