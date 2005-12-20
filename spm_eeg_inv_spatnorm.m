function varargout = spm_eeg_inv_spatnorm(varargin)

%=======================================================================
% Spatial Normalization (using VBM preprocessing - SPM5)
% transforms individual sMRI into MNI T1 space
% save the inverse deformations (...vbm_inv_sn.mat)
% that will be needed for computing the individual mesh
%
% FORMAT D = spm_eeg_inv_spatnorm(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the inverse deformation .mat file
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_spatnorm.m 389 2005-12-20 12:17:18Z jeremie $

spm_defaults

load('defaults_eeg_mesh.mat');

if nargin == 0
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
	D = spm_eeg_ldata(D);
elseif nargin == 1
    D = varargin{1};
else
	error(sprintf('Trouble reading the data file\n'));
    return
end

val = length(D.inv);

if isempty(D.inv{val}.mesh.sMRI)
    D.inv{val}.mesh.sMRI = spm_select(1,'image','Select subject T1 MRI');
end

jobs{1}.spatial{1}.preproc.data = D.inv{val}.mesh.sMRI;

% Spatial Transformation into MNI space
res = spm_preproc(jobs{1}.spatial{1}.preproc.data);
[sn,isn] = spm_prep2sn(res);

[pth,nam,ext] = spm_fileparts(jobs{1}.spatial{1}.preproc.data);

def_name = [nam '_vbm_sn_' num2str(val) '.mat'];
isndef_name = [nam '_vbm_inv_sn_' num2str(val) '.mat'];
D.inv{val}.mesh.def = fullfile(pth,def_name);
D.inv{val}.mesh.invdef = fullfile(pth,isndef_name);
savefields(D.inv{val}.mesh.def,sn);
savefields(D.inv{val}.mesh.invdef,isn);

% Writing the images
spm_preproc_write(sn,jobs{1}.spatial{1}.preproc.output);
nobias_name = ['m' nam ext];
D.inv{val}.mesh.nobias = fullfile(pth,nobias_name);

% Downsampling the segmented images (by a factor 2 in each direction)
for i = 1:3
    PI = fullfile(pth,['c' num2str(i) nam ext]);
    PO = fullfile(pth,['dc' num2str(i) nam ext]);
    VI = spm_vol(PI);
    dim = round(VI.dim/2);
    mat = [2*VI.mat(:,1:3) VI.mat(:,4)];
    reslice(PI,PO,dim,mat,1);
    [success] = copyfile(PO,PI);
    if success
        delete(PO);
    else
        error(sprintf('Impossible to write file downsampled segmented volume\n'));
    end
    if strcmp(ext,'.img')
        HI = fullfile(pth,['c' num2str(i) nam '.hdr']);
        HO = fullfile(pth,['dc' num2str(i) nam '.hdr']);
        [success] = copyfile(HO,HI);
        if success
            delete(HO);
        else
            error(sprintf('Impossible to write file header of downsampled segmented volume\n'));
        end
    end
end
PI = D.inv{val}.mesh.nobias;
PO = fullfile(pth,['temp_downsampled_mask' ext]);
VI = spm_vol(PI);
dim = round(VI.dim/2);
mat = [2*VI.mat(:,1:3) VI.mat(:,4)];
reslice(PI,PO,dim,mat,1);
[success] = copyfile(PO,PI);
if success
    delete(PO);
else
    error(sprintf('Impossible to write file downsampled segmented volume\n'));
end
[pth,nam,ext] = spm_fileparts(PI);
if strcmp(ext,'.img')
    HI = fullfile(pth,[nam '.hdr']);
    HO = fullfile(pth,'temp_downsampled_mask.hdr');
    [success] = copyfile(HO,HI);
    if success
        delete(HO);
    else
        error(sprintf('Impossible to write file header of downsampled segmented volume\n'));
    end
end

clear jobs

save(D.fname, 'D');
varargout{1} = D;

spm('Pointer','Arrow');
D.inv{val}.mesh

%=======================================================================
function savefields(fnam,p)
if length(p)>1
    error('Can''t save fields.')
end
fn = fieldnames(p);
for i=1:length(fn)
    eval([fn{i} '= p.' fn{i} ';']);
end
if str2num(version('-release'))>=14
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end
%=======================================================================


%=======================================================================
function reslice(PI,PO,dim,mat,hld)
% FORMAT reslice(PI,PO,dim,mat,hld)
%   PI - input filename
%   PO - output filename
%   dim - 1x3 matrix of image dimensions
%   mat - 4x4 affine transformation matrix mapping
%         from vox to mm (for output image).
%         To define M from vox and origin, then
%             off = -vox.*origin;
%              M   = [vox(1) 0      0      off(1)
%                     0      vox(2) 0      off(2)
%                     0      0      vox(3) off(3)
%                     0      0      0      1];
%
%   hld - interpolation method.
%___________________________________________________________________________
% %W% John Ashburner %E%

VI          = spm_vol(PI);
VO          = VI;
VO.fname    = deblank(PO);
VO.mat      = mat;
VO.dim(1:3) = dim;

VO = spm_create_vol(VO); % changed from spm_create_image, KCR
for x3 = 1:VO.dim(3)
    M  = inv(spm_matrix([0 0 -x3 0 0 0 1 1 1])*inv(VO.mat)*VI.mat);
    v  = spm_slice_vol(VI,M,VO.dim(1:2),hld);
    VO = spm_write_plane(VO,v,x3);
end
%=======================================================================
