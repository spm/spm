function varargout = spm_eeg_inv_getmasks(varargin);

%=======================================================================
% Generate the binary images (masks) of the inner skull and scalp surfaces
%
% FORMAT [D] = spm_eeg_inv_getmasks(S,flags);
% Input:
% S		    - input data struct (optional)
% flags     - useful flags
%               * img_norm: normalised image (1) or not (0)
%               * ne,ng   : # of iterations for erosion and growing
%               * thr_im  : threhold applied to the binary images
%           NOTE: these flags are optional
% Output:
% D			- same data struct including the new files and parameters
%
%
% FORMAT [FNbin,flags] = spm_eeg_inv_getmasks(Ivol,flags);
% Input :
% Ivol      - filename of the structural image
% flags     - useful flags
%               * img_norm: normalised image (1) or not (0)
%               * ne,ng   : # of iterations for erosion and growing
%               * thr_im  : threhold applied to the binary images
%           NOTE: these flags are optional
% Output :
% FNbin     - file names of the generated binary images, scalp and iskull
% flags     - flag values
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_getmasks.m 389 2005-12-20 12:17:18Z jeremie $

spm_defaults

def_flags = struct('img_norm',0,'ne',1,'ng',2,'thr_im',[.5 .05]);

% Input arguments
if nargout == 1 
    
    try
        D = varargin{1};
    catch
        D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
        D = spm_eeg_ldata(D);
    end
    
    if ~isfield(D,'inv')
        disp('Error: no inverse structure has been created for this data set');
        return
    end

    val = length(D.inv);

    if isempty(D.inv{val}.mesh.nobias)
        Ivol = spm_select(1,'image','Select nobias file');
        D.inv{val}.mesh.nobias = Ivol;
    else
        Ivol = D.inv{val}.mesh.nobias;
    end
    
elseif nargout == 2
    
    Ivol = varargin{1};
    
else
    
    disp('Error: wrong number of arguments');
    return;
    
end
    
if nargin < 2
    flags = def_flags;
else
    flags = varargin{2};
    fnms  = fieldnames(def_flags);
    for i=1:length(fnms),
		if ~isfield(flags,fnms{i})
            flags = setfield(flags,fnms{i},getfield(def_flags,fnms{i}));
        end
	end
end

fprintf('\n\n'), fprintf('%c','='*ones(1,80)), fprintf('\n')
fprintf(['\tGenerate binary volumes from structural image.\n']);

[pth,nam,ext] = spm_fileparts(D.inv{val}.mesh.sMRI);
if isempty(D.inv{val}.mesh.invdef)
    Misn = fullfile(pth,[nam '_vbm_inv_sn_1.mat']);
else
    Misn = D.inv{val}.mesh.invdef;
end
if isempty(D.inv{val}.mesh.def)
    Msn = fullfile(pth,[nam '_vbm_sn_1.mat']);
else
    Msn = D.inv{val}.mesh.def;
end

% Adds up GM/WM/CSF to produce the inner skull surface
Iin     = strvcat(fullfile(pth,['c1',nam,ext]), ...
                  fullfile(pth,['c2',nam,ext]), ...
                  fullfile(pth,['c3',nam,ext]));
Iisk    = fullfile(pth,[nam,'_iskull',ext]);
fl_ic   = {[],[],'uint8',[]};
Iisk    = spm_imcalc_ui(Iin,Iisk,'i1+i2+i3',fl_ic);

% Applies some erosion/growing to refine
% the mask and write *_iskull.img
Iarg    = strvcat(Iisk,Iisk);
ne      = flags.ne(1);
ng      = flags.ng(1);
thr_im  = flags.thr_im(1);
Iout    = spm_eeg_inv_ErodeGrow(Iarg,ne,ng,thr_im);

% Generate the outer-scalp volume, if possible
Iscl       = spm_vol(Ivol);
Iscl.dat   = spm_loaduint8(Iscl);
ne         = flags.ne(end);
ng         = flags.ng(end);
thr_im     = flags.thr_im(end);
Iscl.dat    = spm_eeg_inv_ErodeGrow(Iscl.dat,ne,ng,thr_im);

% The bottom part of the image is masked out (in MNI space)
p1      = [0 110 -109]'; % point below on the nose
p2      = [0 -55 -110]'; % point below the bottom of cerebellum (mm)
c       = [p2(1:2)' ((p1(2)-p2(2))^2+p1(3)^2-p2(3)^2)/(2*(-p2(3)+p1(3)))];
R2      = (p2(3)-c(3))^2; % center and radius of sphere used to chop off bottom of the head
		
% X,Y,X coordinates in vox, in original image
X       = (1:Iscl.dim(1))'*ones(1,Iscl.dim(2));
X       = X(:)';
Y       = ones(Iscl.dim(1),1)*(1:Iscl.dim(2));
Y       = Y(:)';
Z       = zeros(Iscl.dim(1),Iscl.dim(2));
Z       = Z(:)'; 
Unit    = ones(size(Z));
isn_tmp = Misn; % use of the simple affine transform

% for pp = 1:Iscl.dim(3)
% end
		
% Write files
Iscl.fname = fullfile(pth,[nam,'_scalp',ext]);
Iscl       = spm_create_vol(Iscl);
spm_progress_bar('Init',Iscl.dim(3),'Writing scalp','planes completed');
for pp=1:Iscl.dim(3)
    pp
    spm_progress_bar('Set',pp);
    XYZ             = Iscl.mat*[X ; Y ; Z+pp ; Unit]; % coord in mm of voxels of scalp img
    XYZ             = spm_get_orig_coord(XYZ(1:3,:)', isn_tmp)';           
    lz              = find(((XYZ(1,:)-c(1)).^2+(XYZ(2,:)-c(2)).^2+(XYZ(3,:)-c(3)).^2-R2)>0);
    if ~lz
        break
    end
    val_pp          = Iscl.dat(:,:,pp);
    val_pp(lz)      = 0;
    Isc.dat(:,:,pp) = val_pp;
    Iscl = spm_write_plane(Iscl,double(Iscl.dat(:,:,pp)),pp);
end
spm_progress_bar('Clear');

% Output arguments
if nargout == 1
    D.inv{val}.mesh.msk_iskull = Iisk;
    D.inv{val}.mesh.msk_scalp  = Iscl.fname;
    D.inv{val}.mesh.msk_flags  = flags;
    if str2num(version('-release'))>=14
        save(fullfile(D.path, D.fname), '-V6', 'D');
    else
        save(fullfile(fpath3, D.fname), 'D');
    end
    varargout{1} = D;
elseif nargout == 2
    varargout{1} = strvcat(Pisk,Iscl.fname);
    varargout{2} = flags;
end

