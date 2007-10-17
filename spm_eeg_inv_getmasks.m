function D = spm_eeg_inv_getmasks(varargin);

%=======================================================================
% Generate the binary images (masks) of skull and scalp surfaces
%
% FORMAT [D] = spm_eeg_inv_getmasks(D,val);
% Input:
% D		    - input data struct (optional)
% Output:
% D			- same data struct including the new files and parameters
% D.inv{val}.mesh.msk_iskull;
% D.inv{val}.mesh.msk_cortex;
% D.inv{val}.mesh.msk_scalp;
% D.inv{val}.mesh.msk_flags;
%
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_getmasks.m 954 2007-10-17 15:12:26Z rik $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
fprintf('\n\n'), fprintf('%c','='*ones(1,80)), fprintf('\n')
fprintf(['\tGenerate binary volumes from structural image.\n']);

% flags
%==========================================================================
%  - img_norm: normalised image (1) or not (0)
%  - ne,ng   : number of iterations for erosion and growing
%  - thr_im  : threhold applied to the binary images

% Modified by Rik to allow other masking flags to be passed
try 
  flags = D.inv{val}.mesh.msk_flags;
  if isempty(flags)
     flags = struct('img_norm',0,'ne',1,'ng',2,'thr_im',[.5 .05]);
  end
catch
  flags = struct('img_norm',0,'ne',1,'ng',2,'thr_im',[.5 .05]);
end

%==========================================================================

% flags and filenames
%--------------------------------------------------------------------------
[pth, nam, ext] = spm_fileparts(D.inv{val}.mesh.sMRI);
ne      = flags.ne(1);
ng      = flags.ng(1);
fl_ic   = {[],[],'uint8',[]};

% Add up GM+WM+CSF to produce the skull surface and write *_iskull.img
%--------------------------------------------------------------------------
Iin     = strvcat(fullfile(pth,['c1',nam,ext]), ...
                  fullfile(pth,['c2',nam,ext]), ...
                  fullfile(pth,['c3',nam,ext]));
Iisk    = fullfile(pth,[nam,'_iskull',ext]);
Iisk    = spm_imcalc_ui(Iin,Iisk,'i1+i2+i3',fl_ic);
Iout    = spm_eeg_inv_ErodeGrow(strvcat(Iisk,Iisk),ne,ng,flags.thr_im(1));

% Use GM+WM to produce the cortical surface and write *_cortex.img
%--------------------------------------------------------------------------
Iin     = strvcat(fullfile(pth,['c1',nam,ext]), ...
                  fullfile(pth,['c2',nam,ext]));
Ictx    = fullfile(pth,[nam,'_cortex',ext]);
Ictx    = spm_imcalc_ui(Iin,Ictx,'i1+i2',fl_ic);
Iout    = spm_eeg_inv_ErodeGrow(strvcat(Ictx,Ictx),ne,ng,flags.thr_im(1));

% Use sMRI to produce the scalp surface and write *_scalp.img
%--------------------------------------------------------------------------
Iin     = D.inv{val}.mesh.nobias;
Iscl    = fullfile(pth,[nam,'_scalp',ext]);
Iscl    = spm_imcalc_ui(Iin,Iscl,'i1',fl_ic);
Iout    = spm_eeg_inv_ErodeGrow(strvcat(Iscl,Iscl),ne,ng,flags.thr_im(2));

% Output arguments
%--------------------------------------------------------------------------
D.inv{val}.mesh.msk_iskull = Iisk;
D.inv{val}.mesh.msk_cortex = Ictx;
D.inv{val}.mesh.msk_scalp  = Iscl;
D.inv{val}.mesh.msk_flags  = flags;


