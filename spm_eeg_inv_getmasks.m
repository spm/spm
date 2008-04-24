function mesh = spm_eeg_inv_getmasks(mesh)
% Generate the binary images (masks) of skull and scalp surfaces
%
% FORMAT [mesh] = spm_eeg_inv_getmasks(mesh);
% Input:
% mesh         - input mesh struct (optional)
% Output:
% mesh         - same data struct including the new files and parameters
%   mesh.msk_cortex     - cortex mask built from GM+WM
%   mesh.msk_iskull     - inner skull volume, GM+WM+CSF
%   mesh.msk_oskull     - not built here, 'cos no skull segmentation
%   mesh.msk_scalp      - thersholded MRI
%   mesh.msk_flags      - a few standard values
%
% NOTE: so far, SPM segmentation does not allow us to properly segment the
% skull volume, so the "outer skull" volume is left empty. But this volume
% does exist in the template volume (coming from FT)...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_getmasks.m 1477 2008-04-24 14:33:47Z christophe $

% initialise
%--------------------------------------------------------------------------
fprintf('\n\n'), fprintf('%c','='*ones(1,80)), fprintf('\n')
fprintf(['\tGenerate binary volumes from structural image.\n']);

% flags
%==========================================================================
%  - img_norm: normalised image (1) or not (0)
%  - ne,ng   : number of iterations for erosion and growing
%  - thr_im  : threhold applied to the binary images

% Modified by Rik to allow other masking flags to be passed
try 
  flags = mesh.msk_flags;
  if isempty(flags)
     flags = struct('img_norm',0,'ne',1,'ng',2,'thr_im',[.5 .05]);
  end
catch
  flags = struct('img_norm',0,'ne',1,'ng',2,'thr_im',[.5 .05]);
end

%==========================================================================

% flags and filenames
%--------------------------------------------------------------------------
[pth, nam, ext] = spm_fileparts(mesh.sMRI);
fl_ic   = {[],[],'uint8',[]}; % writing option for ImCalc

% Use GM+WM to produce the cortical surface and write *_cortex.img
%--------------------------------------------------------------------------
Iin     = strvcat(fullfile(pth,['c1',nam,ext]), ...
                  fullfile(pth,['c2',nam,ext]));
Ictx    = fullfile(pth,[nam,'_cortex',ext]);
Ictx    = spm_imcalc_ui(Iin,Ictx,'i1+i2',fl_ic);
Iout    = spm_eeg_inv_ErodeGrow(strvcat(Ictx,Ictx), ...
                                flags.ne(1),0,flags.thr_im(1));

% Add up GM+WM+CSF to produce the inner skull volume and write *_iskull.img
%--------------------------------------------------------------------------
Iin     = strvcat(fullfile(pth,['c1',nam,ext]), ...
                  fullfile(pth,['c2',nam,ext]), ...
                  fullfile(pth,['c3',nam,ext]));
Iisk    = fullfile(pth,[nam,'_iskull',ext]);
Iisk    = spm_imcalc_ui(Iin,Iisk,'i1+i2+i3',fl_ic);
Iout    = spm_eeg_inv_ErodeGrow(strvcat(Iisk,Iisk), ...
                                flags.ne(1),flags.ng(1),flags.thr_im(1));

% Maybe one day we'll produce the outer skull volume
%--------------------------------------------------------------------------

% Use sMRI to produce the scalp surface and write *_scalp.img
%--------------------------------------------------------------------------
Iin     = mesh.nobias;
Iscl    = fullfile(pth,[nam,'_scalp',ext]);
Iscl    = spm_imcalc_ui(Iin,Iscl,'i1',fl_ic);
Iout    = spm_eeg_inv_ErodeGrow(strvcat(Iscl,Iscl), ...
                                flags.ne(1),flags.ng(1),flags.thr_im(2));

% Output arguments
%--------------------------------------------------------------------------
mesh.msk_cortex = Ictx;
mesh.msk_iskull = Iisk;
mesh.msk_oskull = [];
mesh.msk_scalp  = Iscl;
mesh.msk_flags  = flags;

return
