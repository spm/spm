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
% $Id: spm_eeg_inv_getmasks.m 2419 2008-10-30 19:40:32Z vladimir $

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
         flags = struct('img_norm',0,'ne',{3, 3, 3, 3},'ng',{0, 2, 2, 6},'thr_im', {.05, .05, 0.5, 0.02});
    end
catch
    flags = struct('img_norm',0,'ne',{3, 3, 3, 3},'ng',{0, 2, 2, 6},'thr_im', {.05, .05, 0.5, 0.02});
end

%==========================================================================

% flags and filenames
%--------------------------------------------------------------------------
[pth, nam] = spm_fileparts(mesh.sMRI); ext = '.nii';
fl_ic   = {[],[],'uint8',[]}; % writing option for ImCalc

% Use GM+WM to produce the cortical surface and write *_cortex.img
%--------------------------------------------------------------------------
Ictx    = fullfile(pth,[nam,'_cortex',ext]);
if exist(Ictx,'file')
    fprintf('\tCortex binary volume already exists.\n')
else
    Iin     = strvcat(fullfile(pth,['c1',nam,ext]), ...
        fullfile(pth,['c2',nam,ext]));
    Ictx    = spm_imcalc_ui(Iin,Ictx,'i1+i2',fl_ic);
    Iout    = spm_eeg_inv_ErodeGrow(strvcat(Ictx,Ictx), ...
        flags(1).ne, flags(1).ng, flags(1).thr_im);
end

% Add up GM+WM+CSF to produce the inner skull volume and write *_iskull.img
%--------------------------------------------------------------------------
Iisk    = fullfile(pth,[nam,'_iskull',ext]);
if exist(Iisk,'file')
    fprintf('\tInner skull binary volume already exists.\n')
else
    Iin     = strvcat(fullfile(pth,['c1',nam,ext]), ...
        fullfile(pth,['c2',nam,ext]), ...
        fullfile(pth,['c3',nam,ext]));
    Iisk    = spm_imcalc_ui(Iin,Iisk,'i1+i2+i3',fl_ic);
    Iout    = spm_eeg_inv_ErodeGrow(strvcat(Iisk,Iisk), ...
        flags(2).ne, flags(2).ng ,flags(2).thr_im);
end
%%
% Add up GM+WM+CSF+Skull to produce the outer skull volume and write *_oskull.img
%--------------------------------------------------------------------------
Iosk    = fullfile(pth,[nam,'_oskull',ext]);
if exist(Iosk,'file')
    fprintf('\tOuter skull binary volume already exists.\n')
else
    Iin     = strvcat(fullfile(pth,['c1',nam,ext]), ...
        fullfile(pth,['c2',nam,ext]), ...
        fullfile(pth,['c3',nam,ext]),...
        fullfile(pth,['c4',nam,ext]));
    Iosk    = spm_imcalc_ui(Iin,Iosk,'i1+i2+i3+i4',fl_ic);
    Iout    = spm_eeg_inv_ErodeGrow(strvcat(Iosk,Iosk), ...
        flags(3).ne, flags(3).ng ,flags(3).thr_im);
end


%%
% Use sMRI to produce the scalp surface and write *_scalp.img
%--------------------------------------------------------------------------
Iscp    = fullfile(pth,[nam,'_scalp',ext]);
if exist(Iscp,'file')
    fprintf('\tScalp binary volume already exists.\n')
else
    Iin     = strvcat(fullfile(pth,['c1',nam,ext]), ...
        fullfile(pth,['c2',nam,ext]), ...
        fullfile(pth,['c3',nam,ext]), ...
        fullfile(pth,['c4',nam,ext]), ...
        fullfile(pth,['c5',nam,ext]));

    Iscp    = spm_imcalc_ui(Iin,Iscp,'i1+i2+i3+i4+i5',fl_ic);
    Iout    = spm_eeg_inv_ErodeGrow(strvcat(Iscp,Iscp), ...
        flags(4).ne, flags(4).ng, flags(4).thr_im);
end
%%
% Output arguments
%--------------------------------------------------------------------------
mesh.msk_cortex = Ictx;
mesh.msk_iskull = Iisk;
mesh.msk_oskull = Iosk;
mesh.msk_scalp  = Iscp;
mesh.msk_flags  = flags;

return
