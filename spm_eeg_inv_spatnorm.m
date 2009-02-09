function mesh = spm_eeg_inv_spatnorm(mesh)

%==========================================================================
% Spatial Normalization (using a unified model - SPM5)
% transforms individual sMRI into MNI T1 space and
% saves the [inverse] deformations (...vbm_inv_sn.mat)
% that will be needed for computing the individual mesh
%
% FORMAT mesh = spm_eeg_inv_spatnorm(mesh)
% Input:
% mesh        - input data struct 
% Output:
% mesh        - same data struct including the inverse deformation .mat file
%               and filename of noramlised (bias correctd) sMRI
%==========================================================================
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_spatnorm.m 2720 2009-02-09 19:50:46Z vladimir $

% initialise
%--------------------------------------------------------------------------
sMRI = mesh.sMRI;


disp({'Normalising sMRI and computing mapping from canonical';
      'space to subject''s MRI space'})
spm('Pointer','Watch');
%--------------------------------------------------------------------------
[pth,nam,ext] = spm_fileparts(sMRI);

% Spatial Transformation into MNI space
%--------------------------------------------------------------------------
def  = fullfile(pth,['y_' nam '.nii']);
mat  = fullfile(pth,[nam '_seg8.mat']);

if ~(exist(def, 'file') && exist(mat, 'file'))
    % Assume it is an image, so derive deformation field.
    T = fullfile(spm('dir'),'toolbox','Seg','TPM.nii,');
    p = struct('channel',struct('vols',     {{sMRI}},...
        'biasreg',  0.0001,...
        'biasfwhm', 60,...
        'write',    [0 0]),...
        'tissue', struct('tpm',   {{[T,'1']},{[T,'2']},{[T,'3']},{[T,'4']},{[T,'5']},{[T,'6']}},...
        'ngaus', {2,2,2,3,4,2},...
        'native',{[0 0],[0 0],[0 0],[0 0],[0 0],[0 0]},....
        'warped',{[0 0],[0 0],[0 0],[0 0],[0 0],[0 0]}),...
        'warp',   struct('reg',4, 'affreg', 'mni', 'samp', 3, 'write', [1 0]));

    spm_preproc_run(p);
end

mesh.def = def;
mesh.Affine = getfield(load(mat, 'Affine'), 'Affine');

    % finished
%--------------------------------------------------------------------------
spm('Pointer','Arrow');


