function mesh = spm_eeg_inv_spatnorm(mesh)
% Spatial Normalisation (using Unified Segmentation)
% Transforms individual sMRI into MNI space and save the [inverse] 
% deformations that will be needed for computing the individual mesh
%
% FORMAT mesh = spm_eeg_inv_spatnorm(mesh)
% 
% mesh        - input data struct 
% 
% mesh        - same data struct including the inverse deformation .mat file
%               and filename of noramlised (bias correctd) sMRI
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_spatnorm.m 4389 2011-07-12 11:07:33Z john $

spm('Pointer','Watch');

%-Initialise
%--------------------------------------------------------------------------
sMRI = mesh.sMRI;

%-Spatial Transformation into MNI space
%--------------------------------------------------------------------------
[pth,nam] = spm_fileparts(sMRI);
def       = fullfile(pth,['y_' nam '.nii']);
mat       = fullfile(pth,[nam '_seg8.mat']);

if ~(exist(def, 'file') && exist(mat, 'file'))
    % Assume it is an image, so derive deformation field.
    T = fullfile(spm('dir'),'toolbox','Seg','TPM.nii,');
    p = struct('channel',struct('vols',     {{sMRI}},...
        'biasreg',  0.0001,...
        'biasfwhm', 60,...
        'write',    [0 0]),...
        'tissue', struct('tpm',   {{[T,'1']},{[T,'2']},{[T,'3']},{[T,'4']},{[T,'5']},{[T,'6']}},...
        'ngaus', {2,2,2,3,4,2},...
        'native',{[0 0],[0 0],[0 0],[0 0],[0 0],[0 0]},...
        'warped',{[0 0],[0 0],[0 0],[0 0],[0 0],[0 0]}),...
        'warp',   struct('reg',4, 'affreg', 'mni', 'samp', 3, 'write', [0 1],'mrf',0));

    spm_preproc_run(p);
end

mesh.def    = def;
mesh.Affine = getfield(load(mat, 'Affine'), 'Affine');

%-Clean up
%--------------------------------------------------------------------------
spm('Pointer','Arrow');
