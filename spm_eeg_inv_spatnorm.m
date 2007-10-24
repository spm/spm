function varargout = spm_eeg_inv_spatnorm(varargin)

%==========================================================================
% Spatial Normalization (using a unified model - SPM5)
% transforms individual sMRI into MNI T1 space and
% saves the [inverse] deformations (...vbm_inv_sn.mat)
% that will be needed for computing the individual mesh
%
% FORMAT D = spm_eeg_inv_spatnorm(D,val)
% Input:
% D        - input data struct (optional)
% Output:
% D        - same data struct including the inverse deformation .mat file
%            and filename of noramlised (bias correctd) sMRI
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_spatnorm.m 972 2007-10-24 11:48:35Z stefan $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

try
    sMRI = D.inv{val}.mesh.sMRI;
catch
    sMRI = spm_select([1],'image','Select subject''s structural MRI');
end

disp({'Normalising sMRI and computing mapping from canonical';
      'space to subject''s MRI space'})
%--------------------------------------------------------------------------
[pth,nam,ext] = spm_fileparts(sMRI);

% Spatial Transformation into MNI space
%--------------------------------------------------------------------------
res           = spm_preproc(sMRI);
[sn,isn]      = spm_prep2sn(res);
def_name      = [nam '_sn_'     num2str(val) '.mat'];
isndef_name   = [nam '_inv_sn_' num2str(val) '.mat'];
D.inv{val}.mesh.def    = fullfile(pth,def_name);
D.inv{val}.mesh.invdef = fullfile(pth,isndef_name);
save(D.inv{val}.mesh.def, '-STRUCT', 'sn');
save(D.inv{val}.mesh.invdef, '-STRUCT', 'isn');

% Writing the segments (subject space)
%--------------------------------------------------------------------------
opts = struct('biascor',1,...
              'GM',     [0 0 1],...
              'WM',     [0 0 1],...
              'CSF',    [0 0 1],...
              'cleanup',0);
spm_preproc_write(sn,opts);
D.inv{val}.mesh.nobias = fullfile(pth,['m' nam ext]);

% Writing the wmsMRI (MNI space: bias corrected in 1mm voxels)
%--------------------------------------------------------------------------
flags.vox    = [1 1 1];
spm_write_sn(D.inv{val}.mesh.nobias,D.inv{val}.mesh.def,flags);
D.inv{val}.mesh.wmMRI = fullfile(pth,['wm' nam ext]);

% finished
%--------------------------------------------------------------------------
varargout{1} = D;
spm('Pointer','Arrow');
disp(D.inv{val}.mesh)


