function varargout = spm_eeg_inv_spatnorm(varargin)
% Spatial Normalization (using a unified model - SPM5) transforms 
% individual sMRI into MNI T1 space and saves the [inverse] deformations 
% (..._inv_sn.mat) that will be needed for computing the individual mesh
%
% FORMAT D = spm_eeg_inv_spatnorm(D,ival)
% Input:
% D        - input data struct (optional)
% Output:
% D        - same data struct including the inverse deformation .mat file
%            and filename of noramlised (bias correctd) sMRI
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_spatnorm.m 1437 2008-04-17 10:34:39Z christophe $

% initialise
%--------------------------------------------------------------------------
[D,ival] = spm_eeg_inv_check(varargin{:});

try
    sMRI = D.inv{ival}.mesh.sMRI;
catch
    sMRI = spm_select(1,'image','Select subject''s structural MRI');
    D.inv{ival}.mesh.sMRI = sMRI;
end

fprintf(['\n\tNormalising sMRI and computing mapping from canonical\n', ...
        '\tspace to subject''s MRI space.\n\n'])
%--------------------------------------------------------------------------
[pth,nam,ext] = spm_fileparts(sMRI);

% Spatial Transformation into MNI space
%--------------------------------------------------------------------------
if ~isfield(D.inv{ival}.mesh,'def') || isempty(D.inv{ival}.mesh)
    res           = spm_preproc(sMRI);
    [sn,isn]      = spm_prep2sn(res); %#ok<NASGU>
    def_name      = [nam '_sn_'     num2str(ival) '.mat'];
    isndef_name   = [nam '_inv_sn_' num2str(ival) '.mat'];
    D.inv{ival}.mesh.def    = fullfile(pth,def_name);
    D.inv{ival}.mesh.invdef = fullfile(pth,isndef_name);
    save(D.inv{ival}.mesh.def, '-STRUCT', 'sn');
    save(D.inv{ival}.mesh.invdef, '-STRUCT', 'isn');
else
    fprintf('\tNormalisation parameters already exist.\n')
end

% Writing the segments (subject space)
%--------------------------------------------------------------------------
if ~exist(fullfile(pth,['c1',nam,ext]),'file')
    opts = struct('biascor',1,...
                  'GM',     [0 0 1],...
                  'WM',     [0 0 1],...
                  'CSF',    [0 0 1],...
                  'cleanup',0);
    spm_preproc_write(sn,opts);
    D.inv{ival}.mesh.nobias = fullfile(pth,['m' nam ext]);
else
    fprintf('\tSegmentation images already exist.\n')
end
    

% Writing the wmsMRI (MNI space: bias corrected in 1mm voxels)
%--------------------------------------------------------------------------
if ~exist(fullfile(pth,['wm',nam,ext]),'file')
    flags.vox    = [1 1 1];
    spm_write_sn(D.inv{ival}.mesh.nobias,D.inv{ival}.mesh.def,flags);
    D.inv{ival}.mesh.wmMRI = fullfile(pth,['wm' nam ext]);
elseif ~isfield(D.inv{ival}.mesh,'wmMRI') | isempty(D.inv{ival}.mesh.wmMRI)
    D.inv{ival}.mesh.wmMRI = fullfile(pth,['wm' nam ext]);
else
    fprintf('\tNormalised structural image already exist.\n')
end    

% finished
%--------------------------------------------------------------------------
varargout{1} = D;
spm('Pointer','Arrow');
disp(D.inv{ival}.mesh)


