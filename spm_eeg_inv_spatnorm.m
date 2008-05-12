function mesh = spm_eeg_inv_spatnorm(mesh,val)
% Spatial Normalization (using a unified model - SPM5) transforms 
% individual sMRI into MNI T1 space and saves the [inverse] deformations 
% (..._inv_sn.mat) that will be needed for computing the individual mesh
%
% FORMAT mesh = spm_eeg_inv_spatnorm(mesh)
% Input:
% mesh     - input mesh data struct (optional)
% Output:
% D        - same data struct including the inverse deformation .mat file
%            and filename of noramlised (bias correctd) sMRI
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_spatnorm.m 1604 2008-05-12 20:34:46Z christophe $

% initialise
%--------------------------------------------------------------------------
try
    sMRI = mesh.sMRI;
catch
    sMRI = spm_select(1,'image','Select subject''s structural MRI');
    mesh.sMRI = sMRI;
end

if nargin<2
    val=1;
end

fprintf(['\n\tNormalising sMRI and computing mapping from canonical\n', ...
        '\tspace to subject''s MRI space.\n\n'])
%--------------------------------------------------------------------------
[pth,nam,ext] = spm_fileparts(sMRI);

% Spatial Transformation into MNI space
%--------------------------------------------------------------------------
def_name      = [nam '_sn_'     num2str(val) '.mat'];
isndef_name   = [nam '_inv_sn_' num2str(val) '.mat'];
if exist(fullfile(pth,def_name),'file') && exist(fullfile(pth,isndef_name),'file')
    % reuse  the files if they're available
    % comment this if you want to recalculte it anyway.
    mesh.def      = fullfile(pth,def_name);
    mesh.invdef   = fullfile(pth,isndef_name);    
end
if ~isfield(mesh,'def') || isempty(mesh)
    res           = spm_preproc(sMRI);
    [sn,isn]      = spm_prep2sn(res); %#ok<NASGU>
    def_name      = [nam '_sn_'     num2str(val) '.mat'];
    isndef_name   = [nam '_inv_sn_' num2str(val) '.mat'];
    mesh.def      = fullfile(pth,def_name);
    mesh.invdef   = fullfile(pth,isndef_name);
    save(mesh.def, '-STRUCT', 'sn');
    save(mesh.invdef, '-STRUCT', 'isn');
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
    mesh.nobias = fullfile(pth,['m' nam ext]);
else
    fprintf('\tSegmentation images already exist.\n')
    if ~isfield(mesh,'nobias')
        mesh.nobias = fullfile(pth,['m' nam ext]);
    end
end
    

% Writing the wmsMRI (MNI space: bias corrected in 1mm voxels)
%--------------------------------------------------------------------------
if ~exist(fullfile(pth,['wm',nam,ext]),'file')
    flags.vox    = [1 1 1];
    spm_write_sn(mesh.nobias,mesh.def,flags);
    mesh.wmMRI = fullfile(pth,['wm' nam ext]);
elseif ~isfield(mesh,'wmMRI') || isempty(mesh.wmMRI)
    mesh.wmMRI = fullfile(pth,['wm' nam ext]);
else
    fprintf('\tNormalised structural image already exist.\n')
end    

% finished
%--------------------------------------------------------------------------
spm('Pointer','Arrow');
disp(mesh)


