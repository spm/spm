function [mesh] = spm_eeg_inv_segment(mesh)
% Spatial Normalization (using SPM8's segment toolbox) segments the
% individual sMRI (if necessary) 
%
% FORMAT mesh = spm_eeg_inv_spatnorm(mesh)
% Input:
% mesh     - input mesh data struct (optional)
% Output:
% D        - same data struct including the inverse deformation .mat file
%            and filename of noramlised (bias correctd) sMRI
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout, Vladimir Litvak
% $Id: spm_eeg_inv_segment.m 2419 2008-10-30 19:40:32Z vladimir $

% initialise
%--------------------------------------------------------------------------
try
    sMRI = mesh.sMRI;
catch
    sMRI = spm_select(1,'image','Select subject''s structural MRI');
    mesh.sMRI = sMRI;
end

%--------------------------------------------------------------------------
[pth, nam] = spm_fileparts(sMRI);

try
    for i = 1:5
        v = spm_vol(fullfile(pth, ['c' num2str(i) nam '.nii']));
    end
    load(fullfile(pth, [nam '_seg8.mat']));
catch
    fprintf(['\n\tSegmenting sMRI and computing mapping from canonical\n', ...
        '\tspace to subject''s MRI space.\n\n']);
    
    % This is necessary to run the segmentation as the new toolbox can only be
    % accessed via the jobman at the moment.
    spm('defaults','eeg');
    spm_jobman('initcfg');
    clear jobs
    jobs{1}.tools{1}.preproc8.channel.vols{1} = sMRI;
    spm_jobman('run',jobs);
end

tmp=load(fullfile(pth, [nam '_seg8.mat']), 'Affine');
mesh.Affine = tmp.Affine;

% % Spatial Transformation into MNI space
% %--------------------------------------------------------------------------
% def_name      = [nam '_sn_'     num2str(val) '.mat'];
% isndef_name   = [nam '_inv_sn_' num2str(val) '.mat'];
% if exist(fullfile(pth, def_name),'file') && exist(fullfile(pth, isndef_name),'file')
%     % reuse  the files if they're available
%     % comment this if you want to recalculte it anyway.
%     mesh.def      = fullfile(pth,def_name);
%     mesh.invdef   = fullfile(pth,isndef_name);    
% end
% if ~isfield(mesh,'def') || isempty(mesh)
%     res           = spm_preproc(sMRI);
%     [sn,isn]      = spm_prep2sn(res); %#ok<NASGU>
%     def_name      = [nam '_sn_'     num2str(val) '.mat'];
%     isndef_name   = [nam '_inv_sn_' num2str(val) '.mat'];
%     mesh.def      = fullfile(pth,def_name);
%     mesh.invdef   = fullfile(pth,isndef_name);
%     save(mesh.def, '-STRUCT', 'sn');
%     save(mesh.invdef, '-STRUCT', 'isn');
% else
%     fprintf('\tNormalisation parameters already exist.\n')
% end
% 
% % Writing the segments (subject space)
% %--------------------------------------------------------------------------
% if ~exist(fullfile(pth,['c1',nam,ext]),'file')
%     opts = struct('biascor',1,...
%                   'GM',     [0 0 1],...
%                   'WM',     [0 0 1],...
%                   'CSF',    [0 0 1],...
%                   'cleanup',0);
%     spm_preproc_write(sn,opts);
%     mesh.nobias = fullfile(pth,['m' nam ext]);
% else
%     fprintf('\tSegmentation images already exist.\n')
%     if ~isfield(mesh,'nobias')
%         mesh.nobias = fullfile(pth,['m' nam ext]);
%     end
% end
%     
% 
% % Writing the wmsMRI (MNI space: bias corrected in 1mm voxels)
% %--------------------------------------------------------------------------
% if ~exist(fullfile(pth,['wm',nam,ext]),'file')
%     flags.vox    = [1 1 1];
%     spm_write_sn(mesh.nobias,mesh.def,flags);
%     mesh.wmMRI = fullfile(pth,['wm' nam ext]);
% elseif ~isfield(mesh,'wmMRI') || isempty(mesh.wmMRI)
%     mesh.wmMRI = fullfile(pth,['wm' nam ext]);
% else
%     fprintf('\tNormalised structural image already exist.\n')
% end    

% finished
%--------------------------------------------------------------------------
spm('Pointer','Arrow');



