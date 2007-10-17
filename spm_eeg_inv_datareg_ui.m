function D = spm_eeg_inv_datareg_ui(varargin)

%==========================================================================
% Data registration user-interface routine
% commands the EEG/MEG data co-registration within original sMRI space
%
% FORMAT D = spm_eeg_inv_mesh_ui(D,[val])
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the new required files and variables
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_datareg_ui.m 956 2007-10-17 15:19:58Z rik $

% Set-up specfic parameters (POLHEMIUS)
%==========================================================================
pol_skip = 2;                    % first 2 channels are reference channels
pol_fig  = 0;                    % supress graphic when reading pol file
%==========================================================================

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
try
    D.inv{val}.mesh.template;
catch
    D.inv{val}.mesh.template = 0;
end

% For some data formats, fids and headshape encoded in native data formats, 
% i.e, put into D.channels during data conversion	RH 12/9/07
try
    fid_eeg   = D.channels.fid_eeg;
    headshape = D.channels.headshape;
    sensors   = D.channels.Loc';

catch

% Specify sensor space (fiducials sensors and scalp surface)
%==========================================================================
  if strcmp(questdlg('Use Polhemus file?'),'Yes');

    % get fiduicials and headshape
    %----------------------------------------------------------------------
    pol_file            = spm_select(1,'.pol','Select Polhemus file');
    [fid_eeg,headshape] = spm_eeg_inv_ReadPolhemus(pol_file,pol_skip,pol_fig);

    % get sensor locations
    %----------------------------------------------------------------------
    if strcmp(D.modality,'EEG')
        sensors = headshape;
    else
        sensors = load(spm_select(1,'.mat','Select MEG sensor file'));
        name    = fieldnames(sensors);
        sensors = getfield(sensors,name{1});
    end
  else

    % get sensor and fiducial locations
    %----------------------------------------------------------------------
    sensors = load(spm_select(1,'.mat','Select EEG/MEG sensor file'));
    fid_eeg = load(spm_select(1,'.mat','Select EEG/MEG fiducial file'));
    name    = fieldnames(sensors);
    sensors = getfield(sensors,name{1});
    name    = fieldnames(fid_eeg);
    fid_eeg = getfield(fid_eeg,name{1});


    % get headshape and fiducial locations
    %----------------------------------------------------------------------
    if strcmp(questdlg('Use headshape file?'),'Yes');
        headshape = load(spm_select(1,'.mat','Select headshape file (N x 3)'));
        name      = fieldnames(headshape);
        headshape = getfield(headshape,name{1});
        
        % ensure fiducials correspond
        %------------------------------------------------------------------
        % hsp_ind   = inputdlg('Indices of sensor fiducials','SPM',1,{'1 2 3'});
        % hsp_ind   = str2num(hsp_ind{1});
        % headshape = [headshape(hsp_ind,:); headshape];
        
    else
        headshape = sparse(0,3);
    end
  end
end
        
% sensor orientations (MEG)
%--------------------------------------------------------------------------
if strcmp(D.modality,'MEG')
try
    megorient   = D.channels.Orient';
    catch
        megorient = load(spm_select(1,'.mat','Select MEG sensor orients (N x 3)'));
        name      = fieldnames(megorient);
        megorient = getfield(megorient,name{1});
    end
else
    megorient = sparse(0,3);
end


% Speccify anatomical space (fiducials and scalp surface)
%==========================================================================

% fiducial in sMRI space
%--------------------------------------------------------------------------
try
    fid_mri = D.inv{val}.datareg.fid_mri;
catch
    fid_mri = load(spm_select(1,'.mat','Select sMRI fiducials {nasion, left & right ear)'));
    name    = fieldnames(fid_mri);
    fid_mri = getfield(fid_mri,name{1});
end

% get scalp locations if headshape is specified
%--------------------------------------------------------------------------
try
    D.inv{val}.datareg.scalpvert = D.inv{val}.mesh.tess_scalp.vert;
    scalpvert                    = D.inv{val}.datareg.scalpvert;
catch
    if length(headshape)
        scalpvert = spm_select(1,'.mat','Select scalp vertices (N x 3)');
    else
        scalpvert = sparse(0,3);
    end
end


% register
%==========================================================================
D.inv{val}.datareg.sensors   = sensors;
D.inv{val}.datareg.fid_eeg   = fid_eeg;
D.inv{val}.datareg.fid_mri   = fid_mri;
D.inv{val}.datareg.headshape = headshape;
D.inv{val}.datareg.scalpvert = scalpvert;
D.inv{val}.datareg.megorient = megorient;

%--------------------------------------------------------------------------
D = spm_eeg_inv_datareg(D);

% check and display registration
%--------------------------------------------------------------------------
spm_eeg_inv_checkdatareg(D);
