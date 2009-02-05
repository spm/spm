function D = spm_eeg_invert_fuse_ui(S)
% EEG-MEG fusion
% FORMAT D = spm_eeg_invert_fuse_ui(S)
%
% S  - matrix of names of EEG/MEG mat files for inversion
%__________________________________________________________________________
%
% This is the user interface for spm_eeg_invert_fuse and will save a new
% D.mat structure with the fused MAP estimates of source activity.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_invert_fuse_ui.m 2696 2009-02-05 20:29:48Z guillaume $

% Load data
%==========================================================================

% Gile file names
%--------------------------------------------------------------------------
if ~nargin
    S = spm_select(Inf, 'mat','Select EEG/MEG mat files');
end
Ns    = size(S,1);
PWD   = pwd;
val   = 1;

% Load data and set method
%==========================================================================
for i = 1:Ns
    [p f]                = fileparts(S(i,:));
    D{i}                 = spm_eeg_ldata(fullfile(p,f));
    D{i}.path            = p;
    D{i}.fname           = f;
end

% Check for existing forward models and consistent Gain matrices
%--------------------------------------------------------------------------
for i = 1:Ns
    cd(D{i}.path)
    try
        [L, D{i}]  = spm_eeg_lgainmat(D{i});
        Nd(i)  = size(L,2);                             % number of dipoles
    catch
        Nd(i)  = 0;
    end
end

% use template head model where necessary
%======================================================================
NS    = find(Nd ~= 7204);
for i = NS

    cd(D{i}.path)

    % specify cortical mesh size (1 to 4; 1 = 3004, 4 = 7204 dipoles)
    %----------------------------------------------------------------------
    D{i}.inv{val}.mesh.Msize  = 4;

    % use a template head model and associated meshes
    %======================================================================
    D{i} = spm_eeg_inv_template(D{i});

    % get fiducials and sensor locations
    %----------------------------------------------------------------------
    try
        sensors = D{i}.inv{val}.datareg.sensors;
    catch
        if strcmp(questdlg('Use Polhemus file?'),'Yes');

            % get fiduicials and headshape
            %--------------------------------------------------------------
            pol_skip            = 2;
            pol_file            = spm_select(1,'.*\.pol$','Select Polhemus file');
            [fid_eeg,headshape] = spm_eeg_inv_ReadPolhemus(pol_file,pol_skip);
            D{i}.inv{val}.datareg.fid_eeg = fid_eeg;
            
            % get sensor locations
            %--------------------------------------------------------------
            if strcmp(D{i}.modality,'EEG')
                sensors = headshape;
            else
                sensors = load(spm_select(1,'mat','Select MEG sensor file'));
                name    = fieldnames(sensors);
                sensors = getfield(sensors,name{1});
            end
            
        else
            [f p]   = uigetfile('*sens*.mat','select sensor locations');
            sensors = load(fullfile(p,f));
            name    = fieldnames(sensors);
            sensors = getfield(sensors,name{1});
        end
    end
    
    try
        fid_eeg = D{i}.inv{val}.datareg.fid_eeg;
    catch
        [f p]   = uigetfile('*fid*.mat','select fiducial locations');
        fid_eeg = load(fullfile(p,f));
        name    = fieldnames(fid_eeg);
        fid_eeg = getfield(fid_eeg,name{1});
    end

    % get sensor locations
    %----------------------------------------------------------------------
    if strcmp(D{i}.modality,'EEG')
        headshape  = sensors;
        orient     = sparse(0,3);
    else
        headshape  = sparse(0,3);
        try
            orient = D{i}.inv{val}.datareg.megorient;
        catch
            [f p]  = uigetfile('*or*.mat','select sensor orientations');
            orient = load(fullfile(p,f));
            name   = fieldnames(orient);
            orient = getfield(orient,name{1});
        end
    end
    D{i}.inv{val}.datareg.sensors   = sensors;
    D{i}.inv{val}.datareg.fid_eeg   = fid_eeg;
    D{i}.inv{val}.datareg.headshape = headshape;
    D{i}.inv{val}.datareg.megorient = orient;

    % specify forward model
    %----------------------------------------------------------------------
    if strcmp(D{i}.modality,'EEG')
        D{i}.inv{val}.forward.method = 'eeg_3sphereBerg';
    else
        D{i}.inv{val}.forward.method = 'meg_sphere';
    end
end

% Get conditions or trials
%==========================================================================
Nt     = length(D{1}.events.types);
trials = 1:Nt;
if Nt > 1
    if spm_input('All trials?','+1','b',{'yes|no'},[0 1],0)
        str    = sprintf('trials (1:%d)',Nt);
        trials = spm_input(str,'+1','r',trials);
    end
end
for i = 1:Ns
    % and specify inversion parameters
    %----------------------------------------------------------------------
    D{i}.inv{val}.inverse.trials = D{i}.events.types(trials);
    D{i}.inv{val}.inverse.type   = 'GS';
end


% Register and compute a forward model
%==========================================================================
for i = NS

    % Forward model
    %----------------------------------------------------------------------
    fprintf('Registering and computing forward model (subject: %i)\n',i)
    D{i} = spm_eeg_inv_datareg(D{i});
    D{i} = spm_eeg_inv_BSTfwdsol(D{i});
    save(D{i});
    
    % save D
    %----------------------------------------------------------------------
    save(D{i});
    
end

% Invert the forward model
%==========================================================================
D      = spm_eeg_invert_fuse(D);
D.path = PWD;

% Save
%==========================================================================
save(D);
cd(PWD)
