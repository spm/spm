function spm_eeg_inv_group(S);
% Source reconstruction for a group ERP or ERF study
% FORMAT spm_eeg_inv_group
%
% S  - matrix of names of EEG/MEG mat files for inversion
%__________________________________________________________________________
%
% spm_eeg_inv_group inverts forward models for a group of subjects or ERPs
% under the simple assumption that the [empirical prior] variance on each
% source can be factorised into source-specific and subject-specific terms.
% These covariance components are estimated using ReML (a form of Gaussian
% process modelling to give empirical priors on sources.  Source-specific
% covariance parameters are estimated first using the sample covariance
% matrix in sensor space over subjects and trials using multiple sparse
% priors (and,  by default, a greedy search).  The subject-specific terms
% are then estimated by pooling over trials for each subject separately.
% The result is a contrast (saved in D.mat) and a 3-D volume of MAP or
% conditional estimates of source activity that are constrained to the
% same subset of voxels.  These would normally be passed to a second-level
% SPM for classical inference about between-trial effects, over subjects.
%--------------------------------------------------------------------------
str = questdlg('this will overwrite previous source reconstructions OK?');
if ~strcmp(str,'Yes'), return, end
 
% Load data
%==========================================================================
 
% Gile file names
%--------------------------------------------------------------------------
if ~nargin
    S = spm_select(Inf, '.mat','Select EEG/MEG mat files');
end
Ns    = size(S,1);
PWD   = pwd;
val   = 1;
 
for i = 1:Ns
 
    % Load data and set method
    %======================================================================
    D{i}                 = spm_eeg_ldata(S(i,:));
    D{i}.val             = val;
    D{i}.inv{val}.method = 'Imaging';
    cd(D{i}.path);
 
    % specify cortical mesh size (1 tp 4; 1 = 3004, 4 = 7204 dipoles)
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
        [f p]   = uigetfile('*sens*.mat','select sensor locations');
        sensors = load(fullfile(p,f));
        name    = fieldnames(sensors);
        sensors = getfield(sensors,name{1});
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
        headshape = sensors;
        orient    = sparse(0,3);
    else
        headshape = sparse(0,3);
        try
            D{i}.inv{val}.datareg.megorient;
        catch
            [f p]   = uigetfile('*or*.mat','select sensor orientations');
            orient  = load(fullfile(p,f));
            name    = fieldnames(orient);
            orient  = getfield(orient,name{1});
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
if length(D{1}.events.types) > 1
    if spm_input('All conditions or trials','+1','b',{'yes|no'},[1 0],1)
        trials = D{1}.events.types;
    else
        trials = [];
        for  i = 1:length(D{1}.events.types)
            str = sprintf('invert %i',D{1}.events.types(i))
            if spm_input(str,'+1','b',{'yes|no'},[1 0],1);
                trials(end + 1) = D{1}.events.types(i);
            end
        end
    end
else
    trials = D{1}.events.types;
end
 
% specify inversion parameters
%--------------------------------------------------------------------------
inverse.trials = trials;                      % Trials or conditions
inverse.type   = 'GS';                        % Priors; GS, MSP, LOR or IID
 
% specify time-frequency window contrast
%==========================================================================
 
% get time window
%--------------------------------------------------------------------------
woi              = spm_input('Time window (ms)','+1','r',[100 200]);
woi              = sort(woi);
contrast.woi     = round([woi(1) woi(end)]);
 
% get frequency window
%--------------------------------------------------------------------------
fboi             = spm_input('Frequency [band] of interest (Hz)','+1','r',0);
fboi             = sort(fboi);
contrast.fboi    = round([fboi(1) fboi(end)]);
contrast.display = 0;
 
% Register and compute a forward model
%==========================================================================
for i = 1:Ns
 
    fprintf('Registering and computing forward model (subject: %i)\n',i)
    D{i}.inv{val}.inverse = inverse;
 
    % Forward model
    %----------------------------------------------------------------------
    D{i} = spm_eeg_inv_datareg(D{i});
    D{i} = spm_eeg_inv_BSTfwdsol(D{i});
    
end
 
% Invert the forward model
%==========================================================================
D    = spm_eeg_invert(D);

if ~iscell(D)
    D = {D};
end
 
 
% Compute conditional expectation of contrast and produce image
%==========================================================================
for i = 1:Ns
 
    D{i}.inv{val}.contrast = contrast;
 
    % evaluate and write image
    %----------------------------------------------------------------------
    D{i} = spm_eeg_inv_results(D{i});
    D{i} = spm_eeg_inv_Mesh2Voxels(D{i});
 
end
 
% Save
%==========================================================================
S     = D;
for i = 1:Ns
 
    % save D
    %----------------------------------------------------------------------
    D = S{i};
    if spm_matlab_version_chk('7.1') >= 0
        save(fullfile(D.path, D.fname), '-V6', 'D');
    else
        save(fullfile(D.path, D.fname), 'D');
    end
end
