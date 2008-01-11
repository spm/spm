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
% process modelling) to give empirical priors on sources.  Source-specific
% covariance parameters are estimated first using the sample covariance
% matrix in sensor space over subjects and trials using multiple sparse
% priors (and,  by default, a greedy search).  The subject-specific terms
% are then estimated by pooling over trials for each subject separately.
% The result is a contrast (saved in D.mat) and a 3-D volume of MAP or
% conditional estimates of source activity that are constrained to the
% same subset of voxels.  These would normally be passed to a second-level
% SPM for classical inference about between-trial effects, over subjects.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_group.m 1082 2008-01-11 12:50:15Z karl $


% check if to proceed
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

% Load data and set method
%==========================================================================
for i = 1:Ns
    [p f]                = fileparts(S(i,:));
    D{i}                 = spm_eeg_ldata(fullfile(p,f));
    D{i}.path            = p;
    D{i}.fname           = f;
    D{i}.val             = val;
    D{i}.inv{val}.method = 'Imaging';
end

% Check for existing forward models and consistent Gain matrices
%--------------------------------------------------------------------------
for i = 1:Ns
    cd(D{i}.path)
    try
        gainmat   = D{i}.inv{D{i}.val}.forward.gainmat;
        try
            G     = load(gainmat);
        catch
            [p f] = fileparts(gainmat);
            G     = load(f);
            D{i}.inv{D{i}.val}.forward.gainmat = fullfile(pwd,f);
        end
        name   = fieldnames(G);
        L      = sparse(getfield(G, name{1}));
        Nd(i)  = size(L,2);                             % number of dipoles
    catch
        Nd(i)  = 0;
    end
end


% use template head model where necessary
%==========================================================================
NS    = find(Nd ~= 7204);               % subjects to compute forward model
for i = NS

    cd(D{i}.path)

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
        if strcmp(questdlg('Use Polhemus file?'),'Yes');

            % get fiduicials and headshape
            %--------------------------------------------------------------
            pol_skip            = 2;
            pol_file            = spm_select(1,'.pol','Select Polhemus file');
            [fid_eeg,headshape] = spm_eeg_inv_ReadPolhemus(pol_file,pol_skip);

            % get sensor locations
            %--------------------------------------------------------------
            if strcmp(D{i}.modality,'EEG')
                sensors = headshape;
            else
                sensors = load(spm_select(1,'.mat','Select MEG sensor file'));
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
        fid_eeg;
    catch
        try
            fid_eeg = D{i}.inv{val}.datareg.fid_eeg;
        catch
            [f p]   = uigetfile('*fid*.mat','select fiducial locations');
            fid_eeg = load(fullfile(p,f));
            name    = fieldnames(fid_eeg);
            fid_eeg = getfield(fid_eeg,name{1});
        end
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

% Restrictions (not implemented)
%--------------------------------------------------------------------------



% specify time-frequency window contrast
%==========================================================================
str = questdlg('Would you like to specify a time-frequency contrast?');
if strcmp(str,'Yes')

    % get time window
    %----------------------------------------------------------------------
    woi              = spm_input('Time window (ms)','+1','r',[100 200]);
    woi              = sort(woi);
    contrast.woi     = round([woi(1) woi(end)]);

    % get frequency window
    %----------------------------------------------------------------------
    fboi             = spm_input('Frequency [band] of interest (Hz)','+1','r',0);
    fboi             = sort(fboi);
    contrast.fboi    = round([fboi(1) fboi(end)]);
    contrast.display = 0;

else
    contrast = [];
end

% Register and compute a forward model
%==========================================================================
for i = NS

    fprintf('Registering and computing forward model (subject: %i)\n',i)

    % Forward model
    %----------------------------------------------------------------------
    D{i} = spm_eeg_inv_datareg(D{i});
    D{i} = spm_eeg_inv_BSTfwdsol(D{i});

end

% Invert the forward model
%==========================================================================
for i = 1:Ns
    D{i}.inv{val}.inverse = inverse;
end
D     = spm_eeg_invert(D);
if ~iscell(D), D = {D}; end


% Compute conditional expectation of contrast and produce image
%==========================================================================
if length(contrast)

    % evaluate contrast and write image
    %----------------------------------------------------------------------
    for i = 1:Ns
        D{i}.inv{val}.contrast = contrast;
        D{i} = spm_eeg_inv_results(D{i});
        D{i} = spm_eeg_inv_Mesh2Voxels(D{i});

    end
end

% Save
%==========================================================================
S     = D;
for i = 1:Ns
    D = S{i};
    if spm_matlab_version_chk('7.1') >= 0
        save(fullfile(D.path, D.fname), '-V6', 'D');
    else
        save(fullfile(D.path, D.fname), 'D');
    end
end
