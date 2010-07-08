function spm_eeg_inv_group(S)
% Source reconstruction for a group ERP or ERF study
% FORMAT spm_eeg_inv_group(S)
%
% S  - string array of names of M/EEG mat files for inversion (optional)
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
% All trials in D.events.types will be inverted in the order specified.
% The result is a contrast (saved in D.mat) and a 3-D volume of MAP or
% conditional estimates of source activity that are constrained to the
% same subset of voxels.  These would normally be passed to a second-level
% SPM for classical inference about between-trial effects, over subjects.
%__________________________________________________________________________
%
% References:
% Electromagnetic source reconstruction for group studies. V. Litvak and
% K.J. Friston. NeuroImage, 42:1490-1498, 2008.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_inv_group.m 3979 2010-07-08 14:53:46Z vladimir $
 
SVNrev = '$Rev: 3979 $';
 
%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
 
%-Check if to proceed
%--------------------------------------------------------------------------
str = questdlg({'This will overwrite previous source reconstructions.', ...
    'Do you wish to continue?'},'M/EEG Group Inversion','Yes','No','Yes');
if ~strcmp(str,'Yes'), return, end
 
% Load data
%==========================================================================
 
% Give file names
%--------------------------------------------------------------------------
if ~nargin
    [S, sts] = spm_select(Inf, 'mat','Select M/EEG mat files');
    if ~sts, return; end
end
Ns    = size(S,1);
swd   = pwd;
 
% Load data and set method
%==========================================================================
for i = 1:Ns
    
    fprintf('checking for previous inversions: subject %i\n',i);
    D{i}                 = spm_eeg_load(deblank(S(i,:)));
    D{i}.val             = 1;
    D{i}.inv{1}.method   = 'Imaging';
    
    % clear redundant models
    %----------------------------------------------------------------------
    D{i}.inv = D{i}.inv(1);
    
    
    % clear previous inversions
    %----------------------------------------------------------------------
    try, D{i}.inv{1} = rmfield(D{i}.inv{1},'inverse' ); end
    try, D{i}.inv{1} = rmfield(D{i}.inv{1},'contrast'); end
    
    % save forward model parameters
    %----------------------------------------------------------------------
    save(D{i});
    
end
 
% Check for existing forward models and consistent Gain matrices
%--------------------------------------------------------------------------
Nd = zeros(1,Ns);
for i = 1:Ns
    fprintf('checking for forward models: subject %i\n',i);
    try
        [L, D{i}] = spm_eeg_lgainmat(D{i});
        Nd(i) = size(L,2);               % number of dipoles
    catch
        Nd(i) = 0;
    end
end
 
% use template head model where necessary
%==========================================================================
if max(Nd > 1024)
    NS = find(Nd ~= max(Nd));            % subjects requiring forward model
else
    NS = 1:Ns;
end
for i = NS
 
    cd(D{i}.path);
 
    % specify cortical mesh size (1 to 4; 1 = 5125, 2 = 8196 dipoles)
    %----------------------------------------------------------------------
    Msize  = 2;
 
    % use a template head model and associated meshes
    %======================================================================
    D{i} = spm_eeg_inv_mesh_ui(D{i}, 1, 1, Msize);
 
    % save forward model parameters
    %----------------------------------------------------------------------
    save(D{i});
 
end
 
% Get inversion parameters
%==========================================================================
inverse = spm_eeg_inv_custom_ui(D{1});
 
% Select modality
%==========================================================================
% Modality
%------------------------------------------------------------------
[mod, list] = modality(D{1}, 1, 1);
if strcmp(mod, 'Multimodal')
    [selection, ok]= listdlg('ListString', list, 'SelectionMode', 'multiple' ,...
        'Name', 'Select modalities' , 'InitialValue', 1:numel(list),  'ListSize', [400 300]);
    if ~ok
        return;
    end
    
    inverse.modality  = list(selection);
    
    if numel(inverse.modality) == 1
        inverse.modality = inverse.modality{1};
    end
else
    inverse.modality = mod;
end
 
for i = 2:Ns
    [mod, list] = modality(D{i}, 1, 1);
    if ~all(ismember(inverse.modality, list))
        error([inverse.modality ' modality is missing from ' D{i}.fname]);
    end
end
 
% and save them (assume trials = types)
%--------------------------------------------------------------------------
for i = 1:Ns
    D{i}.inv{1}.inverse = inverse;
end
 
% specify time-frequency window contrast
%==========================================================================
tfwin = spm_input('Time-Frequency contrast?','+1','y/n',[1,0],1);
if tfwin
 
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
    contrast.smooth  = 4;
    
    str  = {'evoked','induced'};
    contrast.type = spm_input('Power of the energy or mean energy','+1','b',str,[],1);    
else
    contrast = [];
end
 
% Register and compute a forward model
%==========================================================================
for i = NS
 
    fprintf('Registering and computing forward model (subject: %i)\n',i);
       
    % Forward model
    %----------------------------------------------------------------------
    D{i} = spm_eeg_inv_datareg_ui(D{i}, 1);
    try
        D{i}.inv{1}.forward.voltype = voltype;
        D{i}    = spm_eeg_inv_forward(D{i});
    catch
        D{i}    = spm_eeg_inv_forward_ui(D{i});
        voltype = D{i}.inv{1}.forward.voltype;
    end
    
    % save forward model
    %----------------------------------------------------------------------
    save(D{i});
 
end
 
% Invert the forward model
%==========================================================================
D     = spm_eeg_invert(D);
if ~iscell(D), D = {D}; end
 
% Save
%==========================================================================
for i = 1:Ns
    save(D{i});
end
clear D
 
 
% Compute conditional expectation of contrast and produce image
%==========================================================================
if ~isempty(contrast)
 
    % evaluate contrast and write image
    %----------------------------------------------------------------------
    for i = 1:Ns
        D     = spm_eeg_load(deblank(S(i,:)));
        D.inv{1}.contrast = contrast;
        D     = spm_eeg_inv_results(D);
        D     = spm_eeg_inv_Mesh2Voxels(D);
        save(D);
    end
end
 
% Cleanup
%==========================================================================
cd(swd);

