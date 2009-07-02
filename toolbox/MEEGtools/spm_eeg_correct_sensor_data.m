function D = spm_eeg_correct_sensor_data(S)
% Function for extracting source data using dipoles.
% FORMAT D = spm_eeg_correct_sensor_data(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.method        
%
% Output:
% sD                   - MEEG object (also written on disk)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% Implements:
%   Berg P, Scherg M.
%   A multiple source approach to the correction of eye artifacts.
%   Electroencephalogr Clin Neurophysiol. 1994 Mar;90(3):229-41.
%
% Vladimir Litvak
% $Id: spm_eeg_correct_sensor_data.m 3246 2009-07-02 17:22:57Z vladimir $

SVNrev = '$Rev: 3246 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Correct sensor data');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

if ~any(D.sconfounds)    
    D = spm_eeg_spatial_confounds(S);
    if ~any(D.sconfounds)    
        return;
    end
end

A = D.sconfounds;

if ~isfield(S, 'correction')
    S.correction = spm_input('Correction method','+1', 'SSP|Berg', strvcat('SSP', 'Berg'));
end

label = D.chanlabels(setdiff(meegchannels(D, 'MEEG'), badchannels(D)));

montage = [];
montage.labelorg = D.chanlabels;
montage.labelnew = montage.labelorg;
montage.tra      = eye(D.nchannels);

if size(A, 1)~=numel(label)
    error('Spatial confound vector does not match the channels');
end

if isequal(lower(S.correction), 'berg')
    [ok, D] = check(D, 'sensfid');
    
    if ~ok
        if check(D, 'basic')
            errordlg(['The requested file is not ready for source reconstruction.'...
                'Use prep to specify sensors and fiducials.']);
        else
            errordlg('The meeg file is corrupt or incomplete');
        end
        return
    end

    %% ============ Find or prepare head model

    if ~isfield(D, 'val')
        D.val = 1;
    end

    if ~isfield(D, 'inv') || ~iscell(D.inv) ||...
            ~(isfield(D.inv{D.val}, 'forward') && isfield(D.inv{D.val}, 'datareg')) ||...
            ~isa(D.inv{D.val}.mesh.tess_ctx, 'char') % detects old version of the struct
        D = spm_eeg_inv_mesh_ui(D, D.val);
        D = spm_eeg_inv_datareg_ui(D, D.val);
        D = spm_eeg_inv_forward_ui(D, D.val);
    end

    [L, D] = spm_eeg_lgainmat(D);
    
    B = spm_svd(L*L', 0.01);
    
    lim = min(0.5*size(L, 1), 45); % 45 is the number of dipoles BESA would use.
    
    if size(B, 2) > lim;
        B = B(:, 1:lim);
    end
    
    SX = full([A B]);
    
    SXi = pinv(SX);
    SXi = SXi(1:size(A, 2), :);
    
    tra = eye(size(A, 1)) - A*SXi; 
else
    tra = eye(size(A, 1)) - A*pinv(A); 
end
   
[sel1, sel2]  = spm_match_str(montage.labelorg, label);
[sel3, sel4]  = spm_match_str(montage.labelnew, label);

montage.tra(sel3, sel1) = tra(sel4, sel2);

%% ============  Use the montage functionality to compute source activity.
S1   = [];
S1.D = D;
S1.montage = montage;
S1.keepothers = 'yes';

Dnew = spm_eeg_montage(S1);

if ~isempty(D.badchannels)
    Dnew = badchannels(Dnew, D.badchannels, 1);
end

D = Dnew;

D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName', 'Correct sensor data: done');

