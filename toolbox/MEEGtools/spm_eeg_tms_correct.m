function D = spm_eeg_tms_correct(S)
% Function for removing TMS artefacts
% FORMAT D = spm_eeg_tms_correct(S)
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file

% Output:
% D                   - MEEG object (also written on disk)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% Implements:
%
% Litvak V, Komssi S, Scherg M, Hoechstetter K, Classen J, Zaaroor M, Pratt H, Kahkonen S.
% Artifact correction and source analysis of early electroencephalographic
% responses evoked by transcranial magnetic stimulation over primary motor cortex.
% Neuroimage. 2007; 37(1):56-70.
%
% Vladimir Litvak
% $Id: spm_eeg_tms_correct.m 3602 2009-11-30 13:38:48Z vladimir $

SVNrev = '$Rev: 3602 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Correct TMS artefact');

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

if strcmp(D.type, 'continuous')
    error('The data should be epoched');
end

if ~isfield(D, 'inv')
    error('Prepare and save a forward model before running this function');
end

[L,D] = spm_eeg_lgainmat(D);

save(D);

oD = D;

S = [];
S.D = D;

if D.timeonset < -5e-3
    S.time = 1e3*[D.timeonset -5e-3];
end

D = spm_eeg_bc(S);

S = [];
S.D = D;
S.time = [-2
    15];
D = spm_eeg_interpolate_artefact(S);

delete(S.D);

DD = {};
ncomp = [];

for i = 1:D.ntrials
    tD = reject(D, [], 1);
    tD = reject(tD, i, 0);
    
    tD = path(tD, pwd);
    
    [p f] = fileparts(tD.fname);    
    
    S = [];
    S.D = tD;
    tD = spm_eeg_remove_bad_trials(S);
    
    S = [];
    S.D = tD;
    S.method = 'SVD';
    S.timewin = [-0.002
        0.02];
    S.svdthresh = 30;
    tD = spm_eeg_spatial_confounds(S);
    
    ncomp(end+1) = sum(any(tD.sconfounds));
    
    S = [];
    S.D = tD;
    S.newname = [f '_trial' num2str(i) '.mat'];
    tD = spm_eeg_copy(S); 
    
    delete(S.D);
    
    if ncomp(end)>0
        S = [];
        S.D = tD;
        S.correction = 'Berg';
        DD{end+1} = spm_eeg_correct_sensor_data(S);
        
        delete(S.D);
    else
        DD{end+1} = tD;
    end
end

delete(D);

DD(cellfun('isempty', DD)) = [];

if numel(DD)>1
    S = [];
    S.D = fname(DD{1});
    for f = 2:numel(DD)
        S.D = strvcat(S.D, fname(DD{f}));
    end
    S.recode = 'same';
    D = spm_eeg_merge(S);
    
    fileind =[];
    for f = 1:numel(DD)
        fileind = [fileind f*ones(1, ntrials(DD{f}))];
        delete(DD{f});
    end
    D.fileind = fileind;
else
    D = DD{1};
    D.fileind = ones(1, ntrials(D));
end

D.ncomp = ncomp;

D = reject(D, [], oD.reject);

D = history(D, oD.history, [], 'reset');

D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName', 'Correct TMS artefact: done');
