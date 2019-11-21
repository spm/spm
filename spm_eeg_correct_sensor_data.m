function D = spm_eeg_correct_sensor_data(S)
% Remove artefacts from the data based on their topography
% FORMAT D = spm_eeg_correct_sensor_data(S)
%
% S      - input structure (optional)
% (optional) fields of S:
%   S.D    - MEEG object or filename of M/EEG mat-file
%   S.mode - 'SSP': simple projection
%          - 'Berg': the method of Berg (see the reference below)
% Output:
% D      - MEEG object (also written on disk)
%
% Implements:
%   Berg P, Scherg M.
%   A multiple source approach to the correction of eye artifacts.
%   Electroencephalogr Clin Neurophysiol. 1994 Mar;90(3):229-41.
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_correct_sensor_data.m 7701 2019-11-21 21:50:17Z vladimir $

SVNrev = '$Rev: 7701 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Correct sensor data');

if ~isfield(S, 'mode') && isfield(S, 'correction'),  S.mode  = S.correction;  end
if ~isfield(S, 'prefix'),                            S.prefix   = 'T';        end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

inputfile = fullfile(D);

if ~any(D.sconfounds)
    D = spm_eeg_spatial_confounds(S);
    if ~any(D.sconfounds)
        return;
    end
end

[mod, list] = modality(D, 1, 1);

A = {};
if isequal(mod, 'Multimodal')
    sconf = getfield(D, 'sconfounds');
    
    for i = 1:numel(list)
        chanind = indchantype(D, list{i}, 'GOOD');
        [sel1, sel2] = spm_match_str(chanlabels(D, chanind), sconf.label);
        
        if any(sconf.bad(sel2))
            error(['Channels ' sprintf('%s ', sconf.label{sel2}) ' should be set to bad.']);
        end
        
        A{i} = sconf.coeff(sel2, :);
    end
else
    A = {D.sconfounds};
    list = {mod};
end

Dorig = D;

for i = 1:numel(A)
    label = D.chanlabels(indchantype(D, list{i}, 'GOOD'))';
    
    montage = [];
    montage.labelorg = label;
    montage.labelnew = label;
    
    montage.chantypeold = lower(D.chantype(D.indchannel(label)))';
    montage.chantypenew  = lower(montage.chantypeold);
   
    montage.chanunitold = D.units(D.indchannel(label))';
    montage.chanunitnew  = montage.chanunitold;
    
    if size(A{i}, 1)~=numel(label)
        error('Spatial confound vector does not match the channels.');
    end
    
    if isequal(lower(S.mode), 'berg')
        % These are the locations taken from the file BR_Brain Regions_LR.bsa
        % in BESA distribution, transformed to Tailarach coordinates using
        % BESA simulator and then to MNI template space using tal2icbm_spm 
        % from http://brainmap.org/icbm2tal/
        sources = [
            -47.3124    6.4922  -10.3381
            -49.1870  -38.7590   -4.3574
            -35.2804   38.7888   21.0944
            -41.0574  -16.4018   46.0553
            -34.9993  -70.8348   20.8450
            0.8786   60.5733   -5.2797
            1.5716   38.5344   44.5853
            1.9792  -16.5380   65.3503
            1.8525  -71.0130   44.3363
            1.2689  -91.6233   -5.6259
            49.1987    6.8326  -12.0156
            51.4597  -38.4040   -6.1068
            37.8063   39.0466   19.8240
            44.5032  -16.1000   44.5681
            38.0874  -70.5770   19.5747
            ];
        
        [D, ok] = check(D, 'sensfid');
        
        if ~ok
            if check(D, 'basic')
                error(['The requested file is not ready for source reconstruction.'...
                    'Use prep to specify sensors and fiducials.']);
            else
                error('The meeg file is corrupt or incomplete');
            end
        end
        
        %-Find or prepare head model
        %==================================================================
        
        if ~isfield(D, 'val')
            D.val = 1;
        end
        
        if ~isfield(D, 'inv') || ~iscell(D.inv) ||...
                ~(isfield(D.inv{D.val}, 'forward') && isfield(D.inv{D.val}, 'datareg')) ||...
                ~isa(D.inv{D.val}.mesh.tess_ctx, 'char') % detects old version of the struct
            D = spm_eeg_inv_mesh_ui(D, D.val);
            D = spm_eeg_inv_datareg_ui(D, D.val);
            D = spm_eeg_inv_forward_ui(D, D.val);
            
            save(D);
        end
        
        fwd = spm_eeg_inv_get_vol_sens(D, D.val, 'MNI-aligned', 'inv', list{i});
        
        [vol, sens] = ft_prepare_vol_sens(fwd.(list{i}).vol, fwd.(list{i}).sens, 'channel', label);
        
        
        L = ft_compute_leadfield(spm_eeg_inv_transform_points(inv(fwd.transforms.toMNI), sources), sens, vol);
        %[L, D] = spm_eeg_lgainmat(D, [], label);
        
        B = spm_svd(L*L', 0.01);
        
        lim = min(0.5*size(L, 1), 45); % 45 is the number of dipoles BESA would use.
        
        if size(B, 2) > lim;
            B = B(:, 1:lim);
        end
        
        SX = full([A{i} B]);
        
        SXi = pinv(SX);
        SXi = SXi(1:size(A{i}, 2), :);
        
        montage.tra = eye(size(A{i}, 1)) - A{i}*SXi;
    else
        montage.tra = eye(size(A{i}, 1)) - A{i}*pinv(A{i});
    end    
     
    S1   = [];
    S1.D = D;
    S1.montage = montage;
    S1.keepothers = true;
    S1.updatehistory  = false;
    
    Dnew = spm_eeg_montage(S1); 
    
    if isfield(D,'inv')
        Dnew.inv = D.inv;
    end
    
    if i>1
        delete(D);
    end
    
    D = Dnew;
end

%-Change the channel order to the original order
%==========================================================================
tra = eye(D.nchannels);
montage = [];
montage.labelorg = D.chanlabels';
montage.labelnew = Dorig.chanlabels';

montage.chantypeold  = lower(D.chantype)';
montage.chantypenew  = lower(Dorig.chantype)';

montage.chanunitold  = D.units';
montage.chanunitnew  = Dorig.units';

[sel1, sel2]  = spm_match_str(montage.labelnew, montage.labelorg);

montage.tra = tra(sel2, :);

S1   = [];
S1.D = D;
S1.montage = montage;
S1.keepothers = true;
S1.updatehistory  = 0;

Dnew = spm_eeg_montage(S1);

delete(D);
D = Dnew;

if ~isempty(badchannels(Dorig))
    D = badchannels(D, badchannels(Dorig), 1);
end

D = D.history(mfilename, S);
save(D);

D = move(D, spm_file(inputfile, 'prefix', S.prefix));

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName', 'Correct sensor data: done');
