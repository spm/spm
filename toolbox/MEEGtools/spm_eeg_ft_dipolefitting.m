% Demo script for interactive artefact rejection using Fieldtrip
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%

% Vladimir Litvak
% $Id: spm_eeg_ft_dipolefitting.m 2333 2008-10-13 13:19:22Z vladimir $

[Finter,Fgraph] = spm('FnUIsetup','Fieldtrip dipole fitting', 0);
%%

%% ============ Load SPM EEG file and verify consistency

D = spm_eeg_load;

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

modality = spm_eeg_modality_ui(D);


%% ============ Find or prepare head model

if ~isfield(D, 'val')
    D.val = 1;
end


try
    vol = D.inv{D.val}.forward.vol;
    datareg = D.inv{D.val}.datareg;
catch
    D = spm_eeg_inv_mesh_ui(D, D.val, [], 1);
    D = spm_eeg_inv_datareg_ui(D, D.val);
    datareg = D.inv{D.val}.datareg;
end

vol = D.inv{D.val}.forward.vol;
sens = D.inv{D.val}.datareg.sensors;

M1 = D.inv{D.val}.datareg.toMNI;
[U, L, V] = svd(M1(1:3, 1:3));
M1(1:3,1:3) =U*V';

vol = forwinv_transform_vol(M1, vol);
sens = forwinv_transform_sens(M1, sens);

[vol, sens] = forwinv_prepare_vol_sens(vol, sens, 'channel', D.inv{D.val}.forward.channels);


%% ============ Select the data and convert to Fieldtrip struct

data = D.ftraw(0); % Convert to Fieldtrip without memory mapping

if D.ntrials > 1
    clb = D.conditions;
    ind = spm_input('Select trial',1, 'm', sprintf('%s|', clb{:}),1:D.ntrials);
    
    data.trial = data.trial(ind);
    data.time =  data.time(ind);
end

data = ft_timelockanalysis([], data);


%% =========== Configure and run Fieldtrip dipolefitting

cfg=[];
cfg.channel = D.inv{D.val}.forward.channels;
cfg.vol = vol;
cfg.inwardshift = 0;
cfg.grid.resolution=20;

if strcmp('EEG', modality)
    cfg.elec = sens;
else
    cfg.grad = sens;
end

cfg.latency  = 1e-3*spm_input('Time ([start end] in ms):', '+1', 'r', '', 2);

if spm_input('What to fit?','+1', 'm', 'dipole|pair', [0 1])
    cfg.numdipoles = 2;
    cfg.symmetry = 'x';
end

source = ft_dipolefitting(cfg, data);

%% =========== Plot the actual and the predicted scalp maps

cfg=[];
cfg.xparam='time';
cfg.xlim=[min(source.time) max(source.time)];
cfg.comment ='xlim';
cfg.commentpos='middlebottom';
cfg.electrodes='on';

if strcmp('EEG', modality)
    cfg.elec = sens;
    cfg.rotate = 0;
else
    cfg.grad = sens;
end

figure;
clf
subplot(1,2,1);
cfg.zparam='Vdata';
ft_topoplotER(cfg, source);
title('Data');
subplot(1,2,2);
cfg.zparam='Vmodel';
ft_topoplotER(cfg, source);
title('Model');

%% =========== Convert dipole position to MNI coordinates
Slocation = source.dip.pos;
Slocation(:,4) = 1;
Slocation = Slocation * (datareg.toMNI * inv(M1))';
Slocation = Slocation(:,1:3);


Nlocations = size(source.dip.pos,1);

%% =========== Display dipole locations using SPM's function
figure(Fgraph); clf

sdip= [];
sdip.n_seeds = 1;
sdip.n_dip   = Nlocations ;
sdip.Mtb     = 1;
sdip.j{1}    = zeros(3*Nlocations, 1);
sdip.loc{1}  = double(Slocation)';
spm_eeg_inv_ecd_DrawDip('Init', sdip)