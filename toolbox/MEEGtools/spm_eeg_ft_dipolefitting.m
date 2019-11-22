% Demo script for dipole fitting using Fieldtrip
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%

% Vladimir Litvak
% $Id: spm_eeg_ft_dipolefitting.m 7702 2019-11-22 11:32:26Z guillaume $

[Finter,Fgraph] = spm('FnUIsetup','Fieldtrip dipole fitting', 0);
%%

%% ============ Load SPM EEG file and verify consistency

D = spm_eeg_load;

[D, ok] = check(D, 'sensfid');

if ~ok
    if check(D, 'basic')
        errordlg(['The requested file is not ready for source reconstruction.'...
            'Use prep to specify sensors and fiducials.']);
    else
        errordlg('The meeg file is corrupt or incomplete');
    end
    return
end

modality = spm_eeg_modality_ui(D, 1, 1);

chanind = indchantype(D, modality, 'GOOD');

modality = modality(1:3);

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

volsens = spm_eeg_inv_get_vol_sens(D, [], [], [], modality);

vol  = volsens.(modality).vol;
sens = volsens.(modality).sens;

if isa(vol, 'char')
    vol = ft_read_headmodel(vol);
end




%% ============ Select the data and convert to Fieldtrip struct
if D.ntrials > 1
    clb = D.conditions;
    ind = spm_input('Select trial',1, 'm', sprintf('%s|', clb{:}),1:D.ntrials);    
else
    ind = 1;
end

data = D.fttimelock(chanind, ':',  ind);

try, data = rmfield(data, 'elec'); end
try, data = rmfield(data, 'grad'); end
%% =========== Configure and run Fieldtrip dipolefitting

cfg=[];
cfg.vol = vol;
cfg.grid.resolution=20*1e-3;

if strcmp('EEG', modality)
    cfg.elec = sens;
    reducerank = 3;
else
    cfg.grad = sens;
    reducerank = 2;
end

cfg.latency  = 1e-3*spm_input('Time ([start end] in ms):', '+1', 'r', '', 2);

if spm_input('What to fit?','+1', 'm', 'dipole|pair', [0 1])
    cfg.numdipoles = 2;
    cfg.symmetry = 'x';
end

source = ft_dipolefitting(cfg, data);

%% =========== Plot the actual and the predicted scalp maps

cfg=[];
cfg.xlim=[min(source.time) max(source.time)];
cfg.comment ='xlim';
cfg.commentpos='middlebottom';
cfg.marker='on';

if strcmp('EEG', modality)
    cfg.elec = sens;
    cfg.rotate = 0;
else
    cfg.grad = sens;
end

figure;
clf
subplot(1,2,1);
cfg.parameter ='Vdata';
ft_topoplotER(cfg, rmfield(source, 'dip'));
title('Data');
subplot(1,2,2);
cfg.parameter ='Vmodel';
ft_topoplotER(cfg, rmfield(source, 'dip'));
title('Model');

%% =========== Convert dipole position to MNI coordinates
Slocation = source.dip.pos;
Slocation(:,4) = 1;
Slocation = Slocation * (volsens.transforms.toMNI)';
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