function spm_eeg_ft_beamformer_smoothness(S)
% Script for making power images using DICS beamformer. Requires as input
% an MEEG file where coregistration has been performed.
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL

% Vladimir Litvak
% $Id: spm_eeg_ft_beamformer_smoothness.m 3263 2009-07-10 11:56:55Z vladimir $

[Finter,Fgraph] = spm('FnUIsetup','Smoothness estimation for beamformer images', 0);
%%

%% ============ Load SPM EEG file and verify consistency
if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
    S.D = D;
end

if ischar(D)
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end

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

modality = spm_eeg_modality_ui(D, 1, 1);
channel = D.chanlabels(strmatch(modality, D.chantype))';
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

for m = 1:numel(D.inv{D.val}.forward)
    if strncmp(modality, D.inv{D.val}.forward(m).modality, 3)
        vol  = D.inv{D.val}.forward(m).vol;
        if isa(vol, 'char')
            vol = fileio_read_vol(vol);
        end
        datareg  = D.inv{D.val}.datareg(m);
    end
end

if isequal(modality, 'EEG')
    sens = datareg.sensors;
else
    % This is to make it possible to use the same 'inv' in multiple files
    sens = D.sensors('MEG');
end

M1 = datareg.toMNI;
[U, L, V] = svd(M1(1:3, 1:3));
M1(1:3,1:3) =U*V';

vol = forwinv_transform_vol(M1, vol);
sens = forwinv_transform_sens(M1, sens);


%% ============ Select the data and convert to Fieldtrip struct

clb = D.condlist;

if numel(clb) > 1

    [selection, ok]= listdlg('ListString', clb, 'SelectionMode', 'multiple' ,'Name', 'Select conditions' , 'ListSize', [400 300]);

    if ~ok
        return;
    end
else
    selection = 1;
end
%%
ind = D.pickconditions(clb(selection));
%%
data = D.ftraw(0);
data.trial = data.trial(ind);
data.time =  data.time(ind);
%%
if ~isfield(S, 'beamformer')
    S.beamformer = spm_input('Which beamformer?','+1', 'LCMV|DICS', strvcat('lcmv', 'dics'));
end
%%
refchan = [];

if strcmpi(S.beamformer, 'dics')
    if isfield(S, 'refchan') && ~isempty(S.refchan)
        refchan = S.refchan;
    end

    if ~isfield(S, 'centerfreq')
        S.centerfreq =  spm_input('Frequency (Hz):', '+1', 'r', '', 1);
    end

    if ~isfield(S, 'tapsmofrq')
        S.tapsmofrq  =  spm_input('Frequency window (Hz):', '+1', 'r', '', 1);
    end
end

cfg = [];
cfg.keeptrials = 'yes';
cfg.channel    = [channel; refchan];
%%
if ~isfield(S, 'timewindows')
    for i = 1:spm_input('Number of time windows:', '+1', 'r', '1', 1)
        S.timewindows{i} = spm_input('Time ([start end] in sec):', '+1', 'r', num2str([D.time(1) D.time(end)]), 2);
    end
end

for i = 1:numel(S.timewindows)
    cfg.latency = S.timewindows{i};
    timelock{i} = ft_timelockanalysis(cfg,data);
end

%%
if ~isfield(S, 'lambda')
    S.lambda = [num2str(spm_input('Regularization (%):', '+1', 'r', '0')) '%'];
end
%%
if strcmpi(S.beamformer, 'dics')
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'powandcsd';
    cfg.tapsmofrq = S.tapsmofrq;
    cfg.foilim    = [S.centerfreq S.centerfreq];
    cfg.keeptrials = 'yes';

    for i = 1:numel(timelock)
        freq{i} = ft_freqanalysis(cfg, timelock{i});
    end

    freqall = freq{1};
    if length(freq) > 1
        for i = 2:length(freq)
            freqall.powspctrm = cat(1,  freqall.powspctrm, freq{i}.powspctrm);
            freqall.crsspctrm = cat(1,  freqall.crsspctrm, freq{i}.crsspctrm);
            freqall.cumsumcnt = cat(1,  freqall.cumsumcnt, freq{i}.cumsumcnt);
            freqall.cumtapcnt = cat(1,  freqall.cumtapcnt, freq{i}.cumtapcnt);
        end
    end
else
    if numel(timelock)>1
        for i = 1:numel(timelock)
            timelock{i}.time = timelock{1}.time;
        end
        timelock = ft_appenddata([], timelock{:});
    else
        timelock = timelock{1};
    end

    cfg = [];
    cfg.channel = modality;
    cfg.covariance = 'yes';
    cfg.covariancewindow = 'maxperlength';
    cfg.keeptrials = 'no';
    timelock = ft_timelockanalysis(cfg, timelock);
end
%%

gridres = 5;

cfg                       = [];

if strcmp('EEG', modality)
    cfg.elec = sens;
else
    cfg.grad = sens;
    cfg.reducerank            = 2;
end

cfg.channel     = D.chanlabels(D.meegchannels(modality));
cfg.vol         = vol;

cfg.grid.xgrid  = -90:gridres:90;
cfg.grid.ygrid  = -120:gridres:100;
cfg.grid.zgrid  = -70:gridres:110;
cfg.inwardshift = 0;
grid            = ft_prepare_leadfield(cfg);


cfg = [];

if strcmp('EEG', modality)
    cfg.elec = D.inv{D.val}.datareg.sensors;
else
    cfg.grad = D.sensors('MEG');
    cfg.reducerank = 2;
end

cfg.channel = D.chanlabels(D.meegchannels(modality));

nchan = numel(cfg.channel);

cfg.keepfilter   = 'yes';
cfg.projectnoise = 'no';
cfg.grid         = grid;
cfg.vol          = vol;
cfg.lambda       = S.lambda;

if strcmpi(S.beamformer, 'dics')
    if ~isempty(refchan)
        cfg.refchan = refchan;
    end

    cfg.method       = 'dics';
    if isfield(S, 'fixedori') && S.fixedori
        cfg.dics.fixedori = 'yes';
    end
    cfg.dics.realfilter = 'yes';
    cfg.dics.powmethod  = 'trace';
    cfg.frequency    = S.centerfreq;
    filtsource   = ft_sourceanalysis(cfg, freqall);
elseif strcmpi(S.beamformer, 'lcmv')
    cfg.lcmv.fixedori = 'yes';
    filtsource   = ft_sourceanalysis(cfg, timelock);
end

filter  = filtsource.avg.filter; % use the filter computed in the previous step
%%
res  = mkdir(D.path, 'smoothness');
sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
rand('twister',sum(100*clock));
filenames = [];
nrep = 10;
for i = 1:nrep
    pow  = nan(1, size(filtsource.pos, 1));
    noise = randn(nchan, 1);
    for  j = 1:length(grid.inside)      
      filt = filter(grid.inside(j));
      filt = filt{1};      
      p    = filt*noise;
      if length(p)>1
          p = norm(p);
      end
      pow(grid.inside(j)) = p;
      fprintf('projecting noise %d/%d\n', j, length(grid.inside));
    end
    filtsource.pow = pow;
    cfg = [];
    cfg.parameter = 'pow';
    %cfg.interpmethod = 'nearest';
    cfg.downsample = 1;
    sourceint = ft_sourceinterpolate(cfg, filtsource, sMRI);

    outvol = spm_vol(sMRI);
    outvol.fname= fullfile(D.path, 'smoothness', ['img_' spm_str_manip(D.fname, 'r') '_rep' num2str(i) '.nii']);
    outvol = spm_create_vol(outvol);
    spm_write_vol(outvol, sourceint.pow);
    filenames = strvcat(filenames, outvol.fname);
end

outvol.fname= fullfile(D.path, 'smoothness', 'mask.nii');
spm_write_vol(outvol, ~isnan(sourceint.pow));

cd(fullfile(D.path, 'smoothness'));

[fwhm,VRpv] = spm_est_smoothness(filenames, outvol.fname, [nrep nrep-2]);

[resels_pervox]=spm_read_vols(VRpv);
total_resels=sum(resels_pervox(find(isfinite(resels_pervox))));
disp(fwhm);

