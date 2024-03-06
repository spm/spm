function BF = bf_stat_evoked_t(S)
% Compute t-stat accross trials for evoked response
% FORMAT BF = bf_stat_evoked_t(S)
%   S               - input structure
%  fields of S:
%   S.BF        - path to BF.mat file
%   S.act       - active timpeoint(ms) - 1 x 1 matrix      -Default: none
%   S.base      - base timpeoint(ms) - 1 x 1 matrix        -Default: none
%   S.condact   - active condition label - string          -Default: 'ALL'
%   S.condbase  - baseline condition label - string        -Default: 'ALL'
%   S.MNI       - flag to output in MNI space - logical    -Default: true
%   S.summary   - output summary statistic  - logical      -Default: false 
% Output:
%   BF               - path to BF.mat file
%__________________________________________________________________________

% Tim Tierney
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


%-Arguments
%--------------------------------------------------------------------------
if ~isfield(S, 'condact'),  S.condact  = 'ALL'; end
if ~isfield(S, 'condbase'), S.condbase = 'ALL'; end
if ~isfield(S, 'MNI'),      S.MNI      = true;  end
if ~isfield(S, 'summary'),  S.summary  = false; end

if S.summary
    error('Summary statistics not yet supported.');
end

%-Channels, labels and filters
%--------------------------------------------------------------------------
BF      = load(S.BF);
D       = BF.data.D;
W       = BF.inverse.MEG.W;
labs    = BF.inverse.MEG.channels;
labinds = indchannel(BF.data.D,labs);

%-Conditions
%--------------------------------------------------------------------------
if strcmp(S.condact,'ALL')
    acttrials = 1:size(D,3);
else
    acttrials = indtrial(D,S.condact);
end

if strcmp(S.condbase,'ALL')
    basetrials = 1:size(D,3);
else
    basetrials = indtrial(D,S.condbase);
end

%-Timepoints and trials
%--------------------------------------------------------------------------
[~,bind] = min(abs(S.base/1000-BF.data.D.time()));
[~,aind] = min(abs(S.act/1000-BF.data.D.time()));

ntrials = length(basetrials);
bY = squeeze(D(labinds,bind,basetrials));
aY = squeeze(D(labinds,aind,acttrials));

%-Sensor mean and trial-wise covariance
%--------------------------------------------------------------------------
del       = aY - bY;
sensorMu  = mean(del,2);
del       = del - repmat(sensorMu,1,size(del,2));
sensorC   = del * del' / (ntrials-1);
image.val = zeros(size(W,2),1);

%-Source t-stat accross trials
%--------------------------------------------------------------------------
for i=1:size(image.val)

    brainMu   = W{i} * sensorMu;
    sourceVar = W{i} * sensorC * W{i}';
    sourceSD  = sqrt(diag(sourceVar));
    sourceSE  = sourceSD / sqrt(ntrials);
    t3        = brainMu ./ sourceSE;

    % orientation of max SNR
    or(i,:)   = t3 ./ sqrt(sum(t3.^2));

    % new t-stat in derived orientation
    brainMu   = or(i,:) * W{i} * sensorMu;
    sourceVar = or(i,:) * W{i} * sensorC * (or(i,:)*W{i})';
    sourceSD  = sqrt(diag(sourceVar));
    sourceSE  = sourceSD / sqrt(ntrials);
    t         = brainMu ./ sourceSE;
    image.val(i) = t;

end

%-Source space
%--------------------------------------------------------------------------
sMRI       = BF.data.mesh.sMRI;
source     = BF.sources.grid;
source.pos = BF.sources.grid.allpos;
if S.MNI
    source = ft_transform_geometry(BF.data.transforms.toMNI, source);
    sMRI   = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
end

source.pow = image.val;
source.pow = source.pow(:);
pow        = source.pow;
source.pow = NaN(size(source.pos, 1), 1);
source.pow(source.inside) = pow;

%-Interpolate
%--------------------------------------------------------------------------
cfg = [];
cfg.parameter = 'pow';
cfg.downsample = 1;
cfg.showcallinfo = 'no';
cfg.interpmethod ='linear';

hdr = spm_vol(sMRI);
mri = [];
mri.dim = hdr.dim;
mri.anatomy = double(hdr.private.dat);
mri.hdr = hdr;
mri.transform = hdr.mat;
mri.unit = 'mm';
            
sourceint = ft_sourceinterpolate(cfg, source, mri);

%-Write output image
%--------------------------------------------------------------------------
Y = sourceint.pow;
outvol = spm_vol(sMRI);
outvol.dt(1) = spm_type('float32');
lab = ['t_' S.condact '_' num2str(S.act) '_' S.condbase '_' num2str(S.base)];
lab = [lab '_' spm_file(D.fname, 'basename')];
outvol.fname = fullfile(pwd, [lab '.nii']);
outvol = spm_create_vol(outvol);
spm_write_vol(outvol, Y);
