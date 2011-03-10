function Dsource = spm_eeg_ft_beamformer_source(S)
% Extract source activity from specific points in the brain using LCMV
% beamformer
%
% FORMAT Dsource = spm_eeg_ft_beamformer_source(S)
%
% S            - struct (optional)
% (optional) fields of S:
% S.D          - meeg object or filename
% S.conditions - cell array of strings - condition labels to use
% S.sources    - specification of source locations. Struct with two fields:
%                .label - cell array of strings - labels for the sources
%                .pos - N x 3 matrix of MNI positions
% S.lambda     -  regularization constant (either numeric or string like '10%')
% S.outfile    - output file name (default 'B' + input name)
% S.appendchannels - cell array of strings - labels of channels from the
%                input dataset to be appended to the source data
%
% Output:
% Dsource      - MEEG object containing the source data and (optionally)
%                appended channel data.
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Robert Oostenveld
% $Id: spm_eeg_ft_beamformer_source.m 4238 2011-03-10 19:48:37Z vladimir $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'Beamformer source activity extraction',0);

if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select(1, 'mat', 'Select MEEG mat file');
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

if ~isfield(S, 'sources')
    if spm_input('Regional sources?','+1','yes|no',[1 0], 1);
        S.sources.pos = [];
        S.sources.label = {};
        while isempty(S.sources.pos) || spm_input('Add another source?','+1','yes|no',[1 0], 1);
            S.sources.label = [S.sources.label {spm_input('Source label', '+1', 's')}];
            S.sources.pos = [S.sources.pos;  spm_input('Source MNI coordinates', '+1', 'r', '', 3)'];
        end
    else
        S.sources = spm_eeg_dipoles_ui;
        S.sources.pos = S.sources.pnt;
    end
end

if ~isfield(S.sources, 'ori')
    if ~isfield(S, 'voi')
        if spm_input('Define VOI?','+1','yes|no',[1 0], 1);
            S.voi.radius =   spm_input('VOI radius (mm)', '+1', 'r', '10', 1);
            S.voi.resolution = spm_input('Resolution (mm)', '+1', 'r', '2', 1);
        else
            S.voi = 'no';
        end
    end
end


if ~isfield(S, 'lambda')
    S.lambda = spm_input('lambda (regularization)', '+1', 's',  '0.01%');
end

if ~isfield(S, 'outfile')
    S.outfile = spm_input('Output file name', '+1', 's', ['B' D.fname]);
end

modality = D.modality(1, 1);
if isequal(modality, 'Multimodal')
    if isfield(S, 'modality') && ~isempty(S.modality)
        modality = S.modality;
    else
        modality = spm_eeg_modality_ui(D, 1, 1);
    end
end

%% ============ Select the data and convert to Fieldtrip struct
if ~isfield(S, 'conditions')
    clb = D.condlist;
    
    if numel(clb) > 1
        
        [selection, ok]= listdlg('ListString', clb, 'SelectionMode', 'multiple' ,'Name', 'Select conditions' , 'ListSize', [400 300]);
        
        if ~ok
            return;
        end
    else
        selection = 1;
    end
    S.conditions = clb(selection);
end
%%
trialind = D.pickconditions(S.conditions);
%%
if isempty(trialind)
    error('No data was selected.');
end

%%
if ~isfield(S, 'appendchannels')
    if spm_input('Append other channels?','+1','yes|no',[1 0])
        selection = listdlg('ListString', D.chanlabels, 'SelectionMode', 'multiple' ,'Name', 'Select channels' , 'ListSize', [400 300]);
        if ~isempty(selection)
            S.appendchannels = D.chanlabels(selection);
        else
            S.appendchannels = {};
        end
    else
        S.appendchannels = {};
    end
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

for m = 1:numel(D.inv{D.val}.forward)
    if strncmp(modality, D.inv{D.val}.forward(m).modality, 3)
        vol  = D.inv{D.val}.forward(m).vol;
        if isa(vol, 'char')
            vol = ft_read_vol(vol);
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

vol = ft_transform_vol(M1, vol);
sens = ft_transform_sens(M1, sens);

pos = spm_eeg_inv_transform_points(M1*datareg.fromMNI, S.sources.pos);

channel = D.chanlabels(setdiff(meegchannels(D, modality), D.badchannels));

[vol, sens] = ft_prepare_vol_sens(vol, sens, 'channel', channel);

%%
spm('Pointer', 'Watch');drawnow;
%%

data = D.ftraw(0);
data.trial = data.trial(trialind);
data.time = data.time(trialind);

cfg = [];
cfg.channel = D.chanlabels(setdiff(D.meegchannels(modality), D.badchannels))';
cfg.covariance = 'yes';
cfg.keeptrials = 'no';
timelock1 = ft_timelockanalysis(cfg, data);
cfg.keeptrials = 'yes';
timelock2 = ft_timelockanalysis(cfg, data);
%%

nsources = numel(S.sources.label);

cfg = [];

if ismember(modality, {'MEG', 'MEGPLANAR'})
    cfg.reducerank = 2;
end

if ~isfield(S, 'voi') || isequal(S.voi, 'no')
    nvoi = 0;
    cfg.grid.pos     = pos;
    if isfield(S.sources, 'ori')
        cfg.grid.mom    = S.sources.ori;
    else
        cfg.lcmv.fixedori = 'yes';
    end
else
    cfg.lcmv.fixedori = 'yes';
    vec = -S.voi.radius:S.voi.resolution:S.voi.radius;
    [X, Y, Z]  = ndgrid(vec, vec, vec);
    sphere   = [X(:) Y(:) Z(:)];
    sphere(sqrt(X(:).^2 + Y(:).^2 + Z(:).^2)>S.voi.radius, :) = [];
    nvoi = size(sphere, 1);
    cfg.grid.pos = [];
    for s = 1:size(pos, 1)
        cfg.grid.pos = [cfg.grid.pos; sphere+repmat(pos(s, :), nvoi, 1)];
    end
end

cfg.grad = sens;
cfg.inwardshift = -30;
cfg.vol = vol;
cfg.channel = modality;
cfg.method = 'lcmv';
cfg.keepfilter = 'yes';
cfg.keepleadfield = 'yes';
cfg.lambda =  S.lambda;
source1 = ft_sourceanalysis(cfg, timelock1);

if isfield(S, 'makecorrimage') &&  S.makecorrimage
    if size(pos, 1) ~= numel(source1.avg.filter)
        error('Can only make correlation images for point sources');
    end
    
    res = mkdir(pwd, 'corrimages');
         
    cfg.grid = [];
    cfg.grid.xgrid = -90:10:90;
    cfg.grid.ygrid = -120:10:100;
    cfg.grid.zgrid = -70:10:110;
    cfg.inwardshift = -10;
    
    fsource = ft_sourceanalysis(cfg, timelock1);
    
    
    sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
    
    pow = nan(size(fsource.avg.pow));
    for p = 1:numel(source1.avg.filter)
        for q = 1:length(fsource.inside);
            cc = corrcoef(source1.avg.filter{p}, fsource.avg.filter{fsource.inside(q)});
            pow(fsource.inside(q)) = (cc(1,2))^2;
            fprintf('Correlation image %d/%d\n', q, length(fsource.inside));
        end
        
        fsource.pow = pow;
        
        cfg1 = [];
        cfg1.sourceunits   = 'mm';
        cfg1.parameter = 'pow';
        cfg1.downsample = 1;
        sourceint = ft_sourceinterpolate(cfg1, fsource, sMRI);
        %%           
        outvol = spm_vol(sMRI);
        outvol.dt(1) = spm_type('float32');
        
        if length(trialind) == 1
            suff = num2str(trialind);
        else
            suff = '';
        end
        
        outvol.fname= fullfile(pwd, 'corrimages', ['corrimg_' spm_str_manip(D.fname, 'r') '_' S.sources.label{p} '_' suff '.nii']);
        outvol = spm_create_vol(outvol);
        spm_write_vol(outvol, sourceint.pow);
    end
end

cfg = [];
cfg.inwardshift = -30;
cfg.vol = vol;
cfg.grad = sens;
cfg.grid = ft_source2grid(source1);
cfg.channel = modality;
cfg.lambda =  S.lambda;
cfg.rawtrial = 'yes';
source2 = ft_sourceanalysis(cfg, timelock2);

crosstalk = [];
for i = 1:nsources
    for j = 1:nsources
        cc = corrcoef([cfg.grid.filter{i}' cfg.grid.leadfield{j}]).^2;
        crosstalk(i, j) = cc(1, 2);
    end
end
%%
if nsources > 1    
    Fgraph = spm_figure('GetWin','Graphics');
    colormap(gray)
    figure(Fgraph)
    clf
    
    % images
    %----------------------------------------------------------------------
    subplot(2, 1 , 1)
    imagesc(crosstalk)
    title('Source crosstalk (R^2)','FontSize',13)
    set(gca,'YTick',[1:nsources],'YTickLabel',S.sources.label,'FontSize',8)
    set(gca,'XTick',[1:nsources],'XTickLabel',S.sources.label,'FontSize',8)
    xlabel('from','FontSize',10)
    ylabel('to','FontSize',10)
    axis square
    
    % table
    %----------------------------------------------------------------------
    subplot(2,1,2)
    text(1/6,1/2,num2str(crosstalk,' %-8.2f'),'FontSize',16)
    axis off,axis square
end
%%
sourcedata=[];
sourcedata.trial=zeros(size(timelock2.trial, 1), nsources, length(timelock2.time));

for i=1:length(source2.trial)
    disp(['Extracting source data trial ' num2str(i) '/'  num2str(length(source2.trial))]);
    for j = 1:nsources
        if nvoi>0
            y = cat(1, source2.trial(i).mom{(j-1)*nvoi+[1:nvoi]});            
        else
            y       = source2.trial(i).mom{j};
        end
        
        % compute regional response in terms of first eigenvariate
        %-----------------------------------------------------------------------
        [m n]   = size(y);
        if m > n && n>1
            [v s v] = svd(y'*y);
            Y       = v(:,1);
        elseif m>1
            [u s u] = svd(y*y');
            u       = u(:,1);
            Y       = y'*u;
        else
            Y = y;
        end
        
        sourcedata.trial(i, j, :)= Y;
    end
end


sourcedata.time = timelock2.time;
sourcedata.dimord = 'rpt_chan_time';
sourcedata.label = S.sources.label;
sourcedata.avg = [];
%%
if ~isempty(S.appendchannels)
    cfg = [];
    cfg.channel= S.appendchannels;
    cfg.keeptrials = 'yes';
    externdata = ft_timelockanalysis(cfg, data);
    sourcedata = ft_appenddata([], sourcedata, externdata);
end
%%

Dsource = spm_eeg_ft2spm(sourcedata, S.outfile);

Dsource = chantype(Dsource, 1:length(S.sources.label), 'LFP');

Dsource = conditions(Dsource, [], D.conditions(trialind));

Dsource = history(Dsource, 'spm_eeg_ft_beamformer_source', S);

Dsource.crosstalk = crosstalk;

save(Dsource);

spm('Pointer', 'Arrow');drawnow;
