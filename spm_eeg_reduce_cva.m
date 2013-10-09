function res = spm_eeg_reduce_cva(S)
% Plugin for data reduction using PCA
% FORMAT res = spm_eeg_reduce_pca(S)
%
% S                     - input structure
% fields of S:
%    S.ncomp            - number of PCA components
%
% Output:
%  res -
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided:
%      montage struct implementing projection to PCA subspace
%______________________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_reduce_cva.m 5675 2013-10-09 14:27:17Z vladimir $


if nargin == 0
    
    cvachan = cfg_branch;
    cvachan.tag = 'cvachan';
    cvachan.name = 'Channels to reduce';
    cvachan.val = {spm_cfg_eeg_channel_selector};
    
    refchan = cfg_branch;
    refchan.tag = 'refchan';
    refchan.name = 'Reference channels';
    refchan.val = {spm_cfg_eeg_channel_selector};
    
    ncomp = cfg_entry;
    ncomp.tag = 'ncomp';
    ncomp.name = 'Number of components';
    ncomp.strtype = 'n';
    ncomp.num = [1 1];
    ncomp.val = {1};
    ncomp.help = {'Number of components to retain'};
    
    
    outlabel = cfg_entry;
    outlabel.tag = 'outlabel';
    outlabel.name = 'Output channel label';
    outlabel.strtype = 's';
    outlabel.num = [1 Inf];
    outlabel.help = {'Label for the output channel(s).',...
        'Numbers are added for multiple channels.'};
    
    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency band of interest';
    foi.strtype = 'r';
    foi.num = [1 2];
    foi.val = {[0 Inf]};
    foi.help = {'Frequency window to optimize for'};
    
    tshiftwin = cfg_entry;
    tshiftwin.tag = 'tshiftwin';
    tshiftwin.name = 'Time shift window';
    tshiftwin.strtype = 'r';
    tshiftwin.num = [1 2];
    tshiftwin.val = {[0 0]};
    tshiftwin.help = {'Time shift window (ms). [0 0] - instantaneous'};
    
    tshiftres = cfg_entry;
    tshiftres.tag = 'tshiftres';
    tshiftres.name = 'Time shift resolution';
    tshiftres.strtype = 'r';
    tshiftres.num = [1 1];
    tshiftres.val = {5};
    tshiftres.help = {'Time shift resolution (ms)'};
    
    chanset = cfg_branch;
    chanset.tag = 'chanset';
    chanset.name = 'Set';
    chanset.val = {cvachan, refchan, ncomp, outlabel, foi, tshiftwin, tshiftres};
    
    chansets = cfg_repeat;
    chansets.tag = 'chansets';
    chansets.name = 'Channel sets';
    chansets.values = {chanset};
    chansets.num = [1 Inf];
    chansets.val = {chanset};
    
    
    cva = cfg_branch;
    cva.tag = 'cva';
    cva.name = 'CVA';
    cva.val = {chansets};
    
    res = cva;
    
    return
end

D = S.D;

nsets = numel(S.chanset);
badind = D.badchannels;

% Assuming projecting to columns
montage = [];
montage.labelorg    = D.chanlabels;
montage.labelnew    = {};
montage.chantypenew = {};
montage.tra         = zeros(0, D.nchannels);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nsets, 'Channel sets processed'); drawnow;
if nsets > 100, Ibar = floor(linspace(1, nsets,100));
else Ibar = 1:nsets; end

for i = 1:nsets
    
    spm_progress_bar('Set','ylabel','preparing data...');
    
    cvaind = setdiff(D.selectchannels(spm_cfg_eeg_channel_selector(S.chanset(i).cvachan.channels)), badind);
    refind = setdiff(D.selectchannels(spm_cfg_eeg_channel_selector(S.chanset(i).refchan.channels)), badind);
    
    if any(S.chanset(i).tshiftwin)
        tshiftind = S.chanset(i).tshiftwin(1):S.chanset(i).tshiftres:S.chanset(i).tshiftwin(2);
        tshiftind = [0 round(1e-3*D.fsample*tshiftind)];
        tshiftind = repmat(1:D.nsamples, length(tshiftind), 1)+repmat(tshiftind(:), 1, D.nsamples);
        tshiftind = tshiftind(:, all(tshiftind>0 & tshiftind<=D.nsamples));
    else
        tshiftind = repmat(1:D.nsamples, 2, 1);
    end
    
    nshift = size(tshiftind, 1)-1;
    nt     = size(tshiftind, 2);
    
    Y  = zeros(length(cvaind), nt, D.ntrials);
    Yr = zeros(length(refind)*nshift, nt, D.ntrials);
    
    for j = 1:D.ntrials
        cY  = D(cvaind, :, j);
        cYr = D(refind, :, j);
        
        if S.chanset(i).foi(1) > 0
            cY  = ft_preproc_highpassfilter(cY, D.fsample, S.chanset(i).foi(1),...
                5, 'but', 'twopass', 'reduce');
            cYr = ft_preproc_highpassfilter(cYr, D.fsample, S.chanset(i).foi(1),...
                5, 'but', 'twopass', 'reduce');
        end
        
        if isfinite(S.chanset(i).foi(2))
            cY  = ft_preproc_lowpassfilter(cY, D.fsample, S.chanset(i).foi(2),...
                5, 'but', 'twopass', 'reduce');
            cYr = ft_preproc_lowpassfilter(cYr, D.fsample, S.chanset(i).foi(2),...
                5, 'but', 'twopass', 'reduce');
        end
        
        Y(:, :, j) = cY(:, tshiftind(1, :));
        
        for k = 1:size(cYr, 1)
            Yr(((k-1)*nshift+1):k*nshift, :, j) = reshape(cYr(k, tshiftind(2:end, :)'),[], nshift)';
        end
    end
    
    Y  = reshape(Y,  size(Y, 1), []);
    Yr = reshape(Yr, size(Yr, 1), []);
    
    if D.fsample>2.5*S.chanset(i).foi(2)
        dec = floor(D.fsample/(2.5*S.chanset(i).foi(2)));
        Y   =  Y(:, 1:dec:end);
        Yr  = Yr(:, 1:dec:end);
    end
    
    spm_progress_bar('Set','ylabel','running CVA...');
    
    CVA = spm_cva(Y', Yr');
    
    ncomp = min(S.chanset(i).ncomp, size(CVA.V, 2));
    for j = 1:ncomp
        if ncomp == 1
            montage.labelnew{end+1, 1} = S.chanset(i).outlabel;
        else
            montage.labelnew{end+1, 1} = [S.chanset(i).outlabel num2str(j)];
        end
        
        montage.tra(end+1, end)  = 0;
        montage.tra(end, cvaind) = CVA.V(:, j)';
        
        montage.chantypenew{end+1}='LFP';
    end
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

if ~isempty(S.chanind)
    montage.labelnew = [montage.labelnew; D.chanlabels(S.chanind)'];
    I = eye(D.nchannels);
    montage.tra = [montage.tra; I(S.chanind, :)];
    montage.chantypenew = [montage.chantypenew, D.chantype(S.chanind)];
end

res = montage;
