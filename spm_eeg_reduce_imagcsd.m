function res = spm_eeg_reduce_imagcsd(S)
% Plugin for data reduction based on the imaginary part of CSD
% with a reference chhannel
% FORMAT res = spm_eeg_reduce_imagcsd(S)
%
% S                     - input structure
% fields of S:
%   
%
% Output:
%  res -
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided:
%      montage struct implementing projection to PCA subspace
%______________________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_reduce_imagcsd.m 6852 2016-08-01 12:48:44Z vladimir $


if nargin == 0
    
    origchan = cfg_branch;
    origchan.tag = 'origchan';
    origchan.name = 'Channels to reduce';
    origchan.val = {spm_cfg_eeg_channel_selector};
    
    refchan = cfg_branch;
    refchan.tag = 'refchan';
    refchan.name = 'Reference channels';
    refchan.val = {spm_cfg_eeg_channel_selector};    
    
    
    outlabel = cfg_entry;
    outlabel.tag = 'outlabel';
    outlabel.name = 'Output channel label';
    outlabel.strtype = 's';
    outlabel.num = [1 Inf];
    outlabel.help = {'Label for the output channel(s).'};
    
    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency band of interest';
    foi.strtype = 'r';
    foi.num = [1 2];
    foi.val = {[0 Inf]};
    foi.help = {'Frequency window to optimize for'};       
    
    chanset = cfg_branch;
    chanset.tag = 'chanset';
    chanset.name = 'Set';
    chanset.val = {origchan, refchan, outlabel, foi};
    
    chansets = cfg_repeat;
    chansets.tag = 'chansets';
    chansets.name = 'Channel sets';
    chansets.values = {chanset};
    chansets.num = [1 Inf];
    chansets.val = {chanset};
    
    
    imagcsd = cfg_branch;
    imagcsd.tag = 'imagcsd';
    imagcsd.name = 'Imaginary part of CSD';
    imagcsd.val = {chansets};
    
    res = imagcsd;
    
    return
end

flag_tbx = license('checkout','signal_toolbox') && ~isempty(ver('signal'));
if flag_tbx
    taper = 'dpss';
else
    taper = 'sine';
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
    
    origind = setdiff(D.selectchannels(spm_cfg_eeg_channel_selector(S.chanset(i).origchan.channels)), badind);
    refind = setdiff(D.selectchannels(spm_cfg_eeg_channel_selector(S.chanset(i).refchan.channels)), badind);     
    
    if length(refind)~=1
        error(sprintf('One reference chhannel is necessary, found %d.', length(refind)));
    end
    
    Y  = D(origind, :, :);   
    Yr = spm_squeeze(D(refind, :, :), 1)';    
    
    if isequal(D.type, 'continuous')
        fs   = floor(D.fsample);
        ind = repmat((1:fs:D.nsamples)', 1, fs);
        ind = ind + repmat(0:(fs-1), size(ind, 1), 1);
        
        Yr = Yr(ind);
        
        time = (0:(fs-1))/fs;
    else
        time = D.time;
    end
    
    spm_progress_bar('Set','ylabel','computing CSD...');
    
    foi =  S.chanset(i).foi;
    foi(isinf(foi)) = 0.5*D.fsample;
    
    fY =[];
    for c = 1:size(Y, 1)
        Yc = spm_squeeze(Y(c, :, :), 1)';
        if size(Yc, 1) == 1
            Yc = Yc(ind);
        end        
        [spectrum,ntaper,freqoi] = ft_specest_mtmfft(Yc, time, 'taper', taper, 'freqoi', mean(foi), 'tapsmofrq', 0.5*diff(foi), 'verbose', 0);
        fY = [fY spectrum(:)];
    end
    
    fYr = ft_specest_mtmfft(Yr, time, 'taper', taper, 'freqoi', mean(foi), 'tapsmofrq', 0.5*diff(foi), 'verbose', 0);
    
    fYr = fYr(:);
    
    csd = mean(fY.*repmat(conj(fYr), 1, size(fY, 2)));
    
    mom = imag(csd)./norm(imag(csd));
    
    
    montage.labelnew{end+1, 1} = S.chanset(i).outlabel;
    
    montage.tra(end+1, end)  = 0;
    montage.tra(end, origind) = mom;
    
    montage.chantypenew{end+1}='LFP';
    
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
