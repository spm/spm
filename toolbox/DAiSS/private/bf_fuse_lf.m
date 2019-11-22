function [L, channels] = bf_fuse_lf(BF, modality)
% Prepares lead-fields to match channels in covariance
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_fuse_lf.m 7703 2019-11-22 12:06:29Z guillaume $

D = BF.data.D;

L = cell(size(BF.sources.L.(modality(1:3))));
channels = cell(length(BF.features.(modality).chanind), 1);

[sel1, sel2] = spm_match_str(D.chanlabels(BF.features.(modality).chanind),...
    BF.sources.channels.(modality(1:3)));

nlfcolumns = max(cellfun('size', BF.sources.L.(modality(1:3)), 2));

channels(sel1) = BF.sources.channels.(modality(1:3))(sel2);

% This is for MEG-EEG fusion
if isequal(modality, 'MEG') && length(sel1)<length(BF.features.('MEG').chanind)
    fuse = 1;
    [sel3, sel4] = spm_match_str(D.chanlabels(BF.features.('MEG').chanind),...
        BF.sources.channels.('EEG'));
    
    channels(sel3) = BF.sources.channels.EEG(sel4);
else
    fuse = 0;
end

spm_progress_bar('Init', numel(L), ['Preparing ' modality ' leadfields']); drawnow;
if numel(L) > 100, Ibar = floor(linspace(1, numel(L),100));
else Ibar = 1:numel(L); end

for i = 1:numel(L)
    if ~isnan(BF.sources.L.(modality(1:3)){i})
        L{i} = nan(length(BF.features.(modality).chanind), nlfcolumns);
        L{i}(sel1, :) = BF.sources.L.(modality(1:3)){i}(sel2, :);
        
        if fuse
            L{i}(sel3, :) = BF.sources.L.('EEG'){i}(sel4, :);
        end
    else
        L{i} = nan;
    end
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear')