function ind = meegchannels(this, modality)
% Method for getting index vector of m/eeg channels for display
% FORMAT ind = meegchannels(this, modality)
% modality - (optional) - one of EEG, MEG (excluding planar), MEGPLANAR
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: meegchannels.m 2749 2009-02-16 11:30:15Z vladimir $

type = chantype(this);

if nargin == 1
    ind = unique([strmatch('EEG', type, 'exact'); strmatch('MEG', type); strmatch('REF', type); strmatch('LFP', type)]);
else
    switch modality
        case 'EEG'
            ind = strmatch('EEG', type, 'exact');
        case 'LFP'
            ind = strmatch('LFP', type, 'exact');
        case 'MEG'
            ind = sort([strmatch('MEGMAG', type, 'exact'); strmatch('MEGGRAD', type, 'exact'); strmatch('MEG', type, 'exact')]);
        case 'MEGPLANAR'
            ind = strmatch('MEGPLANAR', type, 'exact');
        otherwise
            error('Unsupported modality');
    end
end

ind = ind(:)'; % must be row to allow to use it as loop indices

