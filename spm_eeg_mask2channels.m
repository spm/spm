function chanind = spm_eeg_mask2channels(D, mask)
% Make a list of channel labels based on scalp mask
% FORMAT chanind = spm_eeg_firstlevel(D, mask)
%
% D - M/EEG object (or filename)
% mask - mask (numeric array, nifti object or image file name)
%        if the mask is 3D channels in all the blobs will be returned
%
% Output:
% chanind - indices of channels in D which correspond to blobs in the mask
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_mask2channels.m 3181 2009-06-04 07:52:57Z vladimir $

%-Check inputs
%--------------------------------------------------------------------------

if nargin < 2
    error('At least two inputs are required.');
end

if isa(D, 'char')
    D = spm_eeg_load(D);
end

if isa(mask, 'char')
    mask = nifti(mask);
end

if isa(mask, 'nifti')
    mask = mask.dat(:, :, :);
end

if  ndims(mask) == 3
    mask = squeeze(any(mask, 3));
end

if any(diff(size(mask)))
    error('The mask should be square');
else
    n = size(mask, 1);
end

%-Get channel indices and coordinates
%--------------------------------------------------------------------------
[Cel, Cind] = spm_eeg_locate_channels(D, n, 1);

modality = spm_eeg_modality_ui(D, 1, 1);

goodchan = setdiff(D.meegchannels(modality), D.badchannels);

%-Find channels corresponding to blobs
%--------------------------------------------------------------------------
chanind = [];
for i = 1:length(Cind)
    if ismember(Cind(i), goodchan)
        val = mask(Cel(i, 1), Cel(i, 2));
        if ~isnan(val) && val~=0
            chanind = [chanind; Cind(i)];
        end
    end
end

