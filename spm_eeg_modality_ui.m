function [mod, chanind]  = spm_eeg_modality_ui(D, scalp, planar)
% Attempts to determine the main modality of an MEEG object.
% If confused, asks the user.
%
% FORMAT [mod, chanind]  = spm_eeg_modality_ui(D, scalp, planar)
%
% D          - MEEG object
% scalp      - only look at scalp modalities [default: false]
% planar     - distinguish between MEG planar and other MEG types [default: false]
%
% modality   - the chosen modality
% chanind    - indices of the corresponding channels
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_modality_ui.m 5171 2013-01-03 15:52:31Z vladimir $

if nargin == 1
    scalp = false;
elseif nargin == 2
    planar = false;
end

[mod, list] = modality(D, scalp, planar);

if strcmpi(mod, 'Multimodal')
    qstr = [];
    for i = 1:numel(list)
        if ~isempty(qstr)
            qstr = [qstr '|'];
        end
        qstr = [qstr list{i}];
    end
    mod = list{spm_input('Which modality?','+1', 'm', qstr)};
end

if strcmpi(mod, 'MEG') && ~planar
    chanind = D.indchantype({'MEG', 'MEGPLANAR'});
else
    chanind = D.indchantype(mod);
end
