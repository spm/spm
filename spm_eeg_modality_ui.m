function [mod, chanind]  = spm_eeg_modality_ui(D, scalp, planar)
% spm_eeg_modality_ui - tries to determine the main modality
% of an meeg object. If confused, asks the user.
%
% FORMAT [mod, chanind]  = spm_eeg_modality_ui(D, scalp, planar)
%
% modality - the chosen modality
% chanind  - indices of the corresponding channels
% planar   - distinguish between MEG planar and other MEG types
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_modality_ui.m 2720 2009-02-09 19:50:46Z vladimir $

if nargin == 1
    scalp = 0;
elseif nargin ==2
    planar = 0;
end

[mod, list] = modality(D, scalp, planar);

if isequal(mod, 'Multimodal')
    qstr = [];
    for i = 1:numel(list)
        if ~isempty(qstr)
            qstr = [qstr '|'];
        end
        qstr = [qstr list{i}];
    end
    mod = list{spm_input('Which modality?','+1', 'm', qstr)};
end

if isequal(mod, 'MEG') && planar
    chanind = setdiff(strmatch('MEG', D.chantype) , strmatch('MEGPLANAR', D.chantype));
else
    chanind = strmatch(mod, D.chantype);
end