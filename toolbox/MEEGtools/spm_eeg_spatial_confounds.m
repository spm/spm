% This script defines spatial confounds and adds them to MEEG dataset. This
% This functionality will be further developed in the future. For now spatial
% confounds are loaded from a BESA *.bsa file.
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_spatial_confounds.m 2347 2008-10-16 12:58:33Z vladimir $

D = spm_eeg_load;

if spm_input('How to define?','+1','SVD|BESA',[0 1], 0);
    sconf = spm_eeg_read_bsa(spm_select(1, '\.bsa$', 'Select BESA *.bsa file'));
else
    timewin = spm_input('PST window (ms)', '+1', 'r', '', 2)/1000;
    [U S V] = spm_svd(mean(D.selectdata(D.chanlabels(setdiff(D.meegchannels, D.badchannels)), timewin, []), 3));
    n = spm_input('How many components?', '+1', 'n', '1', 1);
    [sel1, sel2] = spm_match_str(D.chanlabels(D.meegchannels), D.chanlabels(setdiff(D.meegchannels, D.badchannels)));
    sconf = [];
    sconf.label = D.chanlabels(D.meegchannels);
    sconf.coeff = nan(length(sconf.label), n);
    sconf.coeff(sel1, :) = U(sel2, 1:n);
    sconf.bad = ones(length(sconf.label), 1);
    sconf.bad(sel1, :) = 0;    
end

D = sconfounds(D, sconf);

save(D);