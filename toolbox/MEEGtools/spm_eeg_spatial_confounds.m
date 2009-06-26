function D = spm_eeg_spatial_confounds(S)
% This function defines spatial confounds and adds them to MEEG dataset. 
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_spatial_confounds.m 3228 2009-06-26 17:43:19Z vladimir $


SVNrev = '$Rev: 3228 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Define spatial confounds');


%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

if ~isfield(S, 'method')
    S.method = spm_input('How to define?','+1','SVD|SPMEEG|BESA|Clear', strvcat('SVD', 'SPMEEG', 'BESA', 'Clear'));
end

switch upper(S.method)
    case 'BESA'
        if ~isfield(S, 'conffile')
            S.conffile = spm_select(1, '\.bsa$', 'Select BESA *.bsa file');
        end
        sconf = spm_eeg_read_bsa(S.conffile );
        D = sconfounds(D, sconf);
    case 'SVD'
        if ~isfield(S, 'timewin')
            S.timewin  = spm_input('PST window (ms)', '+1', 'r', '', 2)/1000;
        end
        cl = D.condlist;
        svdinput = [];
        for i = 1:numel(cl)
            svdinput = [svdinput mean(D.selectdata(D.chanlabels(setdiff(D.meegchannels, D.badchannels)), S.timewin, cl{i}), 3)];
        end
        U = spm_svd(svdinput);

        if ~isfield(S, 'ncomp')
            S.ncomp = spm_input('How many components?', '+1', 'n', '1', 1);
        end
        [sel1, sel2] = spm_match_str(D.chanlabels(D.meegchannels), D.chanlabels(setdiff(D.meegchannels, D.badchannels)));
        sconf = [];
        sconf.label = D.chanlabels(D.meegchannels);
        sconf.coeff = nan(length(sconf.label), S.ncomp);
        sconf.coeff(sel1, :) = U(sel2, 1:S.ncomp);
        sconf.bad = ones(length(sconf.label), 1);
        sconf.bad(sel1, :) = 0;
        D = sconfounds(D, sconf);
    case 'SPMEEG'
        if ~isfield(S, 'conffile')
            S.conffile =  spm_select(1, 'mat', 'Select M/EEG mat file with sconfounds');
        end
        Ds = spm_eeg_load(S.conffile);
        sconf = getfield(Ds, 'sconfounds');
        D = sconfounds(D, sconf);
    case 'CLEAR'
        D = rmfield(D, 'sconfounds');
end

D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName', 'Define spatial confounds: done');
