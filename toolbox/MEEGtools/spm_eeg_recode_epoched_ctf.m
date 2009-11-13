function D = spm_eeg_recode_epoched_ctf(S)
% Temporary solution for using trial labels in epoched CTF dataset
% FORMAT  D = spm_eeg_recode_epoched_ctf(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%          .D       - converted epoched CTF dataset
%
% Output:
% D                 - MEEG object relabeled trials (also saved to disk)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_recode_epoched_ctf.m 3566 2009-11-13 12:37:38Z vladimir $

SVNrev = '$Rev: 3566 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName', 'CTF trial recode'); spm('Pointer','Watch');

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

spm_progress_bar('Init', D.ntrials, 'Trials processed'); drawnow;
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
else Ibar = [1:D.ntrials]; end

for i = 1:D.ntrials
    ev = D.events(i);
    types = {ev{1}(:).type};
    values = {ev{1}(:).value};
    
    stimchanev = union(strmatch('UPPT', types), strmatch('FIL_', types), strmatch('trial', types, 'exact'));
    
    otherev = setdiff(1:numel(types), stimchanev);
    
    if length(otherev) > 1
        otherev = otherev(cellfun(@isempty, values(otherev)));
    end
    
    if length(otherev) ~= 1
        otherev = otherev(1);
        warning(sprintf('There are several events in trial %d that can possibly define trial labels. Using the first one.', i));
    end
    
    if ~isempty(otherev)
        D = conditions(D, i, types{otherev});
    else
        warning(sprintf('There are no events in trial %d that can possibly define trial labels. Skipping.', i));
    end
    
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

D = D.history(mfilename, S);

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','CTF trial recode: done'); spm('Pointer','Arrow');
