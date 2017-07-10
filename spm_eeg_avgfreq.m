function D = spm_eeg_avgfreq(S)
% Average a TF-dataset over frequency to get a time-domain dataset
% FORMAT D = spm_eeg_avgfreq(S)
%
% S        - input struct
%  fields of S:
%   D        - MEEG object or filename of M/EEG mat-file with epoched data
%   freqwin  - frequency window to average over [default: [-Inf, Inf]]
%   prefix   - prefix for the output file [default: 'P']
%
% Output:
% D        - MEEG object
%__________________________________________________________________________
% Copyright (C) 2012-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_avgfreq.m 7132 2017-07-10 16:22:58Z guillaume $

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Average over frequency'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix   = 'P';           end
if ~isfield(S, 'freqwin'),      S.freqwin  = [-Inf Inf];    end

D = spm_eeg_load(S.D);

if ~strncmpi(D.transformtype,'TF',2)
    error('This function only works on TF datasets.');
end

freqind = D.indfrequency(min(S.freqwin)):D.indfrequency(max(S.freqwin));
if isempty(freqind) || any(isnan(freqind))
    error('Selected frequency window is invalid.');
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dnew = clone(D, [S.prefix fname(D)], [D.nchannels D.nsamples D.ntrials]);

%-Averaged data
%--------------------------------------------------------------------------
spm_progress_bar('Init', D.ntrials, 'Trials processed');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end

for i = 1:D.ntrials
    
    Dnew(:, :, i) = spm_squeeze(mean(D(:, freqind, :, i), 2), 2);
    
    if any(Ibar == i), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

%-Save the new M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history(mfilename, S);
save(Dnew);

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','Average over frequency: done'); spm('Pointer','Arrow');
