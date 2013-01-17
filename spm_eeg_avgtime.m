function D = spm_eeg_avgtime(S)
% Average a TF-dataset over time to get a spectrum dataset
% FORMAT D = spm_eeg_avgtime(S)
%
% S        - input struct
%  fields of S:
%   D        - MEEG object or filename of M/EEG mat-file with epoched data  
%   timewin  - time window to average over (in PST ms)
%   prefix   - prefix for the output file (default - 'S')
%
%
% Output:
% D        - MEEG object 
%
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_avgtime.m 5190 2013-01-17 15:32:45Z vladimir $

SVNrev = '$Rev: 5190 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Average over time'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix   = 'P';           end
if ~isfield(S, 'timewin'),      S.timewin  = [-Inf Inf];    end

D = spm_eeg_load(S.D);

if ~strncmpi(D.transformtype,'TF',2);
    error('This function only works on TF datasets');
end

timeind = D.indsample(1e-3*(min(S.timewin))):D.indsample(1e-3*(max(S.timewin)));
if isempty(timeind) || any(isnan(timeind))
    error('Selected time window is invalid.');
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dnew = clone(D, [S.prefix fname(D)], [D.nchannels D.nfrequencies 1 D.ntrials]);

%-Averaged data
%--------------------------------------------------------------------------
spm_progress_bar('Init', D.ntrials, 'Trials processed');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end

for i = 1:D.ntrials
    
    Dnew(:, :, :, i) =  spm_squeeze(mean(D(:, :, timeind, i), 3), 3);
    
    if D.trialonset(i) ~= 0
        Dnew = trialonset(Dnew, i,  D.trialonset(i)+ mean(D.time([timeind(1) timeind(end)])));
    end
    
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end  %

Dnew = timeonset(Dnew, mean(D.time([timeind(1) timeind(end)])));

spm_progress_bar('Clear');

%-Save the new M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history(mfilename, S);
save(Dnew);

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','Average over time: done'); spm('Pointer','Arrow');
