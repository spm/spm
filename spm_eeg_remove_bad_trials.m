function D = spm_eeg_remove_bad_trials(S)
% Physically removes trials marked as bad from the dataset
% FORMAT D = spm_eeg_remove_bad_trials(S)
%
% S        - optional input struct
%  fields of S:
% D        - MEEG object or filename of M/EEG mat-file with epoched data
% prefix   - prefix for the output file (default - 'r')
%
% Output:
% D        - MEEG object (also written on disk)
%
% The function also changes the physical order of trials to conform to
% condlist.
%
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_remove_bad_trials.m 6526 2015-08-20 10:28:36Z vladimir $

SVNrev = '$Rev: 6526 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Remove bad trials'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix = 'r';              end

D = spm_eeg_load(S.D);

%-Check that there is any good data available
%--------------------------------------------------------------------------
if ntrials(D)==0 || isempty(indtrial(D, D.condlist, 'GOOD'))
    warning('No good trials were found. Nothing to do.');
    return;
end

cl   = D.condlist;

goodtrials = [];
for i = 1:numel(cl)
    goodtrials  = [goodtrials indtrial(D, cl{i}, 'GOOD')];
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
if strncmpi(D.transformtype,'TF',2) % TF and TFphase
    Dnew = clone(D, [S.prefix fname(D)], [D.nchannels D.nfrequencies D.nsamples length(goodtrials)]);
else
    Dnew = clone(D, [S.prefix fname(D)], [D.nchannels D.nsamples length(goodtrials)]);
end

%-Copy data
%--------------------------------------------------------------------------
spm_progress_bar('Init', length(goodtrials), 'Trials copied');
if length(goodtrials) > 100, Ibar = floor(linspace(1, length(goodtrials), 100));
else Ibar = [1:length(goodtrials)]; end

for i = 1:length(goodtrials)
    
    if strncmpi(D.transformtype,'TF',2) % TF and TFphase
        Dnew(:, :, :, i) =  D(:, :, :, goodtrials(i));
    else
        Dnew(:, :, i) =  D(:, :, goodtrials(i));
    end    
    
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end  %

spm_progress_bar('Clear');

%-Copy trial-specific data.
%--------------------------------------------------------------------------
Dnew = conditions(Dnew, ':', conditions(D, goodtrials));
Dnew = repl(Dnew, ':', repl(D, goodtrials));
Dnew = events(Dnew, ':', events(D, goodtrials));
Dnew = trialonset(Dnew, ':', trialonset(D, goodtrials));
Dnew = trialtag(Dnew, ':', trialtag(D, goodtrials));

%-Save the new M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history(mfilename, S);
save(Dnew);

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','Remove bad trials: done'); spm('Pointer','Arrow');
