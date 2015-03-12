function D = spm_eeg_split_conditions(S)
% Splits a file into different conditions in order to facilitate TF
% processing. The idea is to create several smaller files, run TF, then
% aveage within the condition files using spm_eeg_average_tf and lastly,
% merge again.
% FORMAT D = spm_eeg_split_conditions(S)
%
% S        - optional input struct
% (optional) fields of S:
% D        - MEEG object or filename of M/EEG mat-file with epoched data
%
% Output:
% D        - MEEG object (also written on disk)
%
% The function also physically removes bad trials.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Dominik R Bach
% based on spm_eeg_remove_bad_trials

SVNrev = '$Rev: 6374 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Remove bad trials'); spm('Pointer','Watch');

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

%-Check that there is any good data available
%--------------------------------------------------------------------------
if ntrials(D)==0 || all(badtrials(D, ':'))
    warning('No good trials were found. Nothing to do.');
    return;
end

%-select trials
%--------------------------------------------------------------------------
cl   = D.condlist;

goodtrials = [];
for i = 1:numel(cl)
    goodtrials{i}  = indtrial(D, cl{i}, 'GOOD');
end

%-Generate new files & copy data
%--------------------------------------------------------------------------

for i = 1:length(goodtrials)
    spm_progress_bar('Init', length(goodtrials{i}), 'Trials copied');
    if length(goodtrials{i}) > 100, Ibar = floor(linspace(1, length(goodtrials{i}), 100));
    else Ibar = [1:length(goodtrials{i})]; end
    
    fn = fname(D);
    fn = ['r', fn(1:(end-4)), sprintf('_%02.0f.dat', i)];
    if strncmpi(D.transformtype,'TF',2) % TF and TFphase
        Dnew = clone(D, fn, [D.nchannels D.nfrequencies D.nsamples numel(goodtrials{i})]);
    else
        Dnew = clone(D, fn, [D.nchannels D.nsamples numel(goodtrials{i})]);
    end
    
    for trl = 1:numel(goodtrials{i})
        if strncmpi(D.transformtype,'TF',2) % TF and TFphase
            Dnew(:, :, :, trl) =  D(:, :, :, goodtrials{i}(trl));
        else
            Dnew(:, :, trl) =  D(:, :, goodtrials{i}(trl));
        end
    end;
    
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
   
    %-Copy trial-specific data.
    %--------------------------------------------------------------------------
    Dnew = conditions(Dnew, ':', conditions(D, goodtrials{i}));
    Dnew = repl(Dnew, ':', repl(D, goodtrials{i}));
    Dnew = events(Dnew, ':', events(D, goodtrials{i}));
    Dnew = trialonset(Dnew, ':', trialonset(D, goodtrials{i}));

    %-Save the new M/EEG dataset
    %--------------------------------------------------------------------------
    Dnew = Dnew.history(mfilename, S);
    save(Dnew);

    DD{i} = Dnew;
end  %

spm_progress_bar('Clear');

D = DD;

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','Split conditions: done'); spm('Pointer','Arrow');
