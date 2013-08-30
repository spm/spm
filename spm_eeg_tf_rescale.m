function D = spm_eeg_tf_rescale(S)
% Rescale (avg) spectrogram with nonlinear and/or difference operator
% FORMAT [D] = spm_eeg_tf_rescale(S)
%
% S                    - input structure (optional)
% fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.method           - 'LogR', 'Diff', 'Rel', 'Log', 'Sqrt', 'None'
%   S.timewin          - 2-element vector: start and stop of baseline (ms)
%                        (need to specify this for LogR and Diff)
%   S.Db               - MEEG object or filename of M/EEG mat-file to use
%                        for the baseline (if different from the input dataset).
%   prefix             - prefix for the output file (default - 'r')
%
% Output:
%   D                  - MEEG object with rescaled power data (also
%                        written to disk with prefix r)
%
% For 'Log' and 'Sqrt', these functions are applied to spectrogram
% For 'LogR', 'Rel' and 'Diff' this function computes power in the baseline
% p_b and outputs (i) p-p_b for 'Diff' (ii) 100*(p-p_b)/p_b for 'Rel'
%                 (iii) log (p/p_b) for 'LogR'
%__________________________________________________________________________
% Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_eeg_tf_rescale.m 5626 2013-08-30 14:03:16Z vladimir $

SVNrev = '$Rev: 5626 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Time-frequency rescale'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix   = 'r';           end
if ~isfield(S, 'timewin'),      S.timewin  = [-Inf 0];      end

D = spm_eeg_load(S.D);

Dnew  = clone(D, [S.prefix fname(D)]);

needbaseline = ismember(lower(S.method), {'logr','diff', 'rel', 'zscore'});

if needbaseline    
    if isfield(S, 'Db') && ~isempty(S.Db)
        Db = spm_eeg_load(S.Db);
    else
        Db = D;
    end
    
    if any(abs(D.frequencies - Db.frequencies) > 0.1) || ~isequal(Db.chanlabels, D.chanlabels) ||...
            (Db.ntrials > 1 && (Db.ntrials ~= D.ntrials))
        error('The input dataset and the baseline dataset should have the same frequencies, channels and trial numbers');
    end
    
    timeind = Db.indsample(1e-3*(min(S.timewin))):Db.indsample(1e-3*(max(S.timewin)));
    if isempty(timeind) || any(isnan(timeind))
        error('Selected baseline time window is invalid.');
    end
end


spm_progress_bar('Init', D.ntrials, 'Trials processed');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end

for i = 1:D.ntrials
    
    if needbaseline
        x = spm_squeeze(D(:,:,:,i), 4);
        if Db.ntrials > 1
            xbase = spm_squeeze(Db(:,:,:,i), 4);
        else
            xbase = spm_squeeze(Db(:,:,:,1), 4);
        end
    end
    
    switch lower(S.method)
        
        case 'logr'
            xbase         = mean(log10(xbase(:,:,timeind)),3);
            Dnew(:,:,:,i) = 10*(log10(x) - repmat(xbase,[1 1 Dnew.nsamples 1]));
            
        case 'diff'
            xbase         = mean(xbase(:, :, timeind),3);
            Dnew(:,:,:,i) =  (x - repmat(xbase,[1 1 Dnew.nsamples 1]));
            
        case 'zscore'
            stdev            =  std(xbase(:, :, timeind), [], 3);
            xbase            =  mean(xbase(:,:, timeind),3);            
            Dnew(:, :, :, i) = (x - repmat(xbase,[1 1 Dnew.nsamples 1]))./repmat(stdev,[1 1 Dnew.nsamples 1]);
            
        case 'rel'
            xbase            = mean(xbase(:, :, timeind), 3);
            Dnew(:, :, :, i) = 100*((x./repmat(xbase,[1 1 Dnew.nsamples 1]) - 1));
            
        case 'log'
            Dnew(:, :, :, i) = log(D(:, :, :, i));
            
        case 'logeps'
            Dnew(:, :, : ,i) = log(D(:, :, :, i)+(1/64)*mean(spm_vec(D(:, :, :, i))));
            
        case 'sqrt'
            Dnew(:, :, :, i) = sqrt(D(:, :, :, i));
            
        case 'none'
            Dnew(:, :, :, i) = D(:, :, :, i);
            
        otherwise
            error('Unknown rescaling method');
            
    end
        
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

switch lower(S.method)
    
    case 'logr'
        Dnew         = units(Dnew, ':', 'dB');
        
    case 'rel'
        Dnew         = units(Dnew, ':', '%');
end

D = Dnew;

% Save
D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Time Frequency Rescale: done');
spm('Pointer','Arrow');
