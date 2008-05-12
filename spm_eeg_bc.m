function D = spm_eeg_bc(D, time)
% 'baseline correction' for data in D: subtract average baseline from all
% M/EEG and EOG channels
% FORMAT d = spm_eeg_bc(D, time)
%
% D:    meeg object
% time: 2-element vector with start and end of baseline period [ms]
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_bc.m 1598 2008-05-12 12:06:54Z stefan $

t(1) = D.indsample(time(1));
t(2) = D.indsample(time(2));

indchannels = [D.meegchannels D.eogchannels];

spm_progress_bar('Init', D.ntrials, 'trials baseline-corrected'); drawnow;
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = [1:D.ntrials]; end

switch(transformtype(D))
    case 'TF'
        for k = 1: D.ntrials
            tmp = mean(D(:, :, t(1):t(2), k), 3);
            D(:, :, :, k) = D(:, :, :, k) - repmat(tmp, [1, 1, D.nsamples, 1]);

            if ismember(k, Ibar)
                spm_progress_bar('Set', k);
                drawnow;
            end
        end

    case 'time'
        for k = 1: D.ntrials
            tmp = mean(D(indchannels, t(1):t(2), k), 2);
            D(indchannels, :, k) = D(indchannels, :, k) - repmat(tmp, 1, D.nsamples);
            
            if ismember(k, Ibar)
                spm_progress_bar('Set', k);
                drawnow;
            end

        end


    otherwise
        error('Unknown transform type');

end

spm_progress_bar('Clear');

