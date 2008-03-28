function D = spm_eeg_bc(D, time)
% 'baseline correction' for data in D: subtract average baseline from data
% FORMAT d = spm_eeg_bc(D, time)
%
% D:    meeg object
% time: 2-element vector with start and end of baseline period [ms]
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_bc.m 1278 2008-03-28 18:38:11Z stefan $

t(1) = D.indsample(time(1));
t(2) = D.indsample(time(2));

switch(transformtype(D))
    case 'TF'
        for i = 1 : D.nchannels
            for j = 1 : D.nfrequencies
                for k = 1: D.ntrials
                    tmp = mean(D(i, j, t(1):t(2), k), 3);
                    D(i, j, :, k) = D(i, j, :, k) - tmp;
                end
            end
        end
    case 'time'
        for i = 1 : D.nchannels
            for k = 1: D.ntrials
                tmp = mean(D(i, t(1):t(2), k), 2);
                D(i, :, k) = D(i, :, k) - tmp;
            end
        end
    otherwise
        error('Unknown transform type');

end
