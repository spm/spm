function d = spm_eeg_bc(D, d)
% 'baseline correction' for D: subtract average baseline energy of the 
% samples (start:stop) per epoch.
% FORMAT d = spm_eeg_bc(D, d)
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_bc.m 1143 2008-02-07 19:33:33Z spm $


for i = 1 : length(D.tf.channels)
    for j = 1 : D.Nfrequencies
        tmp = mean(d(i, j, D.tf.Sbaseline(1):D.tf.Sbaseline(2)), 3);
        d(i, j, :) = d(i, j, :) - tmp;
    end
end
