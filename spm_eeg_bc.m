function d = spm_eeg_bc(D, d)
% 'baseline correction' for D: subtract average baseline energy of the 
% samples (start:stop) per epoch.
% FORMAT d = spm_eeg_bc(D, d)
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_bc.m 133 2005-05-09 17:29:37Z guillaume $


for i = 1 : length(D.tf.channels)
	for j = 1 : D.Nfrequencies
		tmp = mean(d(i, j, D.tf.Sbaseline(1):D.tf.Sbaseline(2)), 3);
		d(i, j, :) = d(i, j, :) - tmp;
	end
end
