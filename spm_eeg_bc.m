function d = spm_eeg_bc(D, d)
% 'baseline correction' for D: subtract average baseline energy of the 
% samples (start:stop) per epoch.
% FORMAT d = sjk_eeg_bc(D, d)
%
%_______________________________________________________________________
% @(#)spm_eeg_bc.m	1.1 Stefan Kiebel 04/06/28

for i = 1 : length(D.tf.channels)
	for j = 1 : D.Nfrequencies
		tmp = mean(d(i, j, D.tf.Sbaseline(1):D.tf.Sbaseline(2)), 3);
		d(i, j, :) = d(i, j, :) - tmp;
	end
end
