function [G] = spm_global(V)
% returns the global mean for a memory mapped volume image
% FORMAT [G] = spm_global(V)
% V   - memory mapped volume
% G   - mean global activity
%____________________________________________________________________________
%
% spm_global returns the mean counts integrated over all the  
% slices from the volume
%
% The mean is estimated after discounting voxels outside the object
% using a criteria of greater than > (global mean)/8
%
% see also spm_box.m
%
%__________________________________________________________________________
% %W% %E%

%----------------------------------------------------------------------------
X	= zeros(V(1)*V(2),V(3));
for   i = 1:V(3);
	d    = spm_matrix([0 0 i]);
	d    = spm_slice_vol(V,d,[V(1) V(2)],0);
	X(:,i) = d(:);
end

G	= mean(X(X(:) > mean(X(:))/8));
