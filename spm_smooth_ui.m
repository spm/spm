function spm_smooth_ui
% user interface routine for spm_conv_vol
% FORMAT spm_smooth_ui
%____________________________________________________________________________
%
% spm_smooth_ui sets up a list of images to smooth or convolve with a 
% Gaussian kernel of specified width. Convolved images are prefixed with 's'
% (headers are created automatically)
%
%__________________________________________________________________________
% %W% %E%

% get filenames and kernel width
%----------------------------------------------------------------------------
set(2,'Name','Smoothing')

s     = spm_input('smoothing {FWHM in mm}',1);
P     = spm_get(Inf,'.img','select scans');
n     = size(P,1);
% implement the convolution
%---------------------------------------------------------------------------
set(2,'Name','executing','Pointer','watch')
for i = 1:n
	Q = P(i,:);
	Q = Q(Q ~= ' ');
	d = max([find(Q == '/') 0]);
	U = [Q(1:d) 's' Q((d + 1):length(Q))];
	if ~strcmp(U([1:4] + length(U) - 4),'.img'); U = [U '.img']; end
	spm_smooth(Q,U,s);
end
set(2,'Name','','Pointer','arrow'); figure(2); clf
