function [K] = spm_make_filter(k,RT,filterHF,filterLF)
% returns a convolution matrice for bandpass filtering
% FORMAT [K] = spm_make_filter(k,RT,filterHF,filterLF)
%
% k 		: number of scans
% RT		: Repetition time in seconds
% filterHF	: High Frequency specification (Low  pass)
% filterLF	: Low  Frequency specification (High pass)
% norm 		: normalization
%__________________________________________________________________________
% Routine to construct band pass filter convolution matrices;
% First, a high pass filter is designed (specified in filterLF),
% then a low pass (filterHF). These are combined in the temporal domain.
% See spm_fmri_spm_ui for filter structure.
% norm = 'norm' for normalisation such that the signal variance
% is unchanged by the filter.
%------------------------------------------------------------------
% %W% Jean-Baptiste Poline %E%


% Low frequencies   '64 second cut-off'|'specify'|'no'
%   '64 second cut-off'
%   'specify'
%   '64 second cut-off -FIR'
%   'specify -FIR'
%   'no'
%------------------------------------------------------------------
switch filterLF.Choice

	case 'no'
		KLF = speye(k);
	
	case {'64 second cut-off -FIR','specify -FIR'}
		freqLF = 2*RT/filterLF.Param;
		n      = 32;	
		FIL    = spm_fir(n,freqLF,[0 1]);
		FIL    = [FIL([1:n/2] + n/2) zeros(1,k - n/2)];
		KLF    = sparse(toeplitz(FIL));

	case {'64 second cut-off','specify'}

		n      = fix(2*(k*RT)/filterLF.Param + 1);
		X      = spm_dctmtx(k,n);
		X      = X(:,[2:n]);
		KLF    = eye(k) - X*X';

	otherwise

		warning('High pass Filter option unknown');

end
	
% High frequencies 
%   'smooth with hrf'
%   'smooth with hrf derivative'
%   'smooth with Gaussian kernel'
%   'none'
%   'specify'
%------------------------------------------------------------------
switch filterHF.Choice

	case 'none'

		KHF = eye(k);

	case 'specify'
		freqHF = 2*RT/filterHF.Param;
		n      = 32;	
		FIL    = spm_fir(n,freqLF,[1 0]);
		FIL    = [FIL([1:n/2] + n/2) zeros(1,k - n/2)];
		KHF    = sparse(toeplitz(FIL));


	case 'smooth with hrf'  	
		h      = spm_hrf(RT);
		h      = [h; zeros(size(h))];
		n      = length(h);
		R      = conv(h,flipud(h));
		R      = R([1:n] + n - 1);
		V      = toeplitz(R);
		K      = sqrtm(V);
		FIL    = K(n/2,[1:n/2] + n/2 -1);
		FIL    = [FIL zeros(1,k - n/2)];
		KHF    = sparse(toeplitz(FIL));

	case 'smooth with hrf derivative'
		h      = gradient(spm_hrf(RT));
		h      = [h; zeros(size(h))];
		n      = length(h);
		R      = conv(h,flipud(h));
		R      = R([1:n] + n - 1);
		V      = toeplitz(R);
		K      = sqrtm(V);
		FIL    = K(n/2,[1:n/2] + n/2 -1);
		FIL    = [FIL zeros(1,k - n/2)];
		KHF    = sparse(toeplitz(FIL));
	
	case 'smooth with Gaussian kernel'
		FIL    = exp(-[0:(k - 1)].^2/(2*(filterHF.Param/RT)^2));
		KHF    = toeplitz(FIL);

	otherwise
		warning('Low pass Filter option unknown');

end


% Combine and normalize
%------------------------------------------------------------------
KHF  = spdiags(1./sum(KHF')',0,k,k)*KHF;
K    = KHF*KLF;
