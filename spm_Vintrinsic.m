function  [V] = spm_Vintrinsic(n,RT,form,Vpar)
%--------------------------------------------------------------------------
% Compute an assumed intrinsic autocorelation function
% FORMAT [V] = spm_Vintrinsic(nscan,form,Vparam)
% nscan		- sessions number of scans
% RT		- Repetition time in seconds
% form		- '1/f' 
% Vpar		- Parameters pertaining to the form
%
% V		- (intrinsic) autoccorrelation matrix
%__________________________________________________________________________

if nargin < 4 error('not enough argument in spm_Vintrinsic'); end;



if rem(n,2) == 0
 	f=1/(n*RT):1/(n*RT):1/(2*RT)+1/(n*RT);
else 
	f=1/((n-1)*RT):1/((n-1)*RT):1/(2*RT)+1/((n-1)*RT);
end



switch form
  case('1/f')
   SD	= (Vpar(1)./f + 1).^2;	% Spectral density form
   SDS= [fliplr(SD(2:length(SD))) SD];
  otherwise, error('unrecognised type')
end % switch form

R	= real(ifft(SDS));
R	= R/max(R); 
V  	= toeplitz(R(1:n));

