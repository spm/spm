function [hrf,p] = spm_hrf(RT,p);
% returns a hemodynamic response function
% FORMAT [hrf,p] = spm_hrf(RT,[p]);
% RT   - scan repeat time
% p    - parameters of the response function (two gamma functions)
%
%							defaults
%							(seconds)
%	p(1) - delay of response (relative to onset)	   6
%	p(2) - delay of undershoot (relative to onset)    16
%	p(3) - dispersion of response			   1
%	p(4) - dispersion of undershoot			   1
%	p(5) - ratio of response to undershoot		   6
%	p(6) - onset (seconds)				   0
%
% hrf  - hemodynamic response function
% p    - parameters of the response function
%_______________________________________________________________________
% %W% Karl Friston %E%

% default parameters
%-----------------------------------------------------------------------
if nargin < 2
	p = [6 16 1 1 6 0];
end

% modelled hemodynamic response function - {mixture of Gammas}
%-----------------------------------------------------------------------
dt    = RT/8;
u     = [0:(32/dt)] - p(6)/dt;
hrf   = spm_Gpdf(u,[p(1) dt]/p(3)) - spm_Gpdf(u,[p(2) dt]/p(4))/p(5);
hrf   = spm_conv(hrf,8);
hrf   = hrf([1:32/RT]*8);
hrf   = hrf/sum(hrf);
