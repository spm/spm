function [y1,y] = spm_ppi(x,phi,rt,p)
% Does the thing mentioned in Methods meeting 20/03/00
% FORMAT [y1,y] = spm_ppi_thingy(x,phi,rt[,p])
%             y1  - processed data
%             y   - estimated deconvolved x
%             x   - raw data
%             phi - the effect of interest
%             rt  - repeat time
%             p   - other assorted parameters - see spm_hrf.m.
%
% Assume that x is H*y, where H is a Toeplitz convolution matrix
% that convolves with the hæmodynamic response, and y is the
% underlying neuronal activity.  The problem tries to deconvolve
% the hrf from x.  Restricted maximum likelihood is used to estimate
% the optimal weighting (w) for y = (w(1)*H'*H + w(2)*I)\H'*x.
% After deconvolution, the output is
%            y1 = H*diag(phi)*y 
%_______________________________________________________________________
% %W% John Ashburner %E%

if nargin<4,
	hrf = spm_hrf(rt);
else,
	hrf = spm_hrf(rt,p);
end;
x         = x(:);
n         = size(x,1);
l         = min(length(hrf),n);
hrf1      = zeros(n,1); hrf2 = zeros(n,1);
hrf1(1:l) = hrf(1:l); hrf2(1) = hrf(1);
H         = toeplitz(hrf1,hrf2);
y         = reml(H,x);
y1        = H*(phi(:).*y);
return;

function [bet,w]=reml(X,y)
% Restricted Maximum Likelihood estimation.
n  = size(X,1);
w  = [1 1];
Xy = X'*y;
XX = X'*X;
obet = Inf*ones(n,1);
for i=1:100
	T    = inv(w(1)*XX + w(2)*eye(n));
	bet  = T*(w(1)*Xy);
	p1   = trace(T*XX*w(1));
	p2   = n-p1;
	w(1) = (n-p1)/sum((y-X*bet).^2);
	w(2) = (n-p2)/sum(bet.^2);
	if sum(((bet-obet)./(bet+obet+1e-20)).^2) < 0.000001*n,
		break;
	end;
	obet = bet;
	fprintf('.');
end;
fprintf('\n');
return;
