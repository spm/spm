function [O] = spm_conv(X,sx,sy)
% Gaussian convolution of a matrix variable
% FORMAT [O] = spm_conv(X,sx,[sy]);
% X    - matrix
% sx   - kernel width (FWHM) in pixels
% sy   - optional non-isomorphic smoothing
%___________________________________________________________________________
%
% spm_conv is a one or two dimensional convolution of a matrix variable
% in working memory.  It capitalizes on the sparisity structure of the
% problem and the separablity of multidimensional convolution with a Gaussian
% kernel by using one-dimensional convolutions and kernels that are
% restricted to non near-zero values
%
%__________________________________________________________________________
% %W% %E%


%---------------------------------------------------------------------------
if nargin == 2
	s = sx;
	if sx == 0; O = X; return; end;
	[lx ly] = size(X);
        s       = s/sqrt(8*log(2));
	E       = round(s*4);
        x       = exp(-[-E:E].^2/(2*s^2));
        x       = x/sum(x);
	if min([lx ly]) == 1
		X = conv(X(:),x);
		O = X((E + 1):max([lx ly]) + E);
		return;
	end
	O    = zeros((lx + 2*E),(ly + 2*E));
	O((E + 1):(E + lx),(E + 1):(E + ly)) = X;
	X    = conv(O(:),x);
        X    = X(E+1:E+(lx+2*E)*(ly+2*E));
        O(:) = X; X = O';
	X    = conv(X(:),x);
        X    = X(E+1:E+(lx+2*E)*(ly+2*E));O=O';
        O(:) = X; O = O';
	O    = O(E+1:E+lx,E+1:E+ly);
else
	if sx == 0 & sy == 0; O=X; return; end;
	[lx ly] = size(X);
	if sx ~= 0;
		sx=sx/sqrt(8*log(2));
		Ex = round(s*4);

		x=exp(-[-Ex:Ex].^2/(2*sx^2));
		xx=x/sum(x);
	end
	if sy ~= 0;
		sy = sy/sqrt(8*log(2));
		Ex = round(s*4);
		x  = exp(-[-Ey:Ey].^2/(2*sy^2));
		xy = x/sum(x);
	end
	E    = max([Ex Ey]);
	O    = zeros(lx+2*E,ly+2*E);
    	O(E+1:E+lx,E+1:E+ly)=X;
	if sx ~= 0;
		X = conv(O(:),xx);
		X = X(Ex+1:Ex+(lx+2*E)*(ly+2*E)); O(:) = X;
	end
	if sy ~= 0;
		X = O'; X = conv(X(:),xy);
		X = X(Ey+1:Ey+(lx+2*E)*(ly+2*E));
		O = O'; O(:) = X; O = O';
	end
	O    = O(E+1:E+lx,E+1:E+ly);
end














% weighting in z for 16 mm FWHM
% 0.5000    0.8409   1.0000    0.8409    0.5000

% weighting in z for 12 mm FWHM
% 0.2916    0.7349   1.0000    0.7349    0.2916

% weighting in z for 24 mm FWHM
% 0.0577    0.0727   0.0785    0.0727    0.0577 
