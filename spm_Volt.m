function [G,I,W] = spm_Volt(s,RT)
% returns second order [design] matrix of explanatory variables
% FORMAT [G I W] = spm_Volt(s,RT);
% G  -  Design Matrix
% I  -  indices of columns pertaining to 1st and 2nd order terms
% W  -  basis functions employed
% s  -  multivariate time-series (one column per component)
% RT -  interscan interval
%___________________________________________________________________________
% %W% Karl Friston %E%


% parameters
%---------------------------------------------------------------------------
T     = size(s,1);				% number of observations
n     = size(s,2);				% n-variate system

% basis fuctions
%---------------------------------------------------------------------------
W     = spm_Volt_W([0:RT:32]);
w     = size(W,2);

% Construct G
%---------------------------------------------------------------------------
m     = n*w + (n*(n + 1)/2)*w*w;
G     = zeros(T,m);
I     = zeros(2,m);

% 1st order terms
%---------------------------------------------------------------------------
g     = 0;
d     = 1:T;
for i = 1:n
	for j = 1:w
		g      = g + 1;
		x      = s(:,i);
		x      = conv(x,W(:,j));
		x      = x(d);
		G(:,g) = x;
		I(1,g) = i;
	end
end


% 2nd order terms
%---------------------------------------------------------------------------
for i = 1:n
	for j = i:n
		for p = 1:w
			for q = 1:w
				g      = g + 1;
				x      = s(d,i);
				y      = s(d,j);
				x      = conv(x,W(:,p));
				y      = conv(y,W(:,q));
				x      = x(d);
				y      = y(d);
				G(:,g) = x.*y;
				I(1,g) = i;
				I(2,g) = j;
			end
		end
	end
end
