function [b] = spm_fir(N,Wn,M)
% designs an N'th order FIR digital filter
% FORMAT [b] = spm_fir(N,Wn,M)
%
% N 		: order
% Wn		: n vector of transition Hz (0 < Wn < 1)
% M      	: n + 1 vector band pass specification for each Hz interval
%		  e.g. [1 0] - low  pass
%		       [0 1] - high pass
%		       [0 1 0] - band pass
%
% b 		: filter
%__________________________________________________________________________
$
% FIR filter design using least-squares error minimization. Based on firls.m
% in the signal processing toolbox. Wn = 1 corresponds to the Nyquist
% frequency or half the sampling frequency.
%--------------------------------------------------------------------------
% %W% Karl Friston %E%

% specify bands
%--------------------------------------------------------------------------
M     = [M(:)'; M(:)'];
F     = [0,Wn(1:length(Wn)); Wn(1:length(Wn)),1];
N     = N + 1;
F     = F(:)/2;
M     = M(:);
L     = (N - 1)/2;
Nodd  = rem(N,2);
if ~Nodd
        m = (0:L) + .5;
else
        m = (0:L);
end
k     = m';
if Nodd
	k  = k(2:length(k));
	b0 = 0;
end;
b     = zeros(size(k));

% Compute filter
%--------------------------------------------------------------------------
for s =1:2:length(F)

        m  =(M(s + 1) - M(s))/(F(s + 1) - F(s));	%  slope
        b1 = M(s) - m*F(s);				%  intercept
        if Nodd
		b0 = b0 + (b1*(F(s+1)-F(s)) + m/2*(F(s+1)*F(s+1)-F(s)*F(s)));
        end
        b  = b + (m/(4*pi*pi)*(cos(2*pi*k*F(s+1))-cos(2*pi*k*F(s)))./(k.*k));
        b  = b + (F(s+1)*(m*F(s+1)+b1)*sinc(2*k*F(s+1)) ...
                   - F(s)*(m*F(s)+b1)*sinc(2*k*F(s)));
end;
if Nodd
	b=[b0; b];
end;
a     = 4*b;
if Nodd
	a(1) = a(1)/2;
end
if Nodd
	h    = [a(L+1:-1:2)/2; a(1); a(2:L+1)/2].';
else
        h    = .5*[flipud(a); a].';
end;

% Apply Hanning window (or not)
%--------------------------------------------------------------------------
W      = 0.54 - 0.46*cos(2*pi*(0:N-1)'/(N-1));
b      = h.*W(:)';
b      = h;

