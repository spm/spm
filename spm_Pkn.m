function [P] = Pnz(s0,h0, s,u,V)
% A JB special
% FORMAT
%_______________________________________________________________________
% %W% Jean-Baptiste Poline %E%

% keyboard;
h0 = h0(h0>=0);
s0 = s0(s0>0);

s0 = s0(:)';
h0 = h0(:)';

if(length(s0) ~= length(h0) & length(h0) )
	disp(['	length(s0) ~= length(h0) in Pnk ']);
 	return;
end;

D	= max(size(s));
L	= diag(1./(2*(s.^2)));
nu	= 4*u^2/D;

Emu	= V * sqrt(det(L))*(u^(D-1))*exp(-(u.^2)/2) / (2*pi)^(D/2 + .5);
EN 	= V * spm_Ncdf(-u);
Enu	= EN/Emu;
a	= Enu*u^(D/2)/gamma(D/2 + 1);
beta    = (gamma(D/2 + 1)/Enu)^(2/D);


Pn = spm_Pn(s0,s,u,V);
Pz = spm_Pz(s,u+h0,V);

Pn(find(Pn>1)) = ones(size(find(Pn>1)));
Pn(find(Pn<0)) = zeros(size(find(Pn<0)));
Pz(find(Pz>1)) = ones(size(find(Pz>1)));
Pz(find(Pz<0)) = zeros(size(find(Pz<0)));

for i=1:length(Pn)

	OPT = zeros(1,14); OPT(2) = .1; OPT(14) = 100;
	if Pn(i) > Pz(i)
		n_opt = fmin('spm_min_Pn', s0(i), 12*s0(i),...
				OPT, Pz(i), s, u, V);
		s0(i) = n_opt;
	else
		OPT(2) = .01;
		z_opt = fmin('spm_min_Pz', h0(i)+u, h0(i)+u+5,...
				OPT, Pn(i), s, V);
		h0(i) = z_opt - u;
	end;

	hh = [h0(i):0.02:5 + 1.5*h0(i)];
	b	= a*hh.^(D/2)*(nu)/s0(i);
%	ps	= spm_Chi2cdf(nu,b);
	ps	= 1 - spm_Xcdf(b,nu);
	theo_s0_h0	= mean(ps.*u.*exp(-u.*hh)).*(max(hh)-min(hh));
	
	theo_1_h0 = exp(-beta.*s0(i).^(2/D));
	theo_s0_0 = exp(-u.*h0(i));

	P(i) =  theo_1_h0 + theo_s0_0 - theo_s0_h0;
	P(i) = 1 - exp(-Emu*P(i));

end

