function [P] = spm_Pkn(s0,h0, s,u,S)
%
% probability based on the expected number of cluster
% with spatial extent s0 and maximum height h0.
%
% FORMAT [P] = Pkn(s0,h0, s,u,S);
%
% s0  - number of voxels in the region
% h0  - maximum height in the region
% s   - smoothness: length(s) = D {dimension}
% u   - height threshold
% S   - Lebesgue measure of S {volume in voxels}
%
% P   - Pkn(s0,h0, s,u,S)
%___________________________________________________________________________
%
% spm_Pkn returns the probability of one or more regions defined
% with the threshold u with maximum height h0 (above u)
% and spatial extent s0 in volume S of a D - dimensional Gaussian
% process of non-isotropic smoothness s, thresholded at u.
%
% spm_Pkn calls : spm_Pn_min, spm_Pz_min, spm_P
%
% Ref: Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
% Ref: Friston et al (1994) Assessing the significance of focal
% activations with their spatial extent. HBM 1:210-220
% Ref: Poline J-B, Worsley K., Evans A., Friston K. NeuroImage,
% submited.
%_______________________________________________________________________
% @(#)spm_Pkn.m	1.3 Jean-Baptiste Poline 96/08/23


%------ Check h0 and s0 are positive numbers 	--------%
h0 = h0(h0>=0);
s0 = s0(s0>0);

s0 = s0(:)';
h0 = h0(:)';

% Check that the vectors sizes are compatible : %
if(length(s0) ~= length(h0) & length(h0) )
	disp(['	length(s0) ~= length(h0) in Pkn ']);
 	return;
end;

%----------------------------------------------------------%
D	= max(size(s));		% Dimension of the process
L	= diag(1./(2*(s.^2)));	% s -> Lambda
nu	= 4*u^2/D;		% Degree of freedom 

Emu	= S * sqrt(det(L))*(u^(D-1))*exp(-(u.^2)/2) / (2*pi)^(D/2 + .5);
EN 	= S * spm_Ncdf(-u);
Enu	= EN/Emu;
a	= Enu*u^(D/2)/gamma(D/2 + 1);
beta    = (gamma(D/2 + 1)/Enu)^(2/D);

for i=1:length(s0)
     Pn(i) = spm_P(1,s,u,s0(i),S);	% Probability based on spatial extent
     Pz(i) = spm_P(1,s,u+h0(i),0,S); 	% Probability based on peak height
end

%---- Make sure the results are between 0 and 1 :-------%
Pn(find(Pn>1)) = ones(size(find(Pn>1)));
Pn(find(Pn<0)) = zeros(size(find(Pn<0)));
Pz(find(Pz>1)) = ones(size(find(Pz>1)));
Pz(find(Pz<0)) = zeros(size(find(Pz<0)));

for i=1:length(Pn)

	%---------------------------------------------------------------%
	% get the minimum probability and find the parameter such   	%
	% that these probabilities equal.				%
	OPT = zeros(1,14); OPT(2) = .1; OPT(14) = 100;
	if Pn(i) > Pz(i)
		n_opt = fmin('spm_min_Pn', s0(i), 12*s0(i),...
				OPT, Pz(i), s, u, S);
		s0(i) = n_opt;
	else
		OPT(2) = .01;
		z_opt = fmin('spm_min_Pz', h0(i)+u, h0(i)+u+5,...
				OPT, Pn(i), s, S);
		h0(i) = z_opt - u;
	end;

	%---------------------------------------------------------------%
	% Computation of Pkn : integration over hh 		   	%
	hh = [h0(i):0.02:5 + 1.5*h0(i)];
	b	= a*hh.^(D/2)*(nu)/s0(i);
	ps	= spm_Xcdf(b,nu);
	theo_s0_h0	= mean(ps.*u.*exp(-u.*hh)).*(max(hh)-min(hh));
	
	theo_1_h0 = exp(-beta.*s0(i).^(2/D));
	theo_s0_0 = exp(-u.*h0(i));
	
	%---------------------------------------------------------------%
	% Probability equals the marginals - the previous integration	%
	P(i) =  theo_1_h0 + theo_s0_0 - theo_s0_h0;
	P(i) = 1 - exp(-Emu*P(i));

end

