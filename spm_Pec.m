function [CPF,EC,R] = spm_Pec(STATISTIC,SPACE,t,L,W,df)
% Returns the corrected P value using unifed EC theory
% FORMAT [CPF,EC,R] = spm_Pec(STATISTIC,SPACE,U,L,W,df)
%
% STATISTIC - Statisical feild
%		'Z' - Gaussian feild
%		'F' - Gaussian feild
% SPACE     - Search space
%		'S' - Sphere
%		'B' - Box
%		'V' - Discrete voxels
% U         - value[s] of statistic
% L         - space definition {in voxels}
%		L = radius        {Sphere}
%		L = [a b c]       {Box}
%		L = XYZ pointlist {Discrete voxels}
% W         - smoothness of the component fields
% df        - [df(interest) df(residuals)] {F}
%
% CPF       - corrected P value
% EC        - Euler characteristic densities
% R         - Resel counts
%
%___________________________________________________________________________
%
% Reference : Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
%_______________________________________________________________________
% %E% Karl Friston %W%


% Dimensionality
%---------------------------------------------------------------------------
if size(W,1) > 1, W = mean(W); end
D     = size(W,2);
t     = t(:)';
N     = length(t);
R     = zeros(4,N);
EC    = zeros(4,N);
w     = W*sqrt(8*log(2));


% Resel Counts (R)
%===========================================================================
if      SPACE == 'S'

	% Sphere
	%-------------------------------------------------------------------
	s  = L(:)./w(:);
	s  = prod(s).^(1/D);
	R  = [1 4*s 2*pi*s^2 (4/3)*pi*s^3];

elseif  SPACE == 'B'

	% Box
	%-------------------------------------------------------------------
	s  = L(:)./w(:);
	R  = [1 sum(s) (s(1)*s(2) + s(2)*s(3) + s(1)*s(3)) prod(s)];

elseif  SPACE == 'V'

	% Voxels
	%-------------------------------------------------------------------
	R  = spm_Pec_resels(L,W);

end


% EC densities (EC}
%===========================================================================
if      STATISTIC == 'Z'

	% Gaussian Field
	%-------------------------------------------------------------------
	EC(1,:) = 1 - spm_Ncdf(t);
	EC(2,:) = (4*log(2))^(1/2)/(2*pi)*exp(-t.^2/2);
	EC(3,:) = (4*log(2))/((2*pi)^(3/2))*exp(-t.^2/2).*t;
	EC(4,:) = (4*log(2))^(3/2)/((2*pi)^2)*exp(-t.^2/2).*(t.^2 - 1);

elseif  STATISTIC == 'F'

	% F Field
	%-------------------------------------------------------------------
	k       = df(1);
	v       = df(2);
	v       = min([v (340 - k)]);
	a       = (4*log(2))/(2*pi);
	b       = gamma(v/2)*gamma(k/2);
	if k + v <= D, error('df too small'); end

	EC(1,:) = 1 - spm_Fcdf(t,df);
	EC(2,:) = a^(1/2)*gamma((v+k-1)/2)*2^(1/2)/b*(k*t/v).^(1/2*(k-1))...
	        .*(1+k*t/v).^(-1/2*(v+k-2));

	EC(3,:) = a*gamma((v+k-2)/2)/b*(k*t/v).^(1/2*(k-2))...
	        .*(1+k*t/v).^(-1/2*(v+k-2))...
	        .*((v-1)*k*t/v - (k-1));

	EC(4,:) = a^(3/2)*gamma((v+k-3)/2)*2^(-1/2)/b*(k*t/v).^(1/2*(k-3))...
	        .*(1+k*t/v).^(-1/2*(v+k-2))...
	        .*((v-1)*(v-2)*(k*t/v).^2 - (2*v*k-v-k-1)*(k*t/v)...
                + (k-1)*(k-2));
end


% probability
%---------------------------------------------------------------------------
d      = [1:(D + 1)];
CPF    = R(d)*EC(d,:);

% p value assuming Poisson behaviour
%---------------------------------------------------------------------------
CPF    = 1 - exp(-CPF); 
