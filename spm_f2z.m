function [z,f1u,f1l] = spm_f2z(f,df)
% Students f to standard Normal (z-score) distribution
% FORMAT [z,f1u,f1l] = spm_f2z(f,df);
% f   - f values 
% df  - (vector of) degrees of freedom - [df1,df2]
% f1u - upper f-value where linear extrapolation starts
%       empty if no extrapolation
% f1l - lower f-value where linear extrapolation starts
%       empty if no extrapolation
%__________________________________________________________________________
%
% spm_f2z implements a distributional transformation from the F distribution
% to the the standard Gaussian, using incomplete Beta functions and the 
% inverse error function.
%
% Returns z as deviates from the standard Normal (Gaussian) distribution,
% with lower tail probability equal to that of the supplied F statistics,
% on df degrees of freedom.
%
% For very large F deviates with very small upper tail probabilities,
% the corresponding z is computed by linear extrapolation of the f2z
% relationship. For F deviates near zero with very small lower tail
% probabilities, the corresponding z is computed by linear extrapolation
% of the f2z relationship of the reciprocal F statistic, giving a z value
% which must be negated.
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-version control-%
% V1a	- 03/02/95 - 

%===========================================================================

% p-value tolerance: F-values with tail probabilities less than this are
% converted to z by linear extrapolation
%---------------------------------------------------------------------------
Tol         = 10^(-14);

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end
df          = df(:)';
if (length(df)~=2) error('df must be a two-vector'), end
if any(df(:)<=0) error('df must be strictly positive'), end
if any(floor(df(:))~=ceil(df(:))) error('df must be integer'), end

%-Computation
%===========================================================================
z        = zeros(size(f));

% mask out f <= 0 where z = -Inf and betainc(0,a,b) involves log of zero
%---------------------------------------------------------------------------
z(f<=0)  = -Inf*ones(sum(f<=0),1);
Q        = find(f>0);

if ~length(Q); return; end

%-Compute lower tail probability
%---------------------------------------------------------------------------
p        = 1-betainc(df(2)./(df(2)+df(1).*f(Q)),df(2)/2,df(1)/2);

%-Compute standard normal deviate with lower tail prob. equal to p
%---------------------------------------------------------------------------
p_ok = abs(p-0.5) < 0.5-Tol;
if any(p_ok)
	z(Q(p_ok)) = sqrt(2)*erfinv(2*p(p_ok) - 1);
end % if


%-Compute standard normal deviates for large f where p-value overflows.
%-Use linear extrapolation, based on line through low f end point,
% with gradient estimated from lowest 0.5 (f) of the f2z relationship.
%---------------------------------------------------------------------------
p_lge = p>1-Tol;
if any(p_lge)
	f1          = spm_fzero('spm_Fcdf',10,[],0,df,1-Tol);
	z1          = sqrt(2)*erfinv(2*(1-Tol) - 1);
	f2          = f1-[1:5]/10;
	z2          = spm_f2z(f2,df);
	%-least squares line through ([f1,f2],[z1,z2]) : z = m*f + c
	mc          = [[f1,f2]',ones(length([f1,f2]),1)] \ [z1,z2]';
	%-adjust c for line through (f1,z1)
	mc(2)       = z1-mc(1)*f1;
	%-Perform extrapolation
	z(Q(p_lge)) = f(Q(p_lge))*mc(1) + mc(2);
	f1u=f1;
end % if p_lge

%-Compute standard normal deviates for small f where p-value underflows.
%-Do this by linear extrapolation of 1/f values, which have an F
% distribution with [df(2),df(1)] degrees of freedom. 
%---------------------------------------------------------------------------
p_sml = p<Tol;
if any(p_sml)
	[z(Q(p_sml)),f1l] = spm_f2z(1./f(Q(p_sml)),[df(2),df(1)]);
	z(Q(p_sml))       = - z(Q(p_sml));
	f1l               = 1/f1l;
end % if p_sml
