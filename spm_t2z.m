function [z,t1,z1] = spm_t2z(t,df,Tol)
% Students t to standard Normal (z-score) distribution
% FORMAT [z,t1] = spm_t2z(t,df,Tol);
% t   - t values 
% df  - degrees of freedom
% Tol - minimum tail probability for direct computation
%       Defaults to 10^(-16), a z of about 8.2
% t1  - (absolute) t-value where linear extrapolation starts
%       empty if no extrapolation
% z1  - Equivalent standard Normal ordinate to t-value t1
%__________________________________________________________________________
%
% spm_t2z implements a distributional transformation from the Student's
% t to the unit Gaussian using incomplete Beta functions and the
% inverse error function.
%
% Returns z as deviates from the standard Normal (Gaussian)
% distribution with lower tail probability equal to that of the
% supplied t statistics with df degrees of freedom.
%
% For t deviates with very small tail probabilities (< Tol), the
% corresponding z is computed by extrapolation of the t2z relationship
% z=f(t). This extrapolation takes the form of z = log(t-t1+l0) +
% (z1-log(l0)). Here (t1,z1) is the t & z ordinates with tail
% probability Tol. l0 is chosen such that at the point where
% extrapolation takes over (t1,z1), continuity of the first derivative
% is maintained. Thus, the gradient of the f(t) at t1 is estimated as m
% using six points equally spaced to t1-0.5, and l0 is then 1/m.
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-version notes-%
% 16/12/93 - 
% 01/06/94 - Stripped down for SPM
% 03/02/95 - Included extrapolation and argument checks
% 20/07/95 - Altered t==0 checking to allow for overflow to 1 in computation 
%            of betainc argument - this was generating "log of zero" 
%            errors - annoying but not serious.
% 23/08/96 - Changed to logarithmic extrapolation
%            - fMRI was giving rediculous Z values!



%-Initialisation
%===========================================================================

% p-value tolerance: t-values with tail probabilities less than this are
% converted to z by linear extrapolation
%---------------------------------------------------------------------------
if nargin<3, Tol = 10^(-16); end

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end
if (length(df)~=1) error('df must be a scalar'), end
if df<=0 error('df must be strictly positive'), end
if df > 32; df = round(df); end



%-Computation
%===========================================================================
z            = zeros(size(t));

%-Mask out t == 0 (z==0) where betainc(1,*,*) involves log of zero
% (betainc(0,*,*) involves log(0) too, but t == +/- Inf to get this!)
%---------------------------------------------------------------------------
tmp          = df./(df + t.^2);
Q            = find(tmp~=1);
if ~length(Q); return; end

%-Compute (smaller) tail probability
%---------------------------------------------------------------------------
p           = betainc(tmp(Q),df/2,.5)/2;


%-Compute standard normal deviate lower tail prob equal to p
%---------------------------------------------------------------------------
p_ok = p > Tol;
if any(p_ok)
	z(Q(p_ok)) = sqrt(2)*erfinv(2*p(p_ok) - 1);
end % if


%-Compute standard normal deviates for large t where p-value under/overflows
%-Use logarithmic function for extrapolation, fitted such that first
% derivative is continuous. Estimate gradient from the last 0.5 (t) of
% the (computable) t2z relationship.
%---------------------------------------------------------------------------
if any(~p_ok)
	t1          =-spm_fzero('spm_Tcdf',-10,[],0,df,Tol);
	z1          =-sqrt(2)*erfinv(2*Tol-1);
	t2          =t1-[1:5]/10;
	z2          =spm_t2z(t2,df);
	%-least squares line through ([f1,t2],[z1,z2]) : z = m*f + c
	mc          = [[t1,t2]',ones(length([t1,t2]),1)] \ [z1,z2]';

	%-------------------------------------------------------------------
	%-Logarithmic extrapolation
	%-------------------------------------------------------------------
	l0=1/mc(1);
	%-Perform logarithmic extrapolation, negate z for positive t-values
	Q    = Q(~p_ok); % positions of t-values left to process
	z(Q) = - ( log( (2*(t(Q)>0)-1).*t(Q) -t1 + l0 ) + (z1-log(l0)) );
	%-------------------------------------------------------------------

%	%-------------------------------------------------------------------
%	%-Linear extrapolation
%	%-------------------------------------------------------------------
%	%-adjust c for line through (t1,z1)
%	mc(2)       = z1-mc(1)*t1;
%
%	%-Perform extrapolation, negate positive t-values
%	Q           = Q(~p_ok); % positions of t-values left to process
%	z(Q)        = - ( (2*(t(Q)>0)-1).*t(Q)*mc(1) + mc(2) );
%	%-------------------------------------------------------------------

end


%-Negate (the negative) z-scores corresponding to positive t-values
%---------------------------------------------------------------------------
z(t>0)=-z(t>0);
