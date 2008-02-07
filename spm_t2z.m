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
% The standard normal distribution approximates Student's
% t-distribution for large degrees of freedom. In univariate
% situations, conventional wisdom states that 30 degrees of freedom is
% sufficient for such an approximation. In the imaging context, the
% multiple comparisons problem places emphasis on the extreme tails of
% the distribution. For PET neuroimaging simulation suggests that 120
% degrees of freedom are required before the distribution of the
% maximal voxel value in a t-statistic image is adequately approximated
% by that of the maxima of a gaussian statistic image (these
% distributions usually being approximated using the theory of
% continuous random fields)  (KJW - private communication). For fMRI
% with it's higher resolution, it is likely that even greater degrees
% of freedom are required for such an approximation.
%
% *No* one-one approximation is made in this code for high df: This is
% because the t2z accuracy reduces as t increases in absolute value
% (particularly in the extrapolation region, underestimating the true
% z. In this case imposing a one-one relationship for df>d say would
% give a jump from df=d-1 to df=d.
%
% For t deviates with very small tail probabilities (< Tol = 10^(-10),
% corresponding to a z of about 6), the corresponding z is computed by
% extrapolation of the t2z relationship z=f(t). This extrapolation
% takes the form of z = log(t-t1+l0) + (z1-log(l0)). Here (t1,z1) is
% the t & z ordinates with tail probability Tol. l0 is chosen such that
% at the point where extrapolation takes over (t1,z1), continuity of
% the first derivative is maintained. Thus, the gradient of the f(t) at
% t1 is estimated as m using six points equally spaced to t1-0.5, and
% l0 is then 1/m.  Experience suggests that this underestimates z,
% especially for ludicrously high t and/or high df, giving conservative
% (though still significant) results.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_t2z.m 1143 2008-02-07 19:33:33Z spm $


%-Initialisation
%===========================================================================

% p-value tolerance: t-values with tail probabilities less than Tol are
%                    `converted' to z by extrapolation
%---------------------------------------------------------------------------
if nargin<3, Tol = 10^(-10); end

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2, error('insufficient arguments'), end
if length(df)~=1, error('df must be a scalar'), end
if df<=0, error('df must be strictly positive'), end

%-Computation
%===========================================================================
z     = zeros(size(t));

%-Mask out t == 0 (z==0) and t == +/- Inf (z==+/- Inf), where
% betainc(1,*,*) and betainc(0,*,*) warn "Log of zero"
%---------------------------------------------------------------------------
Qi    = find(isinf(t));
if length(Qi), z(Qi)=t(Qi); end
tmp   = df./(df + t.^2);
Q     = find(tmp~=1 & ~isinf(t));
if ~length(Q); return; end

%-Mask out at +/- t1 for interpolation
%---------------------------------------------------------------------------
t1    = -spm_invTcdf(Tol,df);
mQb   = abs(t(Q)) > t1;

%-t->z using Tcdf & invNcdf for abs(t)<=t1 
%===========================================================================
if any(~mQb)
    QQnb = find(~mQb);

    %-Compute (smaller) tail probability
    %-Chunk up to avoid convergence problems for long vectors in betacore
    p   = zeros(size(QQnb));
    tmp = [1,[501:500:length(QQnb)],length(QQnb)+1];
    for i = 1:length(tmp)-1
        p(tmp(i):tmp(i+1)-1) = ...
           betainc(df./(df + t(Q(QQnb(tmp(i):tmp(i+1)-1))).^2),df/2,.5)/2;
    end
    
    %-Compute standard normal deviate lower tail prob equal to p
    z(Q(QQnb)) = sqrt(2)*erfinv(2*p - 1);
end


%-Compute standard normal deviates for large t where p-value under/overflows
%-Use logarithmic function for extrapolation, fitted such that first
% derivative is continuous. Estimate gradient from the last 0.5 (t) of
% the (computable) t2z relationship.
%===========================================================================
if any(mQb)
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
    QQ    = Q(mQb); % positions of t-values left to process
    z(QQ) = - ( log( (2*(t(QQ)>0)-1).*t(QQ) -t1 + l0 ) + (z1-log(l0)) );
    %-------------------------------------------------------------------

%   %-------------------------------------------------------------------
%   %-Linear extrapolation
%   %-------------------------------------------------------------------
%   %-adjust c for line through (t1,z1)
%   mc(2)       = z1-mc(1)*t1;
%
%   %-Perform extrapolation, negate positive t-values
%   QQ    = Q(mQb); % positions of t-values left to process
%   z(QQ) = - ( (2*(t(QQ)>0)-1).*t(QQ)*mc(1) + mc(2) );
%   %-------------------------------------------------------------------

end


%-Negate (the negative) z-scores corresponding to positive t-values
%---------------------------------------------------------------------------
z(t>0)=-z(t>0);
