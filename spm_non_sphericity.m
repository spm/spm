function [xVi] = spm_non_sphericity(xVi)
% return error covariance constraints for basic ANOVA designs
% FORMAT [xVi] = spm_non_sphericity(xVi)
%
% required fields:
% xVi.I    - n x 4 matrix of factor level indicators
%              I(n,i) is the level of factor i for observation n
% xVi.var  - 1 x 4 vector of flags
%              var(i) = 1; levels of factor i have unequal error variances
% xVi.dep  - 1 x 4 vector of flags
%              dep(i) = 1; levels of factor i induce random effects
%
% Output:
% xVi.Vi   -  cell of covariance components
% or
% xVi.V    -  speye(n,n)
%
% See also; spm_Ce.m & spm_spm_ui.m
%___________________________________________________________________________
% Non-sphericity specification
% =========================
%
% In some instances the i.i.d. assumptions about the errors do not hold:
%
% Identity assumption:
% The identity assumption, of equal error variance (homoscedasticity), can
% be violated if the levels of a factor do not have the same error variance.
% For example, in a 2nd-level analysis of variance, one contrast may be scaled
% differently from another.  Another example would be the comparison of
% qualitatively different dependant variables (e.g. normals vs. patients).  If
% You say no to identity assumptions, you will be asked whether the error
% variance is the same over levels of each factor.  Different variances
% (heteroscedasticy) induce different error covariance components that
% are estimated using restricted maximum likelihood (see below).
%
% Independence assumption.
% In some situations, certain factors may contain random effects.  These induce
% dependencies or covariance components in the error terms.   If you say no 
% to independence assumptions, you will be asked whether random effects
% should be modelled for each factor.  A simple example of this would be
% modelling the random effects of subject.  These cause correlations among the
% error terms of observation from the same subject.  For simplicity, it is 
% assumed that the random effects of each factor are i.i.d. 
%
% ReML
% The ensuing covariance components will be estimated using ReML in spm_spm
% (assuming the same for all responsive voxels) and used to adjust the 
% statistics and degrees of freedom during inference. By default spm_spm
% will use weighted least squares to produce Gauss-Markov or Maximum
% likelihood estimators using the non-sphericity structure specified at this 
% stage. The components will be found in xX.xVi and enter the estimation 
% procedure exactly as the serial correlations in fMRI models.
% 
%___________________________________________________________________________
% %W% Karl Friston %E%

% create covariance components Q{:}
%===========================================================================
[n f] = size(xVi.I);			% # observations, % # Factors
l     = max(xVi.I);                     % levels
Q     = {};

% if i.i.d
%---------------------------------------------------------------------------
if ~any(xVi.var) & ~any(xVi.dep)
    xVi.V  = speye(n,n);
    return
end


% error components
%---------------------------------------------------------------------------
for i = 1:f
    
    % add variance component for level j of factor i
    %-----------------------------------------------------------------------
    if xVi.var(i) & (l(i) > 1)
        for j = 1:l(i)
            u          = xVi.I(:,i) == j;
            q          = spdiags(u,0,n,n);
            Q{end + 1} = q;
        end
    end
end

% unless all repeated measures are identically distributed
%---------------------------------------------------------------------------
if ~length(Q)
    Q{end + 1} = speye(n,n);
end

% effects (discounting factors with dependencies)
%---------------------------------------------------------------------------
X    = [];
for i = 1:f
    if (l(i) > 1) & ~xVi.dep(i)
        q     = sparse(1:n,xVi.I(:,i),1,n,l(i));
        X     = [X q];
    end
end

% dependencies among repeated measures
%---------------------------------------------------------------------------
for i = 1:f
    if xVi.dep(i) & (l(i) > 1)
        q     = sparse(1:n,xVi.I(:,i),1,n,l(i));
        P     = q*q';
        for j = 1:size(X,2)
            for k = (j + 1):size(X,2)
                Q{end + 1} = (X(:,j)*X(:,k)' + X(:,k)*X(:,j)').*P;
            end
        end
    end
end

% set Q in non-sphericity structure
%---------------------------------------------------------------------------
xVi.Vi = Q;
