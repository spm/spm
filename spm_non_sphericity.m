function [xVi] = spm_non_sphericity(xVi)
% return error covariance constraints for basic ANOVA designs
% FORMAT [xVi] = spm_non_sphericity(xVi)
%
% required fields:
% xVi.I    - n x 4 matrix of factor level indicators
%              I(n,i) is the level of factor i for observation n
% xVi.var  - 1 x 4 vector of flags
%              var(i) = 1; different variance among levels of factor i
% xVi.dep  - 1 x 4 vector of flags
%              dep(i) = 1;      dependencies within levels of factor i
%
% Output:
% xVi.Vi   -  cell of covariance components
% or
% xVi.V    -  speye(n,n)
%
% See also; spm_Ce.m & spm_spm_ui.m
%__________________________________________________________________________
%
% Non-sphericity specification
% =========================================================================
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
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_non_sphericity.m 5219 2013-01-29 17:07:07Z spm $


% create covariance components Q{:}
%==========================================================================
[n,f] = size(xVi.I);                    % # observations, % # Factors
l     = max(xVi.I);                     % levels

% if var(i): add variance component for each level of factor i,
%--------------------------------------------------------------------------
Q     = {};
for i = find(xVi.var)
    for j = 1:l(i)
        u          = xVi.I(:,i) == j;
        q          = spdiags(u,0,n,n);
        Q{end + 1} = q;
    end
end

% effects (discounting factors with dependencies) as defined by interactions
%--------------------------------------------------------------------------
X     = ones(n,1);
for i = find(~xVi.dep & (l > 1))
    Xi    = sparse(1:n,xVi.I(:,i),1,n,l(i));
    Xj    = X;
    X     = sparse(n,0);
    for j = 1:size(Xi,2)
        for k = 1:size(Xj,2)
            X(:,end + 1) = Xi(:,j) & Xj(:,k);
        end
    end
end

% dependencies among repeated measures created by the Hadamard product
%--------------------------------------------------------------------------
for i = find(xVi.dep)
    q     = sparse(1:n,xVi.I(:,i),1,n,l(i));
    P     = q*q';
    for j = 1:size(X,2)
        for k = (j + 1):size(X,2)
            Q{end + 1} = (X(:,j)*X(:,k)' + X(:,k)*X(:,j)').*P;
        end
    end
end

% set Q in non-sphericity structure
%==========================================================================

% if i.i.d nonsphericity (V) is known otherwise there are components {Vi}
%--------------------------------------------------------------------------
if length(Q) > 1
    xVi.Vi = Q;
else
    xVi.V  = speye(n,n);
end
