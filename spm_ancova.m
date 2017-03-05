function [F,df,beta,xX,xCon] = spm_ancova(xX,V,Y,c)
% estimation and inference of a linear model
% FORMAT [F,df,beta,xX,xCon] = spm_ancova(xX,V,Y,c);
%
% xX    - (m x p) Design matrix or structure
% V     - (m x m) error covariance constraint
% Y     - {m x n} matrix of response {m x 1} variables
% c     - {p x q} matrix of (q) contrasts
%
% F     - {t x n} matrix of T or F values
% df    - {1 x 2} vector of degrees of freedom
% beta  - {p x n} matrix of parameter estimates
% xX    - design matrix structure
% xCon  - contrast structure
%__________________________________________________________________________
%
% spm_ancova uses a General Linear Model of the form:
%
%   Y  = X*beta + K*e
%
% to compute the parameter estimates (beta) and make inferences (T or F)
% where V = K*K' represents the correlation structure. If c has only one
% column T statistics are returned, otherwise F ratios are computed.
%__________________________________________________________________________
% Copyright (C) 2002-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ancova.m 7033 2017-03-05 11:19:18Z karl $

% assume V = I (i.i.d.) if empty
%--------------------------------------------------------------------------
if isempty(V)
    V          = speye(size(Y,1));
end

% or scalar
%--------------------------------------------------------------------------
if isscalar(V)
    V          = speye(size(Y,1))*V;
end

% create design matrix structure if necessary
%--------------------------------------------------------------------------
if ~isstruct(xX)
    xX         = spm_sp('Set',xX);
end
if ~isfield(xX,'pX')
    xX.pX      = spm_sp('x-',xX);
end


% estimate parameters and sum of squared residuals
%--------------------------------------------------------------------------
beta           = xX.pX*Y;             %-Parameter estimates
res            = spm_sp('r',xX,Y);    %-Residuals
ResSS          = sum(res.^2);         %-Res sum-of-squares


% default contrast
%--------------------------------------------------------------------------
if nargin < 4
    c = speye(size(beta,1));
end

% contrast
%--------------------------------------------------------------------------
xCon           = spm_FcUtil('Set','','F','c',c,xX);
h              = spm_FcUtil('Hsqr',xCon,xX);
X1o            = spm_FcUtil('X1o',xCon,xX);

% ensure trace(V) = m and get traces
%--------------------------------------------------------------------------
V              = V*length(V)/trace(V);
[trRV, trRVRV] = spm_SpUtil('trRV',xX,V);
[trMV, trMVMV] = spm_SpUtil('trMV',X1o,V);
df             = [trMV^2/trMVMV trRV^2/trRVRV];


% F statistics
%--------------------------------------------------------------------------
F              = sum((h*beta).^2,1)./(ResSS*trMV/trRV);

if size(c,2) == 1

    % T statistics
    %----------------------------------------------------------------------
    F          = sqrt(F).*sign(c'*beta);

end 
