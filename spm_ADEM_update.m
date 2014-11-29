function [DEM] = spm_ADEM_update(DEM,COV)
% Updates ADEM structure using conditional expectations
% FORMAT [DEM] = spm_ADEM_update(DEM,COV)
%
% DEM - DEM structure
% COV - flag for Bayesian belief updating (with covariance)
%
% this routine updates posterior expectations about states and parameters
% by replacing prior expectations with posterior expectations (and
% similarly updating hidden states and causes to the final iteration). It
% called with an extra argument, the posterior variances of the
% parameters are also updated.
%
% COV ranges from 0 to 1 to conrol the degree of covariance updating.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ADEM_update.m 6270 2014-11-29 12:04:48Z karl $


% preliminaries
%--------------------------------------------------------------------------
if nargin < 2, COV = 0; end

% update states and parameters (model)
%--------------------------------------------------------------------------
n     = length(DEM.M);
C     = DEM.qP.C;
for i = 1:(n - 1)
    
    % states
    %----------------------------------------------------------------------
    DEM.M(i).x  = spm_unvec(DEM.qU.x{i}(:,end),DEM.M(i).x);
    
    
    % parameters
    %----------------------------------------------------------------------
    DEM.M(i).pE = DEM.qP.P{i};
    
    % and parameter covariance
    %----------------------------------------------------------------------
    np          = length(DEM.M(i).pC);
    Q           = spm_inv(C(1:np,1:np));
    P           = spm_inv(DEM.M(i).pC);
    DEM.M(i).pC = spm_inv(COV*Q + (1 - COV)*P);
    np          = np + 1;
    C           = C(np:end,np:end);

    
end
for i = 1:n
    if ~isempty(DEM.M(i).v)
        DEM.M(i).v  = spm_unvec(DEM.qU.v{i}(:,end),DEM.M(i).v);
    end
end

% update states and action (process)
%--------------------------------------------------------------------------
n     = length(DEM.G);
for i = 1:(n - 1)
    DEM.G(i).x = spm_unvec(DEM.pU.x{i}(:,end),DEM.G(i).x);
end
for i = 1:n
    if ~isempty(DEM.G(i).v)
        DEM.G(i).v  = spm_unvec(DEM.pU.v{i}(:,end),DEM.G(i).v);
    end
end
DEM.G(n).a = spm_unvec(DEM.qU.a{n}(:,end),DEM.G(n).a);
