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
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ADEM_update.m 6198 2014-09-25 10:38:48Z karl $
 
% update states and parameters (model)
%--------------------------------------------------------------------------
n     = length(DEM.M);
C     = DEM.qP.C;
for i = 1:(n - 1)
    
    % expectations
    %----------------------------------------------------------------------
    DEM.M(i).x  = spm_unvec(DEM.qU.x{i}(:,end),DEM.M(i).x);
    DEM.M(i).pE = DEM.qP.P{i};
    
    % and covariance if required
    %----------------------------------------------------------------------
    if nargin > 1
        np           = length(DEM.M(i).pC);
        DEM.M(i).pC  = C(1:np,1:np);
        np           = np + 1;
        C            = C(np:end,np:end);
    end
    
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
    DEM.G(i).x  = spm_unvec(DEM.pU.x{i}(:,end),DEM.G(i).x);
end
for i = 1:n
    if ~isempty(DEM.G(i).v)
        DEM.G(i).v  = spm_unvec(DEM.pU.v{i}(:,end),DEM.G(i).v);
    end
end
DEM.G(n).a      = spm_unvec(DEM.qU.a{n}(:,end),DEM.G(n).a);
