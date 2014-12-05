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
% $Id: spm_ADEM_update.m 6282 2014-12-05 21:57:47Z karl $


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
    
    if nargin > 1
        
        % parameter covariance
        %------------------------------------------------------------------
        np          = length(DEM.M(i).pC);
        qP          = spm_inv(C(1:np,1:np));
        pP          = spm_inv(DEM.M(i).pC);
        pC          = spm_inv(COV*qP + (1 - COV)*pP);
        np          = np + 1;
        C           = C(np:end,np:end);
        DEM.M(i).pC = pC;
        
        % and parameters
        %------------------------------------------------------------------
        pE          = spm_vec(DEM.M(i).pE);
        qE          = spm_vec(DEM.qP.P{i});
        qE          = pC*(COV*qP*qE + (1 - COV)*pP*pE);
        DEM.M(i).pE = spm_unvec(qE,DEM.M(i).pE);
        
    else
        
        % parameters
        %------------------------------------------------------------------
        qE          = spm_vec(DEM.qP.P{i});
        DEM.M(i).pE = spm_unvec(qE,DEM.M(i).pE);
        
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
    DEM.G(i).x = spm_unvec(DEM.pU.x{i}(:,end),DEM.G(i).x);
end
for i = 1:n
    if ~isempty(DEM.G(i).v)
        DEM.G(i).v  = spm_unvec(DEM.pU.v{i}(:,end),DEM.G(i).v);
    end
end
DEM.G(n).a = spm_unvec(DEM.qU.a{n}(:,end),DEM.G(n).a);
