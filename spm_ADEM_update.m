function [DEM] = spm_ADEM_update(DEM,COV)
% Updates ADEM structure using conditional expectations
% FORMAT [DEM] = spm_ADEM_update(DEM,COV)
%
% DEM - DEM structure
% COV - Covariance of parameter (P) fluctuations (E): P(i + 1) = P(i) + E
%     - where cov(E) = COV*pC
%
% This routine updates posterior expectations about states and parameters
% by replacing prior expectations with posterior expectations (and
% similarly updating hidden states and causes to the final iteration). If
% called with an extra argument, the posterior variances of the
% parameters are also updated.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ADEM_update.m 6506 2015-07-24 10:26:51Z karl $


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
    qE          = DEM.qU.x{i};
    if ~isempty(qE)
        DEM.M(i).x = spm_unvec(qE(:,end),DEM.M(i).x);
    end
    
    % parameters
    %----------------------------------------------------------------------
    qE          = spm_vec(DEM.qP.P{i});
    DEM.M(i).pE = spm_unvec(qE,DEM.M(i).pE);
    
    if nargin > 1
        
        % parameter covariance
        %------------------------------------------------------------------
        pC          = DEM.M(i).pC;
        np          = length(pC);
        qC          = C(1:np,1:np);
        C           = C(np + 1:end,np + 1:end);
        DEM.M(i).pC = qC + COV*pC;
        
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
if isfield(DEM.G,'a')
    DEM.G(n).a = spm_unvec(DEM.qU.a{n}(:,end),DEM.G(n).a);
end
