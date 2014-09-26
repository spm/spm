function [DCM,S] = spm_large_dcm_reduce(DCM)
% Optimises the number of prior connectivity eigenmodes
% FORMAT [DCM,S] = spm_large_dcm_reduce(DCM)
% DCM    - DCM structure
% S      - log-evidences
%
% This routine optimises the number of eigenmodes of the prior covariance
% matrix using the eigenvectors of the functional connectivity matrix. The
% optimisation uses post hoc model reduction.
%__________________________________________________________________________
%
% Reference
%
% M.L. Seghier and K.J. Friston, "Network discovery with large DCMs".
% NeuroImage, 68:181-191, 2013.
%__________________________________________________________________________
% Copyright (C) 2012-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_large_dcm_reduce.m 6204 2014-09-26 18:12:27Z mohamed $


%-Create priors
%==========================================================================

% eigenvector constraints on pC for large models
%--------------------------------------------------------------------------
n     = size(DCM.M.pE.A,1);
j     = 1:(n*n);

% remove confounds and find principal modes
%--------------------------------------------------------------------------
y     = DCM.Y.y - DCM.Y.X0*(pinv(DCM.Y.X0)*DCM.Y.y);
V     = spm_svd(y');

for i = 1:n
    
    % remove minor modes from priors on A
    %----------------------------------------------------------------------
    v       = V(:,1:i);
    v       = kron(v*v',v*v');   
    rC{i}   = v*DCM.M.pC(j,j)*v';
    
end


%-Loop over prior covariances to get log-evidences
%==========================================================================

% Get priors and posteriors
%--------------------------------------------------------------------------
qE    = DCM.Ep;
qC    = DCM.Cp;
pE    = DCM.M.pE;
pC    = DCM.M.pC;

% Remove (a priori) null space
%--------------------------------------------------------------------------
U     = spm_svd(pC);
qE    = U'*spm_vec(qE);
pE    = U'*spm_vec(pE);
qC    = U'*qC*U;
pC    = U'*pC*U;

% model search
%--------------------------------------------------------------------------
S     = zeros(1,n);
pc      = DCM.M.pC;
for i = 1:n
    pc(j,j) = rC{i} ;
    S(i) = spm_log_evidence(qE,qC,pE,pC,pE,U'*pc*U);
end

% model evidence
%--------------------------------------------------------------------------
S     = S - min(S);
p     = exp(S - max(S));
p     = p/sum(p);

% Show results
%--------------------------------------------------------------------------
Fgraph = spm_figure('FindWin','Graphics');

if ~isempty(Fgraph)
    subplot(2,2,1,'Parent',Fgraph)
    bar(S)
    title('log-posterior','FontSize',16)
    xlabel('number of prior eigenmodes','FontSize',12)
    ylabel('log-probability','FontSize',12)
    set(gca,'YLim',[(max(S) - 1000), max(S)])
    axis square
    
    subplot(2,2,2,'Parent',Fgraph)
    bar(p)
    title('model posterior','FontSize',16)
    xlabel('number of prior eigenmodes','FontSize',12)
    ylabel('probability','FontSize',12)
    axis square
    drawnow
end


%-Get posterior density of reduced model
%==========================================================================

% Get full priors and posteriors
%--------------------------------------------------------------------------
[p,i] = max(p);
qE    = DCM.Ep;
qC    = DCM.Cp;
pE    = DCM.M.pE;
pC    = DCM.M.pC;

pc(j,j) = rC{i} ;
[F,Ep,Cp] = spm_log_evidence_reduce(qE,qC,pE,pC,pE,pc);

% Bayesian inference and variance
%--------------------------------------------------------------------------
T        = full(spm_vec(pE));
Pp       = spm_unvec(1 - spm_Ncdf(T,abs(spm_vec(Ep)),diag(Cp)),Ep);
Vp       = spm_unvec(diag(Cp),Ep);

% Store parameter estimates
%--------------------------------------------------------------------------
DCM.M.pC = pc;
DCM.Ep   = Ep;
DCM.Cp   = Cp;
DCM.Pp   = Pp;
DCM.Vp   = Vp;

% and in DEM format
%--------------------------------------------------------------------------
DCM.qP.P{1} = Ep;
DCM.qP.C    = Cp;
DCM.qP.V{1} = spm_unvec(diag(Cp),Ep);

% approximations to model evidence: negative free energy, AIC, BIC
%--------------------------------------------------------------------------
DCM.F    = F;
% evidence = spm_dcm_evidence(DCM);
% DCM.AIC  = evidence.aic_overall;
% DCM.BIC  = evidence.bic_overall;
