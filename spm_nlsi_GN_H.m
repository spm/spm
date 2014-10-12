function [Ep,Cp,Eh,F] = spm_nlsi_GN_H(M,U,Y)
% Bayesian inversion of a nonlinear model with hierarchical optimisation
% FORMAT [Ep,Cp,Eh,F] = spm_nlsi_GN_H(M,U,Y)
%
% Dynamical MIMO models
%__________________________________________________________________________
%
% M.IS - function name f(P,M,U) - generative model
%        This function specifies the nonlinear model:
%        y = Y.y = IS(P,M,U) + X0*P0 + e
%        were e ~ N(0,C).  For dynamic systems this would be an integration
%        scheme (e.g. spm_int). spm_int expects the following:
%
%     M.f  - f(x,u,P,M)
%     M.g  - g(x,u,P,M)
%       x  - state variables
%       u  - inputs or causes
%       P  - free parameters
%       M  - fixed functional forms and parameters in M
%
% M.FS - function name f(y,M)   - feature selection
%        This [optional] function performs feature selection assuming the
%        generalized model y = FS(y,M) = FS(IS(P,M,U),M) + X0*P0 + e
%
% M.P  - starting estimates for model parameters [optional]
%
% M.pE - prior expectation  - E{P}   of model parameters
% M.pC - prior covariance   - Cov{P} of model parameters
%
% M.hE - prior expectation  - E{h}   of log-precision parameters
% M.hC - prior covariance   - Cov{h} of log-precision parameters
%
% U.u  - inputs (or just U)
% U.dt - sampling interval
%
% Y.y  - outputs (samples x observations)
% Y.dt - sampling interval for outputs
% Y.X0 - Confounds or null space      (over size(y,1) bins or all vec(y))
% Y.Q  - q error precision components (over size(y,1) bins or all vec(y))
%
%
% Parameter estimates
%--------------------------------------------------------------------------
% Ep  - (p x 1)         conditional expectation    E{P|y}
% Cp  - (p x p)         conditional covariance     Cov{P|y}
% Eh  - (q x 1)         conditional log-precisions E{h|y}
%
% log evidence
%--------------------------------------------------------------------------
% F   - [-ve] free energy F = log evidence = p(y|f,g,pE,pC) = p(y|m)
%
%__________________________________________________________________________
% This is the same as spm_nlsi_GN but tries to model the free energy as a
% function of conditional expectations using a sparse mixture of scaled
% Gaussians. The objective is to account for local maxima when optimising
% free energy by recasting the problem in terms of a parameterised mapping 
% from conditional expectations to free energy explicitly.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_nlsi_GN_H.m 6233 2014-10-12 09:43:50Z karl $
 
% dimension reduction of parameter space
%--------------------------------------------------------------------------
pC    = M.pC;
pE    = spm_vec(M.pE);
if isstruct(pC)
    pC = spm_vec(pC);
end
if isvector(pC);
    pC = diag(pC);
end
V     = spm_svd(pC);
pC    = V'*pC*V;
pE    = V'*pE;
np    = size(V,2);
 
% priors on Gaussian mixture model
%--------------------------------------------------------------------------
ppE(1).a = 0;
ppE(1).m = spm_vec(pE);
ppE(1).c = pC;
 
for i = 2:4
    ppE(i).a = -2;
    ppE(i).m = spm_vec(pE) + spm_sqrtm(pC)*sparse(1 + rem(i,np),1,1,np,1);
    ppE(i).c = pC;
end
 
% set up second order model
%--------------------------------------------------------------------------
MM.IS = 'spm_ho_gm';
MM.hE = 8;
 
 
% hierarchical search
%--------------------------------------------------------------------------
npp   = 64;
for i = 1:4
    
    
    % generate potential search points (PP)
    %----------------------------------------------------------------------
    MM.pE = ppE;
    for j = 1:npp
        PP(:,j) = ppE(1).m + spm_sqrtm(ppE(1).c)*randn(np,1);
    end
    
    % evaluate F over selected points (S) to generate YY.y
    %----------------------------------------------------------------------
    M.Nmax    = 1;
    M.nograph = 1;
    for j = 1:npp
        M.P                 = spm_unvec(V*PP(:,j),M.pE);
        [Ep,Cp,Eh,F,L,dFdp] = spm_nlsi_GN(M,U,Y);
        YY.y(j,:)           = spm_vec(F,dFdp)';
    end
    
    % adjust F
    %----------------------------------------------------------------------
    YY.y(:,1) = YY.y(:,1) - max(YY.y(:,1));
    YY.y      = YY.y/std(YY.y(:,1));
    
    % invert second level model to optimise PP
    %----------------------------------------------------------------------
    ppE       = spm_nlsi_GN(MM,PP,YY);
 
    % sort mixture
    %----------------------------------------------------------------------
    [a,j]     = sort(spm_cat({ppE.a}),'descend');
    j         = j(a > (max(a) - 3));
    ppE       = ppE(j);
    
    
    % show progress results
    %----------------------------------------------------------------------
    spm_figure('GetWin','Hierarchical optimisation')
    x     = ppE(1).m(1);
    d     = 8*sqrt(ppE(1).c(1,1));
    xd    = linspace(x - d,x + d,32);
    y     = ppE(1).m(2);
    d     = 8*sqrt(ppE(1).c(2,2));
    yd    = linspace(y - d,y + d,32);
    [x,y] = meshgrid(xd,yd);
    
    XX      = ppE(1).m*ones(1,32*32);
    XX(1,:) = spm_vec(x)';
    XX(2,:) = spm_vec(y)';
    
    F     = spm_ho_gm(ppE,MM,XX,'eval');
    F     = reshape(F,32,32);
    
    subplot(2,2,1),mesh(x,y,F), axis square,     title('Free energy landscape')
    subplot(2,2,2),imagesc(xd,yd,F), axis square,title('Free energy landscape')
    
    for j = 1:size(XX,2)
        M.P       = spm_unvec(V*XX(:,j),M.pE);
        [~,~,~,F] = spm_nlsi_GN(M,U,Y);
        SF(j)     = F;
    end
    SF    = reshape(SF,32,32);
    subplot(2,2,3),mesh(x,y,SF), axis square,     title('Free energy landscape')
    subplot(2,2,4),imagesc(xd,yd,SF), axis square,title('Free energy landscape')
    
end
 
% final gradient descent
%--------------------------------------------------------------------------
M.Nmax       = 32;
M.nograph    = 0;
M.P          = V*ppE(1).m;
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);
