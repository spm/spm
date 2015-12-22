function [P,f] = spm_gn_fmin(fun,Q,C,varargin)
% objective function minimisation using Gauss-Newton line searches
% FORMAT [P,F] = spm_gn_fmin(fun,Q,C,varargin)
%
% fun - function or inline function f - fun(P,varargin)
% P   - free parameters (prior mean)
% C   - prior covariance
%
% P   - optimised parameters
% f   - optimised value of fun(P)
%
%--------------------------------------------------------------------------
% spm_fmin is a slow but robust function minimiser that uses a Gauss-Newton
% method and successive line searches
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gn_fmin.m 6654 2015-12-22 12:55:36Z spm $
 
 
% stochastic search
%--------------------------------------------------------------------------
P     = spm_vec(Q);
pmin  = P;
n     = length(P);                            % number of parameters
N     = 32;                                   % number of samples
try
    C;                                        % prior covariance
catch
    C = speye(n,n);
end
C     = C + speye(n,n)*exp(-32);
 
% Optimise sampling distribution iteratively
%==========================================================================
spm_figure('GetWin','FMIN');
 
 
% sample objective function using N(P,C)
%----------------------------------------------------------------------
p      = P*ones(1,N) + spm_sqrtm(C)*randn(n,N);
p(:,1) = pmin;
F      = sparse(N,1);
for i = 1:N
    F(i) = feval(fun,spm_unvec(p(:,i),Q),varargin{:});
end
[f,i] = min(F);
pmin  = p(:,i);
clear p F
 
 
t = [-4:4];
for k = 1:8
 
    % gradiants and curvatures
    %----------------------------------------------------------------------
    Q            = spm_unvec(pmin,Q);
    [dfdpp,dfdp] = spm_diff(fun,Q,varargin{:},[1 1]);
    dfdpp        = spm_cat(dfdpp')';
 
    % line search down steepest gradient
    %----------------------------------------------------------------------
    for i = 1:length(t)
        dp     = spm_dx(dfdpp,dfdp',{t(i)});
        p(:,i) = pmin + dp;
        F(i,1) = feval(fun,spm_unvec(p(:,i),Q),varargin{:});
    end
    [f,i]    = min(F);
    pmin     = p(:,i);
    M(k)     = f;
    P(:,k)   = pmin;
 
    % plot objective function and minimisers
    %----------------------------------------------------------------------
    subplot(2,2,1)
    plot(t,F)
    xlabel('rate of decent','FontSize',12)
    title('objective function','FontSize',16)
    axis square
 
    subplot(2,2,2)
    plot(P')
    xlabel('iteration','FontSize',12)
    title('parameter values','FontSize',16)
    axis square
 
    subplot(2,2,3)
    bar(pmin)
    xlabel('parameter','FontSize',12)
    title('minimiser','FontSize',16)
    axis square
 
    subplot(2,2,4)
    bar(M)
    xlabel('iteration','FontSize',12)
    title('minimum','FontSize',16)
    axis square
    drawnow
 
    % convergence
    %----------------------------------------------------------------------
    if k > 1
        if norm(P(:,k) - P(:,k - 1),1)/norm(P(:,k),1) < exp(-8), break, end
    end
 
end
 
% minimiser
%--------------------------------------------------------------------------
P     = spm_unvec(pmin,Q);
