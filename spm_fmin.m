function [P F] = spm_fmin(fun,P,C,varargin)
% objective function minimisation
% FORMAT [P F] = spm_fmin('fun',P,C,varargin)
%
% fun - function or inline function f - fun(P,varargin)
% P   - free parameters (prior mean)
% C   - prior covariance
%
% P   - optimised parameters
% f   - optimised value of fun(P)
%
%--------------------------------------------------------------------------
% spm_fmin is a slow but robust function minimiser that uses a stochastic
% sampling of the objective function to be minimised (supplemented by a line
% search along the principal eigenvariate at the current sampling density.
% The sampling density is approximated with a Gaussian (first and second
% order moments) using that the sampling density is:
%
%           p(P) = (1/Z)*exp(-fun(P)/T)
%
% where the temperature; T is the sample standard deviation of the sampled
% objective function.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fmin.m 1961 2008-07-26 09:38:46Z karl $


% stochastic search
%--------------------------------------------------------------------------
n     = length(spm_vec(P));                   % number of parameters
N     = 128;                                  % number of samples
try
    C;                                        % prior covariance
catch
    C = speye(n,n);
end

% Optimise sampling distribution iteratively
%==========================================================================
spm_figure('GetWin','FMIN');
for k = 1:8

    % report - current sample density
    %----------------------------------------------------------------------
    p1    = linspace(-3*sqrt(C(1,1)),3*sqrt(C(2,2)),64);
    p2    = linspace(-3*sqrt(C(2,2)),3*sqrt(C(2,2)),64);
    iC    = inv(C(1:2,1:2));
    for i = 1:64
        for j = 1:64
            d(i,j) = exp(-[p1(i) p2(j)]*iC*[p1(i) p2(j)]'/2);
        end
    end
    subplot(3,2,1)
    imagesc(p2 + P(2),p1 + P(1),d)
    xlabel('2nd parameter','FontSize',12)
    ylabel('1st parameter','FontSize',12)
    title(sprintf('iteration %i',k - 1),'FontSize',16)
    axis square xy

    % sample objective function using N(P,C)
    %----------------------------------------------------------------------
    p     = spm_vec(P)*ones(1,N) + spm_sqrtm(C)*randn(n,N);
    try
        p(:,1) = pm;
    end
    F     = sparse(N,1);
    for i = 1:N
        F(i) = feval(fun,spm_unvec(p(:,i),P),varargin{:});
    end
    [m i] = min(F);
    P     = spm_unvec(p(:,i),P);

    % supplement with line search along principal eigenvector
    %----------------------------------------------------------------------
    [U S] = spm_svd(C);
    U     = U(:,1);
    S     = S(1);
    x     = linspace(-3*sqrt(S),3*sqrt(S),16 + 1);
    for i = 1:(16 + 1)
        p(:,end + 1) = spm_vec(P) + U*x(i);
        F(end + 1,1) = feval(fun,spm_unvec(p(:,end),P),varargin{:});
    end
    [m i] = min(F);
    pm    = p(:,i);
    M(k)  = m;

    % plot objective function along eigenvariate and sampled values
    %----------------------------------------------------------------------
    subplot(3,2,3)
    plot(x,F(end - 16:end),'.')
    xlabel('first eigenvariate','FontSize',12)
    title('objective function','FontSize',16)
    axis square

    subplot(3,2,5)
    plot(sort(F))
    xlabel('sample','FontSize',12)
    title('ranked samples','FontSize',16)
    axis square


    % Laplace approximation: p(P)
    %======================================================================

    % temperature
    %----------------------------------------------------------------------
    T      = sqrt(2)*std(F);

    % mean (and record in R)
    %----------------------------------------------------------------------
    q      = exp(-(F - mean(F))/T);
    q      = q/sum(q);
    P      = p*q;
    R(:,k) = P;

    % dispersion
    %----------------------------------------------------------------------
    for i = 1:n
        for j = 1:n
            C(i,j) = ((p(i,:) - P(i)).*(p(j,:) - P(j)))*q;
        end
    end

    % report - updated sampling density
    %----------------------------------------------------------------------
    subplot(3,2,2)
    iC      = inv(C(1:2,1:2));
    for i = 1:64
        for j = 1:64
            d(i,j) = exp(-[p1(i) p2(j)]*iC*[p1(i) p2(j)]'/2);
        end
    end
    imagesc(p2 + P(2),p1 + P(1),d)
    title(sprintf('iteration %i',k),'FontSize',16)
    axis square xy

    % superimpose line search and plot means and min(F)
    %----------------------------------------------------------------------
    subplot(3,2,1),hold on
    plot(p(2,end - 16:end),p(1,end - 16:end),'r'), hold off
    subplot(3,2,2),hold on
    plot(p(2,end - 16:end),p(1,end - 16:end),'r'), hold off

    subplot(3,2,4)
    plot(R')
    xlabel('iteration','FontSize',12)
    title('parameter values','FontSize',16)
    axis square

    subplot(3,2,6)
    bar(M)
    xlabel('iteration','FontSize',12)
    title('minimum','FontSize',16)
    axis square
    drawnow

    % convergence
    %----------------------------------------------------------------------
    if k > 1
        if norm(R(:,k) - R(:,k - 1),1)/norm(R(:,k),1) < exp(-8), break, end
    end

end

% minimiser
%---------------------------------=----------------------------------------
[m i] = min(F);
P     = spm_unvec(p(:,i),P);
