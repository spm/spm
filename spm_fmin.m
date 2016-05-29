function [P,F] = spm_fmin(fun,Q,C,varargin)
% Objective function minimisation
% FORMAT [P,F] = spm_fmin('fun',P,C,varargin)
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
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fmin.m 6801 2016-05-29 19:18:06Z karl $


% stochastic search
%--------------------------------------------------------------------------
P     = spm_vec(Q);
pmin  = P;
n     = length(P);                            % number of parameters
try
    C;                                        % prior covariance
catch
    C = speye(n,n);
end
C     = C + speye(n,n)*exp(-16);

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
    if k == 1
        N = 256;                       % number of samples
    else
        N = 32;
    end
    p      = P*ones(1,N) + spm_sqrtm(C)*randn(n,N);
    p(:,1) = pmin;
    F      = sparse(N,1);
    for i = 1:N
        F(i) = feval(fun,spm_unvec(p(:,i),Q),varargin{:});
    end
    [m,i] = min(F);
    pmin  = p(:,i);

    % supplement with line search along principal eigenvector
    %----------------------------------------------------------------------
    [U,S] = spm_svd(C);
    U     = U(:,1);
    S     = S(1);
    x     = linspace(-3*sqrt(S),3*sqrt(S),16 + 1);
    for i = 1:(16 + 1)
        p(:,end + 1) = pmin + U*x(i);
        F(end + 1,1) = feval(fun,spm_unvec(p(:,end),Q),varargin{:});
    end
    [m,i]  = min(F);
    pmin   = p(:,i);
    M(k)   = m;
    R(:,k) = pmin;

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

    % mean (P)
    %----------------------------------------------------------------------
    q      = exp(-(F - mean(F))/T);
    q      = q/sum(q);
    P      = p*q;

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

end

% minimiser
%--------------------------------------------------------------------------
P     = spm_unvec(pmin,Q);
