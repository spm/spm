function [x,pnum,pgauss] = spm_mci_postslices (post,M,U,Y,Nbins)
% Univariate slices through posterior density
% FORMAT [x,pnum,pgauss] = spm_mci_postslices (post,M,U,Y,Nbins)
%
% post      posterior data structure
% M,U,Y     as usual
% Nbins     Number of bins per dimension
%
% x         [Np x Nbins] matrix where x(p,:) is domain for pth variable
% pnum      [Np x Nbins] where pnum(p,j) = p(x(p)=xj|x(\p),Y) ie. the posterior
%           density of variable p conditioned on the posterior mean of the other
%           variables. This is estimated numerically from evaluation of log joint
% pgauss    As pnum but under assumption that posterior is multivariate Gaussian
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_postslices.m 6697 2016-01-27 14:57:28Z spm $

try, Nbins=Nbins; catch, Nbins=50; end

k=5; % Defines width of domain in number of SDs

% For computing log prior term
M = spm_mci_priors (M);
M = spm_mci_minit (M);

% Need high tolerances to avoid drop-outs in log-joint plots
M.reltol=1e-6;
M.abstol=1e-6;

Np=length(post.Ep);
Ep=post.Ep;
s=sqrt(diag(post.Cp));
pE=spm_vec(M.pE);
for p=1:Np,
    xmin=Ep(p)-k*s(p);
    xmax=Ep(p)+k*s(p);
    x(p,:)=linspace(xmin,xmax,Nbins);
    P = Ep;
    for j=1:Nbins,
        % Get parameters in eigenspace of prior
        P(p) = x(p,j);
        par = M.V'*(P-pE);
        eq(j) = spm_mci_joint (par,M,U,Y);
        pg(j) = spm_mvNpdf(P,Ep,post.Cp);
    end
    eq=exp(eq);
    pnum(p,:)=eq/sum(eq);
    pgauss(p,:)=pg/sum(pg);
end

figure
for i=1:Np,
    subplot(Np,1,i);
    plot(x(i,:),pgauss(i,:),'r');
    hold on
    plot(x(i,:),pnum(i,:),'k');
    legend('Gaussian','Numeric');
    xlabel(sprintf('P(%d)',i));
end