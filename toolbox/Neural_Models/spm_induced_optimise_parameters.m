function [Ep M] = spm_induced_optimise_parameters(PARAMS)
% Demo routine that optimises free parameters
%==========================================================================
%
% This exemplar routine illustrates how one can adjust or tune prior
% parameter expectations to produce desired spectral responses as specified
% by the complex eigenvalue spectrum - or a reduced form that considers a
% small number of complex values (roots).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_induced_optimise_parameters.m 5038 2012-11-06 20:25:24Z karl $
 
 
% Model specification
%==========================================================================
 
% options
%--------------------------------------------------------------------------
J     = [2 4];                % indices of hidden states producing outputs
Nm    = 3;                    % number of (complex) eigenmodes to consider
Nc    = 1;
Ns    = 1;
options.spatial  = 'LFP';
options.model    = 'CMM';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
M.J            = J;
M.Nm           = Nm;
 
 
% get priors
%--------------------------------------------------------------------------
[pE pC] = spm_dcm_neural_priors({0 0 0},{},1,options.model);
[pE pC] = spm_L_priors(M.dipfit,pE,pC);
[pE pC] = spm_ssr_priors(pE,pC);
[x,f]   = spm_dcm_x_neural(pE,options.model);
Ep      = pE;
 
 
% orders and model
%==========================================================================
nx      = length(spm_vec(x ));
np      = length(spm_vec(pE));
nu      = Ns;
u       = sparse(1,nu);
 
 
% fix priors if a subset of parameters are specified
%--------------------------------------------------------------------------
try
    V = [];
    for i = 1:length(PARAMS)
        V = [V; spm_fieldindices(pE,PARAMS{i})];
    end
    V   = sort(V);
catch
    V   = 1:np;
end
V   = sparse(V,V,1,np,np);
 
 
% create LFP model
%--------------------------------------------------------------------------
M.f    = f;
M.g    = 'spm_gx_erp';
M.x    = x;
M.n    = nx;
M.pE   = pE;
M.pC   = pC;
M.hE   = 8;
M.hC   = exp(-16);
M.m    = nu;
M.l    = Nc;
 
% solve for steady state
%--------------------------------------------------------------------------
M.x    = spm_dcm_neural_x(pE,M);
M.u    = u;
M.Hz   = 4:64;
M.Nmax = 32;
 
 
% Characterise contributions of parameters to spectral representation
%==========================================================================
[dSdp,s] = spm_diff('spm_ssm2s',pE,M,[],1);
 
% number of complex (oscillatory) eigenmodes
%--------------------------------------------------------------------------
m      = M.Nm;
n      = numel(J);
 
 
% Show results with current (prior) parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Spectral responses'); clf
set(gcf,'UserData',pE);
 
% and plot spectra
%--------------------------------------------------------------------------
subplot(4,2,1)
[g,w]  = spm_s2csd(s,M.Hz);
plot(w,g)
title('Spectral modes','FontSize',16)
xlabel('Frequency')
ylabel('Power')
 
subplot(4,2,3), cla
for i = 1:n
    pE.J       = sparse(1,J(i),1,1,nx);
    G          = spm_csd_mtf(pE,M,[]);
    plot(w,G{1}), hold on
end
title('Spectral respones','FontSize',16)
xlabel('Frequency')
ylabel('Spectral density')
hold off
 
% roots
%--------------------------------------------------------------------------
subplot(4,2,2)
plot(s(1:m),'o'), hold on
title('Spectral density','FontSize',16)
xlabel('Real (amplitude)')
ylabel('Imaginary (frequency)')
 
% wiegths
%--------------------------------------------------------------------------
subplot(4,2,4)
bar(reshape(s(m + 1:end),n,m)')
title('Wiegths','FontSize',16)
xlabel('Relative weight')
xlabel('Mode')
 
 
% and partial derivataives
%--------------------------------------------------------------------------
subplot(4,2,5)
imagesc(real(dSdp(1:m,:)))
title('Contribution (real) (right click)','FontSize',16)
xlabel('Parameter')
ylabel('Mode')
set(get(gca,'Children'),'ButtonDownFcn','spm_opt_bfun')
 
subplot(4,2,7)
imagesc(imag(dSdp(1:m,:)))
title('(Imaginary)','FontSize',16)
xlabel('Parameter')
ylabel('Mode')
set(get(gca,'Children'),'ButtonDownFcn','spm_opt_bfun')
 
subplot(2,2,4)
imagesc(real(dSdp(m + 1:end,:)))
title('Weights','FontSize',16)
xlabel('Parameter')
ylabel('Mode')
set(get(gca,'Children'),'ButtonDownFcn','spm_opt_bfun')
 
 
Ep = pE; return
 
% Optimisation: Target (Y)
%==========================================================================
 
% weights (W): coupling between modes and hidden states (for 4 modes)
%--------------------------------------------------------------------------
W     = [1 0 1/2 1;
         1/2 1 1 1];
 
% Target spectrum - gamma, beta and alpha
%--------------------------------------------------------------------------
S     = -[32 32 48 128]' + 1i*2*pi*[38 16 12 4]';
S     = S(1:M.Nm);
W     = W(:,1:M.Nm);
W     = W/max(W(:));
Y.y   = [S; spm_vec(W)*32];
 
 
% Optimise using Gauss Newton
%==========================================================================
 
% create generative model of the eigenspectrum (M) and invert
%--------------------------------------------------------------------------
M.IS   = 'spm_ssm2s';
M.Nmax = 16;
pC     = V*diag(spm_vec(pC))*V;
M.pC   = pC;
for  i = 1:1
    [Ep Cp] = spm_nlsi_GN(M,[],Y);
    M.pE    = Ep;
    M.pC    = pC/2;
    M.x     = spm_dcm_neural_x(Ep,M);
end
 
 
% Show results with optimised parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Optimisation');
subplot(2,2,1)
sp      = spm_ssm2s(pE,M,[]);
sq      = spm_ssm2s(Ep,M,[]);
gp  = spm_s2csd(sp,w);
gq  = spm_s2csd(sq,w);
plot(w,gp,'r',w,gq,'b')
title('Before (red) and after (blue)','FontSize',16)
xlabel('Frequency')
ylabel('Spectral density of modes')
axis square
 
 
% show results with optimised parameters (measured)
%--------------------------------------------------------------------------
subplot(2,2,2), cla
for i = 1:numel(J)
    
    Ep.J  = sparse(1,J(i),1,1,nx);
    pE.J  = Ep.J;
    Gp    = spm_csd_mtf(pE,M,[]);
    Gq    = spm_csd_mtf(Ep,M,[]);
    
    plot(w,Gp{1},'r',w,Gq{1},'b'),hold on
    
end
title('Before (red) and after (blue)','FontSize',16)
xlabel('Frequency')
ylabel('Spectral responses')
axis square
hold off
 
 
% Show old and new priors
%--------------------------------------------------------------------------
subplot(4,1,3)
E    = spm_vec(Ep) - spm_vec(pE);
C    = diag(Cp);
j    = find(abs(E) > exp(-8) & C > 0);
spm_plot_ci(E(j,:),C(j,:))
 
title('Posterior updates','FontSize',16)
xlabel('Free parameters')
ylabel('Value')
set(gca,'XTick',1:length(j))
set(gca,'XTickLabel',spm_fieldindices(pE,j))
 
subplot(4,1,4)
dSdp = spm_diff('spm_ssm2s',Ep,M,[],1);
dSdp = diag(dSdp'*dSdp);
bar(log(dSdp(j,:)))
title('Overall contribution','FontSize',16)
xlabel('Free parameters')
ylabel('log precision')
set(gca,'XTick',1:length(j))
set(gca,'XLim',[0 (length(j) + 1)])
set(gca,'XTickLabel',spm_fieldindices(pE,j))
 
 
return
 
 
 
% Multi-start scheme
%==========================================================================
np    = length(spm_vec(M.pE));
ppC   = spm_sqrtm(diag(spm_vec(M.pC)));
for j = 1:8
    
    % use model inversion with a small number of iterations
    %----------------------------------------------------------------------
    M.Nmax = 16;
    M.nograph = 0;
    for i = 1:256
        M.P        = spm_vec(M.pE) + ppC*randn(np,1);
        [Ep,~,~,F] = spm_nlsi_GN(M,[],Y);
        FF(i)      = F;
        PP(:,i)    = spm_vec(Ep);
    end
    
    % get approximate Bayesian model average
    %----------------------------------------------------------------------
    [~, i] = max(FF);
    Ep     = spm_unvec(PP(:,i),M.pE);
    M.pE   = Ep;
    ppC    = ppC/2;
    
    % show progress
    %----------------------------------------------------------------------
    M.Nmax    = 32;
    M.nograph = 0;
    M.P       = Ep;
    Ep        = spm_nlsi_GN(M,[],Y);
    
    subplot(2,1,2);
    title(sprintf('Parameter updates: Multi-start - %i/8',j),'FontSize',16)
    drawnow
    
end
