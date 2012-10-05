function Ep = spm_induced_optimise_parameters
% Demo routine that optimises free parameters
%==========================================================================
%
% This exemplar routine illustrates how one can adjust or tune prior
% parameter expectations to produce desired spectral responses as specified by
% the complex eigenvalue spectrum - or a reduced form that considered a
% small number of complex values and a single (unstable) real mode.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_induced_optimise_parameters.m 4989 2012-10-05 19:25:07Z karl $
 
 
% Model specification
%==========================================================================
 
% options
%--------------------------------------------------------------------------
Nc    = 1;
Ns    = 1;
options.spatial  = 'LFP';
options.model    = 'CMM';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
 
 
% get priors
%--------------------------------------------------------------------------
[pE pC] = spm_dcm_neural_priors({0 0 0},{},1,options.model);
[pE pC] = spm_L_priors(M.dipfit,pE,pC);
[pE pC] = spm_ssr_priors(pE,pC);
[x,f]   = spm_dcm_x_neural(pE,options.model);
  
% orders and model
%==========================================================================
nx      = length(spm_vec(x));
nu      = Ns;
u       = sparse(1,nu);
 
% create LFP model
%--------------------------------------------------------------------------
M.f    = f;
M.g    = 'spm_gx_erp';
M.x    = x;
M.n    = nx;
M.pE   = pE;
M.pC   = pC;
M.hE   = [8; 16];
M.hC   = diag(exp(-[16 16]));
M.m    = nu;
M.l    = Nc;
 
% solve for steady state
%--------------------------------------------------------------------------
M.x    = spm_dcm_neural_x(pE,M);
M.u    = u;
M.Hz   = 4:64;
M.Nmax = 32;
 
 
% Target (Y)
%==========================================================================
      
% coupling between modes and hidden states 
%--------------------------------------------------------------------------
M.out  = [0 0 0;
          1 0 2;
          0 0 0;
          0 2 0]*8;
 
 
% Target spectrum - gamma, beta and alpha
%--------------------------------------------------------------------------
S       = -[32 16 32]' + 1i*2*pi*[48 24 12]';
Y.y     = [S; spm_vec(M.out(any(M.out,2),:))];
Y.Q     = {diag([1 1 1 0 0 0 0 0 0]) diag([0 0 0 1 1 1 1 1 1])};
 
 
spm_figure('GetWin','Optimisation');
subplot(2,1,1)
[g,w]  = spm_s2csd(S);
plot(w,g)
title('Target spectral density','FontSize',16)
xlabel('Frequency')
ylabel('Power')
drawnow
 
% create generative model of the eigenspectrum (M) and invert
%--------------------------------------------------------------------------
M.IS  = 'spm_ssm2s';
for i = 1:8
    Ep   = spm_nlsi_GN(M,[],Y);
    M.pE = Ep;
end
 
% Show results with optimised parameters (pooled over eigenmodes)
%--------------------------------------------------------------------------
spm_figure('GetWin','Optimisation');
subplot(3,2,1)
sp      = spm_ssm2s(pE,M,[]);
sq      = spm_ssm2s(Ep,M,[]);
[gp,w]  = spm_s2csd(sp(1:3));
[gq,w]  = spm_s2csd(sq(1:3));
plot(w,g,'r',w,gq,'b')
title('target and predicted','FontSize',16)
xlabel('Frequency')
ylabel('spectral density of 3 principal modes')
legend({'target (red)','new prior (blue)'})
 
subplot(3,2,3)
plot(w,gp,'g')
title('Old prior','FontSize',16)
xlabel('Frequency')
ylabel('spectral density of 3 principal modes')
 
 
% Show results with optimised parameters (measured)
%--------------------------------------------------------------------------
J     = find(any(M.out,2));
subplot(3,2,2), cla
subplot(3,2,4), cla
for i = 1:numel(J)
    
    
    Ep.J       = M.pE.J - M.pE.J;
    Ep.J(J(i)) = 1;
    pE.J       = Ep.J;
    
    [Gp w]     = spm_csd_mtf(pE,M);    
    [Gq w]     = spm_csd_mtf(Ep,M);
    
    subplot(3,2,2)
    plot(w,Gq{1}),hold on
    title('spectral responses (new)','FontSize',16)
    xlabel('Frequency')
    ylabel('spectral density')
    
    subplot(3,2,4)
    plot(w,Gp{1}),hold on
    title('spectral responses (old)','FontSize',16)
    xlabel('Frequency')
    ylabel('spectral density')
    
end
hold off
 
 
% Show old and new priors
%--------------------------------------------------------------------------
subplot(3,1,3)
E    = [spm_vec(pE) spm_vec(Ep)];
C    = [spm_vec(pC) spm_vec(pC)];
j    = find(abs(diff(E,[],2)) > exp(-16));
spm_plot_ci(E(j,:),C(j,:))
 
set(gca,'YLim',[min(min(E(j,:))) - 1,max(max(E(j,:))) + 1]);
title('old and new priors','FontSize',16)
xlabel('free parameters')
ylabel('value')
legend({'old','new'})
 
set(gca,'XTickLabel',spm_fieldindices(pE,j))
