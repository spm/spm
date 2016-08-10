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
% $Id: spm_induced_optimise_parameters.m 6856 2016-08-10 17:55:05Z karl $
 
% Parameters to optimise
%--------------------------------------------------------------------------
if ~nargin, PARAMS = {'G','T','L'}; end


% spectral specification
%==========================================================================
J      = [3 7];                % indices of hidden states producing outputs
Np     = length(J);

% Target spectrum - gamma, beta and alpha
%--------------------------------------------------------------------------
Hz     = 1:128;
s(:,1) = -[128 64 32 64]'   + 1j*2*pi*[4 12 48 64]';
s(:,2) = -[128 64 256 256]' + 1j*2*pi*[8 16 48 64]';



% Model specification
%==========================================================================
Nc    = 1;
Ns    = 1;
options.spatial  = 'LFP';
options.model    = 'TFM';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
M.J            = J;
M.Hz           = Hz;

% get priors
%--------------------------------------------------------------------------
[pE pC] = spm_dcm_neural_priors({0 0 0},{},1,options.model);
[pE pC] = spm_L_priors(M.dipfit,pE,pC);
[pE pC] = spm_ssr_priors(pE,pC);
[x,f]   = spm_dcm_x_neural(pE,options.model);

% suppress measurement noise
%--------------------------------------------------------------------------
pE.a(2) =      - 2;
pE.b    = pE.b - 16;
pE.c    = pE.c - 16;


% a target data
%--------------------------------------------------------------------------
Gu    = spm_csd_mtf_gu(pE,M.Hz);
for i = 1:Np
    csd(:,i,i) = 512*full(Gu.*sum(spm_s2csd(s(:,i),Hz),2));
end
 
% orders and model
%==========================================================================
nx      = length(spm_vec(x));
nu      = Ns;
u       = sparse(1,nu);
 
% fix priors if a subset of parameters are specified
%--------------------------------------------------------------------------
pV    = spm_vec(pC);
V     = pV - pV;
for i = 1:length(PARAMS)
    V(spm_fieldindices(pE,PARAMS{i})) = 1;
end
pC    = spm_unvec(V.*pV,pC);



% create LFP model
%--------------------------------------------------------------------------
M.IS   = 'spm_csd_mtf';
M.FS   = 'spm_diag_array';
M.g    = @(x,u,P,M) P.L*x(M.J(:));

M.f    = f;
M.x    = x;
M.n    = nx;
M.pE   = pE;
M.pC   = pC;
M.hE   = 8;
M.hC   = 1/128;
M.m    = nu;
M.l    = Np;
 
% solve for steady state
%--------------------------------------------------------------------------
M.x    = spm_dcm_neural_x(pE,M);
M.u    = u;
M.Nmax = 32;

% Optimisation: Target (Y)
%==========================================================================
[Ep Cp] = spm_nlsi_GN(M,[],{csd});


% Characterise contributions of parameters to spectral representation
%==========================================================================
 
% Show results with current (prior) parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Spectral responses'); clf
 
% and plot spectra
%--------------------------------------------------------------------------
subplot(2,2,1)
Gp    = spm_csd_mtf(pE,M);
Gq    = spm_csd_mtf(Ep,M);
plot(Hz,abs(spm_diag_array(Gp{1})),'r'), hold on
plot(Hz,abs(spm_diag_array(Gq{1})),'b'),  hold on
plot(Hz,abs(spm_diag_array(csd)),  '--'),  hold off

title({'Spectral responses';'Before (red) and after (blue)'},'FontSize',16)
xlabel('Frequency')
ylabel('Spectral density')
axis square

 
% Show results with optimised parameters
%--------------------------------------------------------------------------
subplot(2,2,2)
sp  = spm_ssm2s(pE,M);
sq  = spm_ssm2s(Ep,M);
gp  = spm_s2csd(sp,Hz);
gq  = spm_s2csd(sq,Hz);
plot(Hz,gp,'r',Hz,gq,'b')

title({'Eigenmodes';'Before (red) and after (blue)'},'FontSize',16)
xlabel('Frequency')
ylabel('Spectral density of modes')
axis square

 
% Show old and new priors
%==========================================================================

% change in parameters (and conditional confidence)
%--------------------------------------------------------------------------
E     = spm_vec(Ep) - spm_vec(pE);
C     = diag(Cp);

% eliminate an interesting parameters and sought on basis of contribution
%--------------------------------------------------------------------------
[c,j] = sort(C,'ascend');
j     = j(find(abs(E(j)) > exp(-8) & C(j) > 0));

subplot(4,1,3)
spm_plot_ci(E(j),C(j))
title('Posterior updates','FontSize',16)
xlabel('Free parameters')
ylabel('Value')
set(gca,'XTick',1:length(j))
set(gca,'XTickLabel',spm_fieldindices(pE,j))
 
subplot(4,1,4)
bar(-log(C(j,:)))
title('Contribution (log precision)','FontSize',16)
xlabel('Free parameters')
ylabel('log precision')
set(gca,'XTick',1:length(j))
set(gca,'XLim',[0 (length(j) + 1)])
set(gca,'XTickLabel',spm_fieldindices(pE,j))

return

% examine Jacobian
%==========================================================================
spm_figure('GetWin','Jacobian'); clf

% get transfer functions and Jacobian
%--------------------------------------------------------------------------
[S,K,s,w,t,dfdx] = spm_dcm_mtf(Ep,M,[]);
 
% and plot
%--------------------------------------------------------------------------
subplot(2,2,1)
imagesc(dfdx)
title('Jacobian','FontSize',16)
xlabel('hidden state')
ylabel('hidden state')
axis square


