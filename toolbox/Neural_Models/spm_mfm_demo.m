% Demo routine for mean-field models
%==========================================================================
% 
% This demo compares and contrasts neural-mass and mean-field models of a 
% single population, using the model described in Marreiros et al 2008). 
% We start by comparing the impulse response of an small ensemble that has 
% some (but not marked) finite size effects) with that of the mean-field
% and neural-mass approximations.  The key difference between these models
% is that the means-field has states that describe the change in dispersion 
% or covariance among the first-order states (current and voltage).
% 
% We then move on to comparing responses to inputs that are transient and 
% sustained, to show that the mean-field model retains key nonlinearities 
% and can show [plausible] bifurcations, as sustained input levels are
% increased.  This is characterised using Fourier transforms, which are %
% plotted alongside spiking responses.  See Marreiros et al:
% 
% Population dynamics under the Laplace assumption.
% 
% A Marreiros, J Daunizeau, S Kiebel, L Harrison & Karl Friston
% 
% Abstract
% In this paper, we describe a generic approach to modelling dynamics in
% neuronal populations.  This approach retains a full density on the states
% of neuronal populations but resolves the problem of solving
% high-dimensional problems by re-formulating density dynamics in terms of
% ordinary differential equations on the sufficient statistics of the
% densities considered.  The particular form for the population density we
% adopt is a Gaussian density (c.f., a Laplace assumption). This means
% population dynamics are described completely by equations governing the
% evolution of the population’s mean and covariance.  We derive these
% equations from the Fokker-Planck formalism and illustrate their
% application to a reasonably simple conductance-based model of neuronal
% exchanges.  One interesting aspect of this formulation is that we can
% uncouple the mean and covariance to furnish a neural-mass model, which
% rests only on the populations mean.  This enables to compare equivalent
% mean-field and neural-mass models of the same populations and evaluate,
% quantitatively, the contribution of population variance to the expected
% dynamics.  The mean-field model presented here will form the basis of a
% dynamic causal model of observed electromagnetic signals in future work.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mfm_demo.m 1277 2008-03-28 18:36:49Z karl $
 
 
% number of regions in coupled map lattice
%--------------------------------------------------------------------------
clear
n     = 1;
 
% extrinsic network connections
%--------------------------------------------------------------------------
if n > 1
A{1}  = diag(ones(n - 1,1),-1);
else
    A{1} = 0;
end
A{2}  = A{1}';
A{3}  = sparse(n,n);
B     = {};
C     = sparse(1,1,1,n,1);
 
 
% get connectivity and other priors
%--------------------------------------------------------------------------
[pE,pC] = spm_nmm_priors(A,B,C);
 
% initialise states
%--------------------------------------------------------------------------
[x M] = spm_x_mfm(pE);
 
 
% create exogenous inputs
%==========================================================================
dt    = 2;
t     = [1:dt:256]';
N     = length(t);
U.dt  = dt/1000;
U.u   = 128*(exp(-(t - 64).^2/8) + rand(N,1)/exp(32));
 
 
% Integrate system to see Transient response - ensemble
%==========================================================================
np    = 8;                                % number of particles in ensemble
 
% initialise ensemble
%--------------------------------------------------------------------------
for i = 1:np
    X(i,:,:,:) = x{1};
end
 
% diffusion
%--------------------------------------------------------------------------
D    = diag([1/32 1 1]);
dfdw = kron(spm_sqrtm(D/2),speye(np*n*3));


 
% integrate ensemble Y(time, particle, source, population, state)
%==========================================================================
fprintf('integrating ensemble\n')

% burn in
%--------------------------------------------------------------------------
for i = 1:16
    [dfdx f]     = spm_diff('spm_fx_mfm_ensemble',X,U.u(1),pE,1);
    dX           = spm_sde_dx(dfdx,dfdw,f,U.dt);
    X            = spm_unvec(dX + spm_vec(X),X);
    fprintf('time-step (burn in) %i (%i)\n',i,N)
end
 
% integrate over peristimulus time
%--------------------------------------------------------------------------
for i = 1:N
 
    [dfdx f]     = spm_diff('spm_fx_mfm_ensemble',X,U.u(i),pE,1);
    dX           = spm_sde_dx(dfdx,dfdw,f,U.dt);
    X            = spm_unvec(dX + spm_vec(X),X);
    Y(i,:,:,:,:) = X;
    
    fprintf('time-step %i (%i)\n',i,N)
end
 
% input
%--------------------------------------------------------------------------
pop   = {'stellate', 'interneurons', 'pyramidal'};
state = {'Voltage', 'excitatory conductance', 'inhibitory conductance'};

p  = 3;                                    % default population (pyramidal)
T  = t;                                    % perstimulus time (ms)
c  = [1 1 1]*.9;                           % default colour for plots
subplot(3,1,1)
plot(T,U.u)
title('input')
axis square
xlabel('time (ms)')
 
% plot
%--------------------------------------------------------------------------
subplot(3,1,2)
plot(t,squeeze(Y(:,:,1,p,1)),':')
title('Pyramidal depolarisation')
axis square
xlabel('time (ms)')
 
% plot moments
%--------------------------------------------------------------------------
m   = mean(Y(:,:,1,p,1),2)';
s   = std(Y(:,:,1,p,1),0,2)'*1.64;
subplot(3,1,3)
plot(t,m*0,'k:'),hold on
patch([t(:)' fliplr(t(:)')],[(m(:)' - s(:)') fliplr(m(:)' + s(:)')],c)
plot(t,m,'r','LineWidth',2), hold off
title('Pyramidal depolarisation')
axis square
xlabel('time (ms)')
drawnow
  
spm_demo_proceed

 
% plot
%--------------------------------------------------------------------------
for i = 1:3
    for j = 1:3
        subplot(3,3,(i - 1)*3 + j)
        plot(t,squeeze(Y(:,:,1,i,j)),'-','color',[1 1 1]*.7)
        title([pop{i} ': ' state{j}])
        xlabel('time (ms)')
    end
end
drawnow

spm_demo_proceed

 
 
% Integrate system to see Transient response - MFM
%==========================================================================
NM    = M;                              % neural mass model
NM.x  = M.x{1};                         % remove second-order moments
MFM   = spm_int_B(pE,M, U);
NMM   = spm_int_B(pE,NM,U);

% LFP - ensemble
%--------------------------------------------------------------------------
m   = mean(Y(:,:,1,p,1),2)';
s   =  std(Y(:,:,1,p,1),0,2)'*1.64;
subplot(3,1,1)
plot(t,m*0,'k:'),hold on
patch([t(:)' fliplr(t(:)')],[(m(:)' - s(:)') fliplr(m(:)' + s(:)')],c)
plot(T,m,'r','LineWidth',2), hold off
title([pop{p} ': Ensemble'])
axis square
xlabel('time (ms)')


% LFP - MFM
%--------------------------------------------------------------------------
subplot(3,1,2)
m   = MFM(:,p)';
s   = sqrt(MFM(:,p*9 + 1)')*1.64;
plot(t,m*0,'k:'),hold on
patch([t(:)' fliplr(t(:)')],[(m(:)' - s(:)') fliplr(m(:)' + s(:)')],c)
plot(t,m,'r','LineWidth',2), hold off
title([pop{p} ': MFM'])
axis square
xlabel('time (ms)')

 
% LFP - NMM
%--------------------------------------------------------------------------
subplot(3,1,3)
m   = NMM(:,p)';
s   = sqrt(64)*1.64;
plot(t,m*0,'k:'),hold on
patch([t(:)' fliplr(t(:)')],[(m(:)' - s(:)') fliplr(m(:)' + s(:)')],c)
plot(t,m,'r','LineWidth',2), hold off
title([pop{p} ': NNM'])
axis square
xlabel('time (ms)')
drawnow

spm_demo_proceed


% create exogenous inputs for responses to transient and sustained input
%==========================================================================
dt    = 2;
t     = [1:dt:512]';
N     = length(t);
U.dt  = dt/1000;
 
 
% responses to different inputs - spikes
%==========================================================================
clear YMF YNM
p     = 3;                                               % pyramidal cells
u     = [1:32];                                          % peak input
for i = 1:length(u)
    
    % create exogenous inputs
    %----------------------------------------------------------------------
    U.u   = u(i)*(exp(-(t - 64).^2/8) + randn(N,1)*exp(-16));
 
    % Integrate systems
    %----------------------------------------------------------------------
    MFM   = spm_int_B(pE,M, U);
    NMM   = spm_int_B(pE,NM,U);
 
    YMF(:,:,i) = MFM(:,1:3);
    YNM(:,:,i) = NMM(:,1:3);
    
    fprintf('input level (spike) %i (%i)\n',i,length(u))
 
end
 
subplot(2,2,1)
imagesc(u,t,squeeze(YMF(:,p,:)))
xlabel('input amplitude (spike)')
ylabel('time (ms)')
title([pop{p} ': MFM'])
 
subplot(2,2,2)
imagesc(u,t,squeeze(YNM(:,p,:)))
xlabel('input amplitude (spike)')
ylabel('time (ms)')
title([pop{p} ': NMM'])
drawnow
 
 
% responses to different inputs - sustained
%--------------------------------------------------------------------------
clear YMF YNM
for i = 1:length(u)
    
    % create exogenous inputs (after 64 ms)
    %----------------------------------------------------------------------
    U.u   = u(i)*(t > 64);
 
    % Integrate systems
    %----------------------------------------------------------------------
    MFM   = spm_int_B(pE,M, U);
    NMM   = spm_int_B(pE,NM,U);
 
    YMF(:,:,i) = MFM(:,1:3);
    YNM(:,:,i) = NMM(:,1:3);
    
    fprintf('input level (sustained) %i (%i)\n',i,length(u))
 
end
 
subplot(2,2,3)
imagesc(u,t,squeeze(YMF(:,p,:)))
xlabel('input amplitude (sustained)')
ylabel('time (ms)')
title([pop{p} ': MFM'])
 
subplot(2,2,4)
imagesc(u,t,squeeze(YNM(:,p,:)))
xlabel('input amplitude (sustained)')
ylabel('time (ms)')
title([pop{p} ': NMM'])
drawnow
 
spm_demo_proceed

 
% frequency responses
%--------------------------------------------------------------------------
clear FMF SMF
p     = 3;
for i = 1:length(u)
    
    % Fourier transform
    %----------------------------------------------------------------------
    FMF(:,i) = abs(fft(YMF(:,p,i)));
    
    % mean spiking
    %----------------------------------------------------------------------
    for j = 1:N
        m        = spm_gx_mfm(YMF(j,:,i),0,pE,M);
        SMF(j,i) = m(p);
    end
    
end
 
 
Hz  = [0:N]/(N*U.dt);
i   = find(Hz > 2 & Hz < 64);
Hz  = Hz(i);
FMF = FMF(i,:);
 
subplot(2,2,1)
imagesc(u,Hz,FMF)
xlabel('input amplitude')
ylabel('Frequency (Hz)')
title([pop{p} ': MFM'])
 
subplot(2,2,2)
plot(Hz,FMF(:,[8 16]))
xlabel('Frequency (Hz)')
ylabel('frequency response')
title([pop{p} ': MFM'])
legend({sprintf('input amplitude %i',u(16)) sprintf('input amplitude %i',u(32))})
 
subplot(2,2,3)
plot(u,mean(SMF),'k')
xlabel('input amplitude')
ylabel('firing response')
title([pop{p} ': MFM'])
axis([u(1) u(end) 0 0.5])
 
subplot(2,2,4)
imagesc(u,t,SMF)
xlabel('input amplitude')
ylabel('time (ms)')
title([pop{p} ': MFM'])
drawnow
 
return
 
 

% NOTES:
% time-frequency
%==========================================================================
W     = 256;
w     = 4:1/4:32;
cpW   = w*W*U.dt;
subplot(2,2,3)
imagesc(t*1000,w,abs(spm_wft(LFP(:,7),cpW,W)).^2);
title('time-frequency response')
axis square xy
xlabel('time (ms)')
 
 
% Kernels
%==========================================================================
 
 
% augment and bi-linearise
%--------------------------------------------------------------------------
[M0,M1,L1,L2] = spm_bireduce(M,M.pE);
 
% compute kernels (over 64 ms)
%--------------------------------------------------------------------------
N          = 64;
dt         = 3/1000;
t          = [1:N]*dt*1000;
[K0,K1,K2] = spm_kernels(M0,M1,L1,L2,N,dt);
 
subplot(2,1,1)
plot(t,K1)
title('1st-order Volterra kernel')
axis square
xlabel('time (ms)')
 
subplot(2,1,2)
imagesc(t,t,K2(1:64,1:64,1))
title('2nd-order Volterra kernel')
axis square
xlabel('time (ms)')
