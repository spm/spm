function DEM_demo_SOC
% Demo for a bird songs: this routine illustrates self organised
% criticality in terms of stimulus induced bifurcations and weak
% synchronisation of recognition (neuronal) dynamics. It uses the birdsong
% example, where stimuli are generated using a Lorentz attractor and
% modelled with the same attractor, with state dependent parameters.
% These control parameters are categorised in terms of a softmax function
% of point attractors in a (two-dimensional) perceptual space. We examine
% the self organised criticality in terms of Lyapunov exponents and the
% free energy - as a function of precision of the motion of hidden states
% (see code after return).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_SOC.m 4804 2012-07-26 13:14:18Z karl $
 
 
% hierarchical non-linear generative model of birdsong
%==========================================================================
 
% timing
%--------------------------------------------------------------------------
N        = 128;                      % length of stimulus (bins)
dt       = 1/64;                     % time bin (seconds)
t        = (1:N)*dt*1000;            % peristimulus time (ms)
 
% setup parameters for generalised filtering
%--------------------------------------------------------------------------
M(1).E.s = 1;
M(1).E.n = 4;
 
 
% level 1
%--------------------------------------------------------------------------
M(1).f  = '[-10 10 0; (v(1) - x(3)) -1 0; x(2) 0 -v(2)]*x/32;';
M(1).g  = 'x([2 3])';
M(1).x  = [0; 0; 1];
M(1).V  = exp(-2);
M(1).W  = exp(4);
 
 
% level 2
%--------------------------------------------------------------------------
M(2).f  = '1 - sum(exp(x)) - sparse(3,1)';
M(2).g  = '[1 15 28;0 8/3 8/3]*spm_softmax(x)';
M(2).x  = -log(3)*ones(3,1);
M(2).x  = [0; -4; -4];
M(2).v  = [0; 0];
M(2).V  = exp(-2);
M(2).W  = exp(0);
 
 
% generative process
%==========================================================================
G(1)    = M(1);
G(2).v  = [0; 0];
 
% create data (onset of a chaotic song)
%==========================================================================
S      = spm_phi(((1:N) - N/4)/(N/64));
U      = [1*(1 - S) + 28*S;0*(1 - S) + 8/3*S];
 
% create innovations & add causes
%--------------------------------------------------------------------------
DEM    = spm_DEM_generate(G,U,{[] []},{0 8},{8});
DEM.M  = M;
DEM.options.eigenvalues = 1;
 
% generalised Bayesian filtering and display
%==========================================================================
DEM    = spm_LAP(DEM);
 
% show song and perceptual inference
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)
 
subplot(3,2,5),colormap('pink')
spm_DEM_play_song(DEM.qU,N*dt);
axis square
title('percept','Fontsize',16)
 
subplot(3,2,6)
plot(t,spm_softmax(DEM.qU.x{2}))
axis square
title('perceptual inference','Fontsize',16)
xlabel('time (ms)')
ylabel('softmax probability')
spm_axis tight, grid on
 
 
% graphical characterisation of self-organised criticality
%==========================================================================
spm_figure('Getwin','Figure 1');
 
% get Lyapunov exponents
%--------------------------------------------------------------------------
k      = 8;
E      = sort(real(DEM.E),1,'descend');
CS     = sum(exp(k*(E - (~E)*1e8)));
 
% graphics
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(t,-DEM.FT,'k')
xlabel('time (ms)')
title('Free energy','fontsize',16)
axis square
 
subplot(2,2,2)
plot(t,(real(E)))
xlabel('time (ms)')
title('Lyapunov exponents','fontsize',16)
axis square
 
subplot(2,2,3)
plot(t,CS,'k')
xlabel('time (ms)')
title('Critical slowing','fontsize',16)
axis square
 
subplot(2,2,4)
plot(t,exp(real(E)))
xlabel('time (ms)')
title('Exponential','fontsize',16)
axis square


return

% repeat over different levels of sensory precision
%==========================================================================
 
% different levels of sensory precision
%--------------------------------------------------------------------------
V     = linspace(0,7,16);
F     = [];
E     = [];
B     = [];
LP    = [];
CI    = [];
for i = 1:length(V)
    
    % generalised filtering
    %----------------------------------------------------------------------
    DEM.M(1).W = exp(V(i));
    DEM        = spm_LAP(DEM);
    spm_DEM_qU(DEM.qU,DEM.pU)
    
    % accumulate statistics
    %----------------------------------------------------------------------
    F(i)       = -DEM.F;
    E          = real(DEM.E);
    E          = sort(E - (~E)*1e8,1,'descend');
    LP(:,i)    = E(1,:)';
    CI(i)      = mean(sum(exp(k*(E - (~E)*1e8))));
    
    % perceptual categorisation (after stimulus onset)
    %----------------------------------------------------------------------
    p          = mean(spm_softmax(DEM.qU.x{2}(:,N/2:N)),2);
    B(i)       = p(3);
    
end
 
 
% graphics showing critical range
%==========================================================================
spm_figure('Getwin','Figure 2'); clf
 
 
% get Lyapunov exponents and determine sensory threshold
%--------------------------------------------------------------------------
i    = N/2:N;
u    = V(find(B > 0.05, 1 ));
v    = V(find(B > 0.05, 1, 'last' ));

subplot(2,2,1)
plot(V,B,[u u],[0 1],'k--',[v v],[0 1],'k--')
xlabel('sensory log-precision')
ylabel('softmax probability')
title('Perceptual categorisation','fontsize',16)
axis square

subplot(2,2,2)
imagesc(V,t(i),LP(i,:)), hold on
plot([u u],[t(1) t(end)],'w--',[v v],[t(1) t(end)],'w--'), hold off
xlabel('sensory log-precision')
ylabel('peristimulus time (ms)')
title('Principal exponent','fontsize',16)
axis square
 
subplot(2,2,4)
plot(V,CI,[u u],[min(CI) max(CI)],'k--',[v v],[min(CI) max(CI)],'k--')
xlabel('sensory log-precision')
title('Critical slowing','fontsize',16)
ylabel('second moment')
axis square

subplot(2,2,3)
F   = F - min(F) + 128;
semilogy(V,F,[u u],[min(F) max(F)],'k--',[v v],[min(F) max(F)],'k--')
xlabel('sensory log-precision')
title('Free energy','fontsize',16)
ylabel('Free energy (nats)')
axis square, grid on, axis tight



return
 
 
% additional graphics for paper
%--------------------------------------------------------------------------
DEM.M(1).W   = exp(4);
DEM          = spm_LAP(DEM);
spm_DEM_qU(DEM.qU,DEM.pU)

% sonogram
%--------------------------------------------------------------------------
subplot(3,2,5),colormap('pink')
spm_DEM_play_song(DEM.qU,N*dt);
axis square
title('percept','Fontsize',16)

% perceptual inference
%--------------------------------------------------------------------------
subplot(3,2,6)
plot(t,spm_softmax(DEM.qU.x{2}))
axis square
title('perceptual inference','Fontsize',16)
xlabel('time (ms)')
ylabel('softmax probability')
spm_axis tight, grid on
 
