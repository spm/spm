function spm_mfa_demo
% Demonstration of mean field approximation for spiking neurons.  This demo
% is just meant to illustrate how one gets from the differential equations
% of a Hodgkin Huxley like neuron to ensemble dynamics through a Fokker
% Planck (ensemble density) formulation.  The key to doing this rests on
% the use of time since last spike as a hidden state (and support of the
% ensemble density).  This means the ensemble dynamics can be expressed as
% modes over time, which effectively converts a spiking model into a rate
% model.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mfa_demo.m 4936 2012-09-18 19:47:55Z karl $
 
% model specification (Parameters and initial states)
%--------------------------------------------------------------------------
P     = [(1 - 1/8) 1/16 1/16]';
x0    = [0.01; 0.01; 0.01; 0.01; -80; 1e-2];
 
% model specification (equations of motion)
%--------------------------------------------------------------------------
M.f   = 'spm_fx_hh';
M.g   = 'spm_gx_hh';
M.pE  = P;
M.x   = x0;
M.m   = 1;
M.n   = length(x0);
M.l   = size(feval(M.g,M.x,M.pE),1);
 
% model dynamics
%==========================================================================
 
% create inputs
%--------------------------------------------------------------------------
dt    = 1e-4;
n     = 100e-3/dt;
u     = 8;
t     = (1:n)*dt*1e3;
 
% integrate with input - output = E{x} and create response
%--------------------------------------------------------------------------
for j = 1:32
    x     = x0;
    for i = 1:n
        dx      = spm_fx_hh(x,u,P) + randn(size(x0));
        x       = x + dx*dt;
        V(i,j)  = x(5);
        T(i,j)  = x(6);
    end
end
 
% single trajectories over time
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
subplot(2,2,1)
plot(t,V)
title('Trajectories over time','FontSize',16)
xlabel('time (ms)','FontSize',12)
ylabel('depolarization','FontSize',12)
axis square
grid on
 
 
% sand in phase (state) space
%--------------------------------------------------------------------------
clear mov
subplot(2,2,2)
V = V(1:2:end,:);
T = T(1:2:end,:);
for i = 1:length(V)
    
    % plot
    %----------------------------------------------------------------------
    plot(T(i,:)*1000,V(i,:),'o','MarkerSize',8)
    axis([0 32 -100 -40])
    title({'Trajectories in state-space';'Click axis for movie'},'FontSize',16)
    xlabel('peri-spike time (ms)','FontSize',12)
    ylabel('depolarization','FontSize',12)
    axis square
    grid on
    drawnow
    
    % addframe
    %----------------------------------------------------------------------
    mov(i) = getframe(gca);
    
end
 
% set ButtonDownFcn
%--------------------------------------------------------------------------
set(gca,'Userdata',{mov,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
 
 
 
% ensemble dynamics (Fokker Planck formulation)
%==========================================================================
 
% event space
%--------------------------------------------------------------------------
x     = {linspace(0,1,2),...
         linspace(0,1,2),...
         linspace(0,1,4),...
         linspace(0,1,2),...
         linspace(-92,-48,16),...
         linspace(0,1e-1,8)};
dx    = zeros(1,length(x0));
nx    = zeros(1,length(x0));
for i = 1:length(x)
    dx(i) = x{i}(2) - x{i}(1);
    nx(i) = length(x{i});
end
M.W   = diag([1/8 1/8 1/8 1/8 4 1e-6].^2);
 
% create ensemble
%--------------------------------------------------------------------------
G     = spm_mfa_G(M,x);
 
% iterate a single ensemble
%--------------------------------------------------------------------------
[M0,M1,L1] = spm_mfa_bi_multi(G);
L          = [G.p0 G.u];
 
dt    = 1e-2;
eJ    = expm(full(M0)*dt);                   % no input
eU    = expm(full(M0 + M1{1}*2)*dt);         % input
 
p     = sparse(1,1,1,length(M0),1);
y     = zeros(3,64);
 
subplot(2,2,4)
clear mov
for i = 1:64
    
    % input-dependent Jacobian
    %----------------------------------------------------------------------
    if i > 16 && i < 32
        J  = eU;
    else
        J  = eJ;
    end
    
    % Gibb's density
    %----------------------------------------------------------------------
    d = real(L*p);
    p = J*p;
    y(:,i) = real(L1*p);
    d = reshape(d,nx(1),nx(2),nx(3),nx(4),nx(5),nx(6));
    d = squeeze(sum(sum(sum(sum(d,1),2),3),4));
    
    imagesc(G.x{6}*1000,G.x{5},spm_interp((1 - d),8))
    title({'Trajectories in state-space';'Click axis for movie'},'FontSize',16)
    xlabel('peri-spike time (ms)','FontSize',12)
    ylabel('depolarization','FontSize',12)
    axis xy square
    grid on
    drawnow
    
    mov(i)  = getframe(gca);
end
 
% set ButtonDownFcn
%--------------------------------------------------------------------------
h = findobj(gca,'type','image');
set(h,'Userdata',{mov,16})
set(h,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
 
% Expected states
%--------------------------------------------------------------------------
% y(1) = V transmembrane potential mV (c.f. LFP)
% y(2) = spike rate (Hz) = 1/PST
% y(3) = dendritic energy = g(1)*x(1).*(V(1) - v).^2 + ... (mV.^2mS)
%--------------------------------------------------------------------------
subplot(2,2,3)
plot((1:length(y))*dt*1000,y(:,:))
title({'Expected satets'},'FontSize',16)
xlabel('time (ms)','FontSize',12)
ylabel('expectded state','FontSize',12)
grid on
axis square
legend({'potential (mV)','spike rate (Hz)','dendritic energy (mV^2mS)'})
 
 
 
% Ensemble dynamics of couples populations
%==========================================================================
 
% define coupling and find bilinear reduction
%--------------------------------------------------------------------------
S(1)   = G;                             % excitatory (1)
S(2)   = G;                             % inhibitory (1)
S(3)   = G;                             % excitatory (2)
S(4)   = G;                             % inhibitory (2)
 
C{2,1} = [0 1 0; 0 0 0; 0 0 0];         % AMPA 1 -> 2
C{1,2} = [0 0 0; 0 1 0; 0 0 0];         % GABA 2 -> 1
C{4,3} = [0 1 0; 0 0 0; 0 0 0];         % AMPA 3 -> 4
C{3,4} = [0 0 0; 0 1 0; 0 0 0];         % GABA 4 -> 3
C{3,1} = [0 1 0; 0 0 0; 0 0 0];         % AMPA 1 -> 3
C{1,3} = [0 0 0; 0 0 0; 0 1 0];         % NMDA 3 -> 1
 
 
% ensemble in bilinear form
%--------------------------------------------------------------------------
[M0,M1,L,M2] = spm_mfa_bi_multi(S,C);
 
% coupling parameters
%--------------------------------------------------------------------------
Pu     = [1 0 0 0];
Pc     = [0 1 1 0;
          1 0 0 0;
          1 0 0 1
          0 0 1 0]/16;
P      = [Pu(:); Pc(:)];
 
% bilinear operators
%--------------------------------------------------------------------------
MB.M0      = M0;
MB.M1      = M1;
MB.M2      = M2;
MB.L       = L(1:3:end,:);                   % only consider LFP outputs
[M0,M1,L1] = spm_mfa_bi(MB,P);
 
 
% solve for (effectively) a rate model using bilinear operators
%--------------------------------------------------------------------------
eJ    = expm(full(M0)*dt);                   % no input
eU    = expm(full(M0 + M1{1}*2)*dt);         % input
p     = sparse(1,1,1,length(M0),1);
y     = sparse(4,64);
for i = 1:64
    
    % input-dependant Jacobian (Fokker Planck operator)
    %----------------------------------------------------------------------
    if i > 16 && i < 32
        J  = eU;
    else
        J  = eJ;
    end
    
    % Gibb's density
    %----------------------------------------------------------------------
    p = J*p;
    y(:,i) = real(L1*p);
    
end
 
% Expected states
%--------------------------------------------------------------------------
% y(1) = V transmembrane potential mV (c.f. LFP)
% y(2) = spike rate (Hz) = 1/PST
% y(3) = dendritic energy = g(1)*x(1).*(V(1) - v).^2 + ... (mV.^2mS)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
 
subplot(2,1,1)
plot((1:length(y))*dt*1000,y)
title({'Expected states for (4) coupled populations';'(simulated ERP)'},'FontSize',16)
xlabel('time (ms)','FontSize',12)
ylabel('depolarization (mV)','FontSize',12)
grid on
axis square
