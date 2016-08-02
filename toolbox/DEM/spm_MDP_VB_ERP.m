function [x] = spm_MDP_VB_ERP(MDP,FACTOR,T)
% auxiliary routine for plotting hierarchical electrophysiological responses
% FORMAT [x] = spm_MDP_VB_ERP(MDP,FACTOR,T)
%
% u - unit rate of change of firing (simulated voltage)
% v - unit responses
%
% T - flag to return cell of expectations (at time T)
%
% MDP - structure (see spm_MDP_VB)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_ERP.m 6853 2016-08-02 08:29:03Z karl $


% defaults
%==========================================================================
try, f = FACTOR; catch, f = 1;      end
xn  = MDP.xn{f};

% dimensions
%--------------------------------------------------------------------------
Nb  = size(xn,1);         % number of time bins per epochs
Nx  = size(xn,2);         % number of states
Ne  = size(xn,3);         % number of epochs

% units to plot
%--------------------------------------------------------------------------
UNITS = [];
if nargin < 3
    for i = 1:Ne
        for j = 1:Nx
            UNITS(:,end + 1) = [j;i];
        end
    end
else
    for j = 1:Nx
        UNITS(:,end + 1) = [j;T];
    end
end

% expected states
%==========================================================================
for k = 1:Ne
    for j = 1:size(UNITS,2)
        x{k,j} = xn(:,UNITS(1,j),UNITS(2,j),k);
    end
    if isfield(MDP,'mdp')
        y{k} = spm_MDP_VB_ERP(MDP.mdp(k),1,1);
    else
        y{k} = [];
    end
end

if nargin > 2, return, end

% synchronise responses
%--------------------------------------------------------------------------
u   = {};
v   = {};
uu  = spm_cat(x(1,:));
for k = 1:Ne
    
    % low-level
    %----------------------------------------------------------------------
    v{end + 1,1}     = spm_cat(y{k});
    if k > 1
        u{end + 1,1} = ones(size(v{end,:},1),1)*u{end,1}(end,:);
    else
        u{end + 1,1} = ones(size(v{end,:},1),1)*uu(1,:);
    end
    
    % high-level
    %----------------------------------------------------------------------
    u{end + 1,1} = spm_cat(x(k,:));
    v{end + 1,1} = ones(size(u{end,:},1),1)*v{end,1}(end,:);
    
end

u  = spm_cat(u);
v  = spm_cat(v);


if nargout, return, end

% phase amplitude coupling
%==========================================================================
dt  = 1/64;                              % time bin (seconds)
t   = (1:(Nb*Ne*Nt))*dt;                 % time (seconds)
Hz  = 4:32;                              % frequency range
n   = 1/(4*dt);                          % window length
w   = Hz*(dt*n);                         % cycles per window

% simulated local field potential
%--------------------------------------------------------------------------
LFP = spm_cat(x);

if Nt == 1, subplot(3,2,1), else subplot(4,1,1),end
imagesc(t,1:(Nx*Ne),spm_cat(z)'),title('Unit responses','FontSize',16)
xlabel('time (seconds)','FontSize',12), ylabel('unit','FontSize',12)
grid on, set(gca,'XTick',(1:(Ne*Nt))*Nb*dt)
grid on, set(gca,'YTick',(1:Ne)*Nx)
if Ne*Nt > 32, set(gca,'XTickLabel',[]), end
if Nt == 1,    axis square,              end

% time frequency analysis and theta phase
%--------------------------------------------------------------------------
wft = spm_wft(LFP,w,n);
csd = sum(abs(wft),3);
lfp = sum(LFP,2);
phi = spm_iwft(sum(wft(1,:,:),3),w(1),n);
lfp = 4*lfp/std(lfp) + 16;
phi = 4*phi/std(phi) + 16;

if Nt == 1, subplot(3,2,3), else subplot(4,1,2),end
imagesc(t,Hz,csd), axis xy, hold on
plot(t,lfp,'w:',t,phi,'w'), hold off
grid on, set(gca,'XTick',(1:(Ne*Nt))*Nb*dt)

title('Time-frequency response','FontSize',16)
xlabel('time (seconds)','FontSize',12), ylabel('frequency','FontSize',12)
if Nt == 1, axis square, end

% local field potentials
%==========================================================================
if Nt == 1, subplot(3,2,4), else subplot(4,1,3),end
plot(t,spm_cat(u)),     hold off, spm_axis tight, a = axis;
plot(t,spm_cat(x),':'), hold on
plot(t,spm_cat(u)),     hold off, axis(a)
grid on, set(gca,'XTick',(1:(Ne*Nt))*Nb*dt),
for i = 2:2:Nt
    h = patch(((i - 1) + [0 0 1 1])*Ne*Nb*dt,a([3,4,4,3]),-[1 1 1 1],'w');
    set(h,'LineStyle',':','FaceColor',[1 1 1] - 1/32);
end
title('Local field potentials','FontSize',16)
xlabel('time (seconds)','FontSize',12)
ylabel('Response','FontSize',12)
if Nt == 1, axis square, end

% firing rates
%==========================================================================
qu   = spm_cat(v);
qx   = spm_cat(z);
if Nt == 1, subplot(3,2,2)
    plot(t,qu),     hold on, spm_axis tight, a = axis;
    plot(t,qx,':'), hold off
    grid on, set(gca,'XTick',(1:(Ne*Nt))*Nb*dt), axis(a)
    title('Firing rates','FontSize',16)
    xlabel('time (seconds)','FontSize',12)
    ylabel('Response','FontSize',12)
    axis square
end


