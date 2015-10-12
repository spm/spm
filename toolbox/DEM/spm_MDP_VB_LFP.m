function u = spm_MDP_VB_LFP(MDP,UNITS)
% auxiliary routine for plotting simulated electrophysiological responses
% FORMAT u = spm_MDP_VB_LFP(MDP,UNITS)
%
% u - selected unit responses {number of trials, number of units}
%
% MDP - structure (see spm_MDP_VB
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_VB_LFP.m 6567 2015-10-12 09:29:04Z karl $
 
 
% deal with a sequence of trials
%==========================================================================
 
% dimensions
%--------------------------------------------------------------------------
Nt    = length(MDP);               % number of trials
NT    = size(MDP(1).V,1) + 1;      % number of transitions
Nx    = size(MDP(1).A,2);          % number of states
Nb    = size(MDP(1).xn,1);         % number of time bins per transition
 
% units to plot
%--------------------------------------------------------------------------
ALL   = [];
for i = 1:NT
    for j = 1:Nx
        ALL(:,end + 1) = [j;i];
    end
end
if nargin < 2;
    UNITS = ALL;
end
    
% summary statistics
%==========================================================================
for i = 1:Nt
    
    % all units
    %----------------------------------------------------------------------
    for j = 1:size(ALL,2)
        for k = 1:NT
            xj{k,j} = MDP(i).xn(:,ALL(1,j),ALL(2,j),k);
        end
    end
    x{i,1} = xj;
    
    % selected units
    %----------------------------------------------------------------------
    for j = 1:size(UNITS,2)
        for k = 1:NT
            uj{k,j} = MDP(i).xn(:,UNITS(1,j),UNITS(2,j),k);
        end
    end
    u{i,1} = uj;
    
    % dopamine or changes in precision
    %----------------------------------------------------------------------
    dn(:,i) = MDP(i).dn;
end

if nargout, return, end
 
% phase amplitude coupling
%==========================================================================
dt  = 1/64;                              % time bin (seconds)
t   = (1:(Nb*NT*Nt))*dt;                 % time (seconds)
Hz  = 4:32;                              % frequency range
n   = 1/(4*dt);                          % window length
w   = Hz*(dt*n);                         % cycles per window
 
% simulated local field potential
%--------------------------------------------------------------------------
LFP = spm_cat(x);
i   = UNITS(1,end) + (UNITS(2,end) - 1)*Nx;
 
if Nt == 1, subplot(3,2,1), else subplot(4,1,1),end
imagesc(t,1:(Nx*NT),LFP'),title('Unit responses','FontSize',16)
xlabel('time (seconds)','FontSize',12), ylabel('unit','FontSize',12)
grid on, set(gca,'XTick',(1:(NT*Nt))*Nb*dt)
grid on, set(gca,'YTick',(1:NT)*Nx)
if NT*Nt > 32, set(gca,'XTickLabel',[]), end
if Nt == 1,    axis square,              end
 
% time frequency analysis and theta phase
%--------------------------------------------------------------------------
wft = spm_wft(LFP,w,n);
csd = sum(abs(wft),3);
lfp = LFP(:,i);
phi = spm_iwft(wft(1,:,i),w(1),n);
lfp = 4*lfp/std(lfp) + 16;
phi = 4*phi/std(phi) + 16;
 
if Nt == 1, subplot(3,2,3), else subplot(4,1,2),end
imagesc(t,Hz,csd), axis xy, hold on
plot(t,lfp,'w:',t,phi,'w'), hold off
grid on, set(gca,'XTick',(1:(NT*Nt))*Nb*dt)

title('Time-frequency response','FontSize',16)
xlabel('time (seconds)','FontSize',12), ylabel('frequency','FontSize',12)
if Nt == 1, axis square, end
 
% local field potentials
%==========================================================================
if Nt == 1, subplot(3,2,2), else subplot(4,1,3),end
plot(t,spm_cat(u)),     hold on, spm_axis tight, a = axis;
plot(t,spm_cat(x),':'), hold off
grid on, set(gca,'XTick',(1:(NT*Nt))*Nb*dt), axis(a)
for i = 2:2:Nt
    h = patch(((i - 1) + [0 0 1 1])*NT*Nb*dt,a([3,4,4,3]),-[1 1 1 1],'w');
    set(h,'LineStyle',':','FaceColor',[1 1 1] - 1/32);
end
title('Local field potentials','FontSize',16)
xlabel('time (seconds)','FontSize',12)
ylabel('Response','FontSize',12)
if Nt == 1, axis square, end

% firing rates
%==========================================================================
qu   = cumsum(spm_cat(u)) + 1/Nx;
qx   = cumsum(spm_cat(x)) + 1/Nx;
if Nt == 1, subplot(3,2,4)
    plot(t,qu),     hold on, spm_axis tight, a = axis;
    plot(t,qx,':'), hold off
    grid on, set(gca,'XTick',(1:(NT*Nt))*Nb*dt), axis(a)
    title('Firing rates','FontSize',16)
    xlabel('time (seconds)','FontSize',12)
    ylabel('Response','FontSize',12)
    axis square
end

% simulated dopamine responses
%==========================================================================
if Nt == 1, subplot(3,1,3), else subplot(4,1,4),end
bar(spm_vec(dn),1,'k'), title('Phasic dopamine responses','FontSize',16)
xlabel('time (updates)','FontSize',12)
ylabel('change in precision','FontSize',12), spm_axis tight
if Nt == 1, axis square, end
 
