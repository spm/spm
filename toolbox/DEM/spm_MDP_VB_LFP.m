function [MDP] = spm_MDP_VB_LFP(MDP,UNITS)
% auxiliary routine for plotting simulated electrophysiological responses
% FORMAT [MDP] = spm_MDP_VB_LFP(MDP,UNITS)
%
% MDP - structure (see spm_MDP_VB
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_LFP.m 6517 2015-08-10 11:21:53Z karl $


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
    
    % ddopamine or changes in precision
    %----------------------------------------------------------------------
    da(:,i) = MDP(i).da;
end

% phase amplitude coupling
%==========================================================================
dt  = 1/64;                              % time bin (seconds)
t   = (1:(Nb*NT*Nt))*dt;                    % ime (seconds)
Hz  = 4:32;                              % frequency range
n   = 1/(2*dt);                          % window length
w   = Hz*(dt*n);                         % cycles per window
K   = exp(-(Hz - 4).^2/4);               % filter (theta)

% simulated local field potential
%--------------------------------------------------------------------------
LFP = spm_cat(x);

if Nt == 1, subplot(2,2,1), else subplot(4,1,1),end
imagesc(t,1:(Nx*NT),LFP'),title('unit responses','FontSize',16)
xlabel('time (seconds)','FontSize',12), ylabel('unit','FontSize',12)
grid on, set(gca,'XTick',(1:(NT*Nt))*Nb*dt)
grid on, set(gca,'YTick',(1:NT)*Nx)
if NT*Nt > 32, set(gca,'XTickLabel',[]), end
if Nt == 1,    axis square,              end

% ttime frequency analysis and theta phase
%--------------------------------------------------------------------------
wft = spm_wft(LFP,w,n);
csd = sum(abs(wft),3);
LFP = spm_iwft(  wft(:,:,1),w,n);
lfp = spm_iwft(diag(K)*wft(:,:,1),w,n);
LFP = 4*LFP/std(LFP) + 16;
lfp = 4*lfp/std(lfp) + 16;

if Nt == 1, subplot(2,2,3), else subplot(4,1,2),end
imagesc(t,Hz,csd), axis xy, hold on
plot(t,lfp,'w',t,LFP,'w:'), hold off
title('time-frequency (and phase)response','FontSize',16)
xlabel('time (seconds)','FontSize',12), ylabel('frequency','FontSize',12)
if Nt == 1, axis square, end

% local field potentials
%==========================================================================
if Nt == 1, subplot(2,2,2), else subplot(4,1,3),end
plot(t,spm_cat(x),':'), hold on
plot(t,spm_cat(u)),     hold off
title('local field potentials','FontSize',16)
xlabel('time (updates)','FontSize',12)
ylabel('Response','FontSize',12), spm_axis tight
grid on, set(gca,'XTick',(1:(NT*Nt))*Nb*dt)
if Nt == 1, axis square, end

% simulated dopamine responses
%==========================================================================
if Nt == 1, subplot(2,2,4), else subplot(4,1,4),end
bar(spm_vec(da),1,'k'), title('phasic dopamine responses','FontSize',16)
xlabel('time (updates)','FontSize',12)
ylabel('change in precision','FontSize',12), spm_axis tight
if Nt == 1, axis square, end


 