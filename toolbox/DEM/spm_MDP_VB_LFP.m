function [u,v] = spm_MDP_VB_LFP(MDP,UNITS,FACTOR)
% auxiliary routine for plotting simulated electrophysiological responses
% FORMAT [u,v] = spm_MDP_VB_LFP(MDP,UNITS,FACTOR)
%
% u - selected unit rate of change of firing (simulated voltage)
% v - selected unit responses {number of trials, number of units}
%
% MDP - structure (see spm_MDP_VB
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_VB_LFP.m 6657 2015-12-31 17:59:31Z karl $
 
 
% defaults
%==========================================================================
try, f = FACTOR; catch, f = 1;      end
try, UNITS;      catch, UNITS = []; end


% dimensions
%--------------------------------------------------------------------------
Nt     = length(MDP);               % number of trials
Ne     = size(MDP(1).V,1) + 1;      % number of epochs
try
    Nx = size(MDP(1).B{f},1);       % number of states
    Nb = size(MDP(1).xn{f},1);      % number of time bins per epochs
catch
    Nx = size(MDP(1).A,2);          % number of states
    Nb = size(MDP(1).xn,1);         % number of time bins per epochs
end

% units to plot
%--------------------------------------------------------------------------
ALL   = [];
for i = 1:Ne
    for j = 1:Nx
        ALL(:,end + 1) = [j;i];
    end
end
if isempty(UNITS)
    UNITS = ALL;
end
    
% summary statistics
%==========================================================================
for i = 1:Nt
    
    % all units
    %----------------------------------------------------------------------
    try
        xn = MDP(i).xn{f};
    catch
        xn = MDP(i).xn;
    end
    for j = 1:size(ALL,2)
        for k = 1:Ne
            zj{k,j} = xn(:,ALL(1,j),ALL(2,j),k);
            xj{k,j} = gradient(zj{k,j}')';
        end
    end
    z{i,1} = zj;
    x{i,1} = xj;
    
    % selected units
    %----------------------------------------------------------------------
    for j = 1:size(UNITS,2)
        for k = 1:Ne
            vj{k,j} = xn(:,UNITS(1,j),UNITS(2,j),k);
            uj{k,j} = gradient(vj{k,j}')';
        end
    end
    v{i,1} = vj;
    u{i,1} = uj;
    
    % dopamine or changes in precision
    %----------------------------------------------------------------------
    dn(:,i) = mean(MDP(i).dn,2);
    
end

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

% simulated dopamine responses
%==========================================================================
if Nt == 1, subplot(3,1,3), else subplot(4,1,4),end
bar(spm_vec(dn),1,'k'), title('Phasic dopamine responses','FontSize',16)
xlabel('time (updates)','FontSize',12)
ylabel('change in precision','FontSize',12), spm_axis tight
if Nt == 1, axis square, end
 
