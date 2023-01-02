function [v] = spm_MDP_VB_LFP(MDP,UNITS,f,SPECTRAL)
% auxiliary routine for plotting simulated electrophysiological responses
% FORMAT [v] = spm_MDP_VB_LFP(MDP,UNITS,FACTOR,SPECTRAL)
%
% MDP        - structure (see spm_MDP_VB_X.m)
%  .xn       - neuronal firing
%  .dn       - phasic dopamine responses
%
% UNITS(1,:) - hidden state                           [default: all]
% UNITS(2,:) - time step
%
% FACTOR     - hidden factor to plot                  [default: 1]
% SPECTRAL   - replace raster with spectral responses [default: 0]
%
% v - selected unit responses {number of trials, number of units}
%
% This routine plots simulated electrophysiological responses. Graphics are
% provided in terms of simulated spike rates (posterior expectations).
%
% see also: spm_MDP_VB_ERP (for hierarchical belief updating)
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
% check for simulated neuronal responses
%--------------------------------------------------------------------------
if ~isfield(MDP(1),'xn')
    warning ('please use another inversion scheme that simulates neuronal responses (e.g., spm_MDP_VB_XX)')
    return
end

% defaults
%==========================================================================
try, f;          catch, f        = 1;  end
try, UNITS;      catch, UNITS    = []; end
try, SPECTRAL;   catch, SPECTRAL = 0;  end
try, MDP = spm_MDP_check(MDP);         end

% dimensions
%--------------------------------------------------------------------------
Nt     = length(MDP);               % number of trials
try
    Ne = size(MDP(1).xn{f},4);      % number of epochs
    Nx = size(MDP(1).B{f}, 1);      % number of states
    Nb = size(MDP(1).xn{f},1);      % number of time bins per epochs
catch
    Ne = size(MDP(1).xn,4);         % number of epochs
    Nx = size(MDP(1).A, 2);         % number of states
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
if  size(ALL,2) > 512
    ii  = round(linspace(1,size(ALL,2),512));
    ALL = ALL(:,ii);
end
if isempty(UNITS)
    UNITS = ALL;
end
ii = 1:size(ALL,2);

% summary statistics: firing rates
%==========================================================================
for i = 1:Nt
    
    % for all units
    %----------------------------------------------------------------------
    str    = {};
    try
        xn = MDP(i).xn{f};
    catch
        xn = MDP(i).xn;
    end
    for j = 1:size(ALL,2)
        for k = 1:Ne
            z{i,1}{k,j} = xn(:,ALL(1,j),ALL(2,j),k);
        end
        str{j} = sprintf('%s: t=%i',MDP(1).label.name{f}{ALL(1,j)},ALL(2,j));
    end
    
    % for selected units
    %----------------------------------------------------------------------
    for j = 1:size(UNITS,2)
        for k = 1:Ne
            v{i,1}{k,j} = xn(:,UNITS(1,j),UNITS(2,j),k);
        end
    end

    
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

% simulated firing rates 
%--------------------------------------------------------------------------
z = spm_cat(z)';                         % firing rates of all units
v = spm_cat(v)';                         % firing rates of selected units

% bandpass filter log rates between 8 and 32 Hz: local field potential
%--------------------------------------------------------------------------
c  = 1/32;
x  = log(z' + c);
u  = log(v' + c);
x  = spm_conv(x,2,0) - spm_conv(x,16,0);
u  = spm_conv(u,2,0) - spm_conv(u,16,0);

% simulated firing rates
%--------------------------------------------------------------------------
if Nt == 1, subplot(3,2,1), else, subplot(4,1,1),end
image(t,ii,64*(1 - z))
title(MDP(1).label.factor{f},'FontSize',16)
xlabel('time (sec)','FontSize',12)

if numel(str) < 16
   grid on, set(gca,'YTick',1:(Ne*Nx))
   set(gca,'YTickLabel',str)
end
grid on, set(gca,'XTick',(1:(Ne*Nt))*Nb*dt)
if Ne*Nt > 32, set(gca,'XTickLabel',{}), end
if Nt == 1,    axis square,              end
 
% time frequency analysis and theta phase
%--------------------------------------------------------------------------
wft = spm_wft(x,w,n);
csd = sum(abs(wft),3);
lfp = sum(x,2);
phi = spm_iwft(sum(wft(1,:,:),3),w(1),n);
lfp = 4*lfp/std(lfp) + 16;
phi = 4*phi/std(phi) + 16;
 
if Nt == 1, subplot(3,2,3), else, subplot(4,1,2),end
imagesc(t,Hz,csd), axis xy, hold on
plot(t,lfp,'w:',t,phi,'w'), hold off
grid on, set(gca,'XTick',(1:(Ne*Nt))*Nb*dt)

title('Time-frequency response','FontSize',16)
xlabel('time (sec)'), ylabel('frequency (Hz)')
if Nt == 1, axis square, end

% spectral responses
%--------------------------------------------------------------------------
if SPECTRAL
    
    % spectral responses (for each unit)
    %----------------------------------------------------------------------
    if Nt == 1, subplot(3,2,1), else, subplot(4,2,1),end
    csd = squeeze(sum(abs(wft),2));
    plot(Hz,log(squeeze(csd)))
    title('Spectral response','FontSize',16)
    xlabel('frequency (Hz)'),
    ylabel('log power')
    spm_axis tight, box off, axis square
    
    % amplitude-to-amplitude coupling (average over units)
    %----------------------------------------------------------------------
    if Nt == 1, subplot(3,2,2), else, subplot(4,2,2),end
    cfc   = 0;
    for i = 1:size(wft,3)
        cfc = cfc + corr((abs(wft(:,:,i)))');
    end
    imagesc(Hz,Hz,cfc)
    title('Cross-frequency coupling','FontSize',16)
    xlabel('frequency (Hz)'),
    ylabel('frequency (Hz)')
    box off, axis square

end
 
% local field potentials
%==========================================================================
if Nt == 1, subplot(3,2,4), else, subplot(4,1,3),end
plot(t,u),     hold off, spm_axis tight, a = axis;
plot(t,x,':'), hold on
plot(t,u),     hold off, axis(a)
grid on, set(gca,'XTick',(1:(Ne*Nt))*Nb*dt), 
for i = 2:2:Nt
    h = patch(((i - 1) + [0 0 1 1])*Ne*Nb*dt,a([3,4,4,3]),-[1 1 1 1],'w');
    set(h,'LineStyle',':','FaceColor',[1 1 1] - 1/32);
end
title('Local field potentials','FontSize',16)
xlabel('time (sec)')
ylabel('response')
if Nt == 1, axis square, end, box off

% firing rates
%==========================================================================
if Nt == 1, subplot(3,2,2)
    plot(t,v),     hold on, spm_axis tight, a = axis;
    plot(t,z,':'), hold off
    grid on, set(gca,'XTick',(1:(Ne*Nt))*Nb*dt), axis(a)
    title('Firing rates','FontSize',16)
    xlabel('time (sec)')
    ylabel('response')
    axis square
end

% simulated dopamine responses (if not a moving policy)
%==========================================================================
if Nt == 1, subplot(3,2,6), else, subplot(4,1,4),end
dn    = spm_vec(dn);
dn    = dn.*(dn > 0);
dn    = dn + (dn + 1/16).*rand(size(dn))/8;
bar(dn,1,'k'), title('Dopamine responses','FontSize',16)
xlabel('time (updates)')
ylabel('change in precision'), spm_axis tight, box off
YLim = get(gca,'YLim'); YLim(1) = 0; set(gca,'YLim',YLim);
if Nt == 1, axis square, end

% simulated rasters
%==========================================================================
Nr    = 16;
if Nt == 1 && numel(ii) < 129
    subplot(3,2,5)
    
    R  = kron(z,ones(Nr,Nr));
    R  = rand(size(R)) > R*(1 - 1/16);
    imagesc(t,1:(Nx*Ne),R),title('Unit firing','FontSize',16)
    xlabel('time (sec)','FontSize',12)
    
    grid on, set(gca,'XTick',(1:(Ne*Nt))*Nb*dt)
    grid on, set(gca,'YTick', 1:(Ne*Nx))
    set(gca,'YTickLabel',str), axis square
    
end

 
