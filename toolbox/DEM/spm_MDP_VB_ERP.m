function [x,y,ind] = spm_MDP_VB_ERP(MDP,FACTOR,T)
% auxiliary routine for plotting hierarchical electrophysiological responses
% FORMAT [x,y] = spm_MDP_VB_ERP(MDP,FACTOR,T)
%
% MDP    - structure (see spm_MDP_VB)
% FACTOR - the hidden factors (at the first and second level) to plot
% T      - flag to return cell of expectations (at time T; usually 1)
%
% x      - simulated ERPs (high-level)
% y      - simulated ERPs (low level)
% ind    - indices or bins at the end of each (synchronised) epoch
%
% This routine combines first and second level hidden expectations by
% synchronising them; such that first level updating is followed by an
% epoch of second level updating - during which updating is suspended 
% (and expectations are held constant). The ensuing spike rates can be
% regarded as showing delay period activity. In this routine, simulated
% local field potentials are band pass filtered spike rates (between eight
% and 32 Hz).
%
% Graphics are provided for first and second levels, in terms of simulated
% spike rates (posterior expectations), which are then combined to show
% simulated local field potentials for both levels (superimposed).
%
% At the lower level, only expectations about hidden states in the first
% epoch are returned (because the number of epochs can differ from trial
% to trial).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_VB_ERP.m 7003 2017-02-02 18:22:56Z karl $
 
 
% defaults: assume the first factor is of interest
%==========================================================================
try, f = FACTOR(1); catch, f = 1; end
xn     = MDP.xn{f};
 
% dimensions
%--------------------------------------------------------------------------
Nb  = size(xn,1);         % number of time bins per epochs
Nx  = size(xn,2);         % number of states
Ne  = size(xn,3);         % number of epochs
 
% units to plot (if T is specified, extract a single time point)
%--------------------------------------------------------------------------
UNITS = [];
if nargin < 3
    for i = 1:1
        for j = 1:Nx
            UNITS(:,end + 1) = [j;i];
        end
    end
else
    for j = 1:Nx
        UNITS(:,end + 1) = [j;T];
    end
end
 
% expected hidden states
%==========================================================================
for k = 1:Ne
    for j = 1:size(UNITS,2)
        x{k,j} = xn(:,UNITS(1,j),UNITS(2,j),k);
    end
    if isfield(MDP,'mdp')
        try, f = FACTOR(2);  catch,  f = 1;  end
        y{k}   = spm_MDP_VB_ERP(MDP.mdp(k),f,1);
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
    v{end + 1,1} = spm_cat(y{k});
    if k > 1
        u{end + 1,1} = ones(size(v{end,:},1),1)*u{end,1}(end,:);
    else
        u{end + 1,1} = ones(size(v{end,:},1),1)*uu(1,:);
    end
    
    % time bin indices
    %----------------------------------------------------------------------
    ind(k) = size(u{end},1);
 
    % high-level
    %----------------------------------------------------------------------
    u{end + 1,1} = spm_cat(x(k,:));
    v{end + 1,1} = ones(size(u{end,:},1),1)*v{end,1}(end,:);
    
    % time bin indices
    %----------------------------------------------------------------------
    ind(k) = ind(k) + size(u{end},1);
    
end
 
% time bin (seconds)
%--------------------------------------------------------------------------
u  = spm_cat(u);
v  = spm_cat(v);
dt = 1/64;
t  = (1:size(u,1))*dt;
 
% bandpass filter between 8 and 32 Hz
%--------------------------------------------------------------------------
c  = 1/32;
x  = log(u + c);
y  = log(v + c);
x  = spm_conv(x,2,0) - spm_conv(x,16,0);
y  = spm_conv(y,2,0) - spm_conv(y,16,0);
 
if nargout > 2, return, end
 
% simulated firing rates and the local field potentials
%==========================================================================
subplot(4,1,1), imagesc(t,1:(size(u,2)),1 - u')
title('Unit responses (high-level)','FontSize',16), ylabel('Unit')
 
subplot(4,1,2), imagesc(t,1:(size(v,2)),1 - v')
title('Unit responses (low-level)' ,'FontSize',16), ylabel('Unit')
 
subplot(4,1,3), plot(t,x',t,y',':')
title('Local field potentials','FontSize',16)
ylabel('Depolarisation'),spm_axis tight
grid on, set(gca,'XTick',(1:(length(t)/Nb))*Nb*dt)

