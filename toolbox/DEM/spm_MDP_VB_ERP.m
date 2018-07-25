function [x,y,ind] = spm_MDP_VB_ERP(MDP,FACTOR,T)
% auxiliary routine for hierarchical electrophysiological responses
% FORMAT [x,y] = spm_MDP_VB_ERP(MDP,FACTOR,T)
%
% MDP    - structure (see spm_MDP_VB)
% FACTOR - the hidden factors (at the second alevel) to plot
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
% $Id: spm_MDP_VB_ERP.m 7382 2018-07-25 13:58:04Z karl $


% defaults: assume the first factor is of interest
%==========================================================================
try, f1 = FACTOR(1); catch, f1 = 1; end
try, f2 = FACTOR(2); catch, f2 = 1; end

% and T = 1
%--------------------------------------------------------------------------
if nargin < 3, T = 1; end

for m = 1:numel(MDP)

    % dimensions
    %----------------------------------------------------------------------
    xn  = MDP(m).xn{f1};      % neuronal responses
    Nb  = size(xn,1);         % number of time bins per epochs
    Nx  = size(xn,2);         % number of states
    Ne  = size(xn,3);         % number of epochs
    
    
    % expected hidden states
    %======================================================================
    x     = cell(Ne,Nx);
    y     = cell(Ne);
    for k = 1:Ne
        for j = 1:Nx
            x{k,j} = xn(:,j,T,k);
        end
        if isfield(MDP,'mdp')
            y{k}   = spm_MDP_VB_ERP(MDP(m).mdp(k),f2,1);
        else
            y{k}   = [];
        end
    end
    
    if nargin > 2, return, end
    
    % synchronise responses
    %----------------------------------------------------------------------
    u   = {};
    v   = {};
    uu  = spm_cat(x(1,:));
    for k = 1:Ne
        
        % low-level
        %------------------------------------------------------------------
        v{end + 1,1} = spm_cat(y{k});
        if k > 1
            u{end + 1,1} = ones(size(v{end,:},1),1)*u{end,1}(end,:);
        else
            u{end + 1,1} = ones(size(v{end,:},1),1)*uu(1,:);
        end
        
        % time bin indices
        %------------------------------------------------------------------
        ind(k) = size(u{end},1);
        
        % high-level
        %------------------------------------------------------------------
        u{end + 1,1} = spm_cat(x(k,:));
        v{end + 1,1} = ones(size(u{end,:},1),1)*v{end,1}(end,:);
        
        % time bin indices
        %------------------------------------------------------------------
        ind(k) = ind(k) + size(u{end},1);

    end
    
    % accumulate over trials
    %----------------------------------------------------------------------
    U{m,1} = u;
    V{m,1} = v;
    
end

% time bin (seconds)
%--------------------------------------------------------------------------
u  = spm_cat(U);
v  = spm_cat(V);
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

% higher-level unit responses
%--------------------------------------------------------------------------
factor = MDP(1).label.factor{f1};
name   = MDP(1).label.name{f1};

subplot(4,1,1), image(t,1:(size(u,2)),64*(1 - u')), ylabel('Unit')
title(sprintf('Unit reponses : %s',factor),'FontSize',16)
if numel(name) < 16
    grid on, set(gca,'YTick',1:numel(name))
    set(gca,'YTickLabel',name)
end

% lower-level unit responses
%--------------------------------------------------------------------------
factor = MDP(1).MDP(1).label.factor{f2};
name   = MDP(1).MDP(1).label.name{f2};

subplot(4,1,2), image(t,1:(size(v,2)),64*(1 - v')), ylabel('Unit')
title(sprintf('Unit reponses : %s',factor),'FontSize',16)
if numel(factor) < 16
    grid on, set(gca,'YTick',1:numel(name))
    set(gca,'YTickLabel',name)
end

% event related responses at both levels
%--------------------------------------------------------------------------
subplot(4,1,3), plot(t,x',t,y','-.')
title('Local field potentials','FontSize',16)
ylabel('Depolarisation'),spm_axis tight
grid on, set(gca,'XTick',(1:(length(t)/Nb))*Nb*dt)

