function FEP_self_entropy
% This demonstration uses an ensemble of particles with intrinsic (Lorentz
% attractor) dynamics and (Newtonian) short-range coupling.  This routine
% illustrates self organisation in terms of the entropy of blanket states
% (and concomitant changes in terms of mutual information (i.e., complexity
% cost or risk). Here, the ensemble average of these entropy measures is
% taken over all (128) particles of macromolecules; where the Markov
% blanket of each particle comprises all but the third (electrochemical)
% hidden state. The graphics produced by this routine simply plot the
% decrease in blanket entropy (and complexity cost) as the system
% approaches its random dynamical attractor. Illustrative trajectories of
% the particles are provided at three points during the (stochastic)
% chaotic transient.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: FEP_self_entropy.m 7679 2019-10-24 15:54:07Z spm $


% default settings (GRAPHICS sets movies)
%--------------------------------------------------------------------------
rng('default')

% Demo of synchronization manifold using coupled Lorenz attractors
%==========================================================================
N    = 128;                         % number of (Lorenz) oscillators
T    = 1024;                        % number of time bins
dt   = 1/32;                        % time interval

% parameters
%--------------------------------------------------------------------------
P.k  = 1 - exp(-rand(1,N)*4);       % variations in temporal scale
P.d  = 1/8;                         % amplitude of random fluctuations

% states
%--------------------------------------------------------------------------
x.p  = randn(2,N)*4;                % microstates (position)
x.v  = zeros(2,N);                  % microstates (velocity)
x.q  = randn(3,N)/32;               % microstates (states)
u    = zeros(1,T);                  % exogenous fluctuations


% generate an dynamics from initial conditions
%==========================================================================
spm_figure('GetWin','Markov blanket');clf

[Q,X,V]  = spm_soup(x,u,P,T,dt,1);

% blanket states - noting the third state is an internal state
%--------------------------------------------------------------------------
% Q    - history of microstates (states)
% X    - history of microstates (position)
% V    - history of microstates (velocity)

for i = 1:size(X,3)
    S(:,:,i) = [Q(1:2,:,i);X(:,:,i);V(:,:,i)];
end

% self organisation in terms of blanket entropy
%==========================================================================

% Sample blanket states for np particles
%--------------------------------------------------------------------------
wt    = 256;                        % sliding window length
it    = 0:32:1024;                  % window spacing
nw    = 16;                         % number of windows
np    = 128;                        % number of particles
for i = 1:np
    
    % time window
    %----------------------------------------------------------------------
    for tt = 1:nw
        t  = it(tt);
        
        % get blanket and external states
        %------------------------------------------------------------------
        clear b e d
        for j = 1:wt
            ni     = 1:size(S,2);
            ni(i)  = [];
            b(j,:) = spm_vec(S(:, i,t + j))';
            e(j,:) = spm_vec(S(:,ni,t + j))';
            
            % Euclidean distance
            %--------------------------------------------------------------
            d(j,:) = (X(1,i,t + j) - X(1,ni,t + j)).^2 + ...
                     (X(2,i,t + j) - X(2,ni,t + j)).^2;
        end
        b     = squeeze(b);
        
        % retain nearest external states
        %------------------------------------------------------------------
        [d,j] = sort(mean(d));
        e     = e(:,j(1:(wt/4)));
        
        % evaluate covariances
        %------------------------------------------------------------------
        Cbe   = cov([b e]);
        Cb    = cov(b);
        Ce    = cov(e);
        
        % and relative entropies
        %------------------------------------------------------------------
        Hbxe  = spm_logdet(Cbe)/2;
        Hb    = spm_logdet(Cb)/2;
        He    = spm_logdet(Ce)/2;
        
        HB(i,tt)  = Hb;
        IBE(i,tt) = Hb + He - Hbxe;
        HBE(i,tt) = Hbxe - He;
        
    end
    
end

% accumulate and plot
%--------------------------------------------------------------------------
subplot(3,1,2)
plot(it(1:nw),HB','b:',it(1:nw),mean(HB),'b')
xlabel('Time','FontSize',12),ylabel('Nats','FontSize',12)
title('Blanket entropy','FontSize',16), spm_axis tight
set(gca,'YLim',[-8 8])

subplot(3,1,3)
plot(it(1:nw),mean(HB),'b',it(1:nw),mean(HBE),'b-.',it(1:nw),mean(IBE),'b:')
xlabel('Time','FontSize',12),ylabel('Nats','FontSize',12)
title('Ensemble entropies and complexity cost','FontSize',16), spm_axis tight
legend({'H(B)','H(B|E)','I(B,E)'})

% illustrate dynamics
%--------------------------------------------------------------------------
tt    = [64,256,512];
for i = 1:3
    t = tt(i) + (1:wt);
    subplot(3,3,i)
    for j = 1:np
        plot(squeeze(X(1,j,t)),squeeze(X(2,j,t))), hold on
    end
    axis([-1 1 -1 1]*8), axis square, hold off
    xlabel('Position','FontSize',12)
    ylabel('Position','FontSize',12)
    str = sprintf('Trajectories (t = %i)',tt(i));
    title(str,'FontSize',16)
end

return
