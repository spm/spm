function DEM_HB_and_LE
%--------------------------------------------------------------------------
% This routine is a numerical examination of the relationship between
% entropy, mutual information and the exponential divergence of
% trajectories as the Rayleigh parameter of a Lorenz attractoris increased
% - through a pitchfork bifurcation and subsequent (subcritical) Hopf
% bifurcation. The (stochastic) Lorentz system is integrated for different
% values of the Rayleigh parameter. The nonequilibrium steady-state density
% is then estimated by embedding into a discrete state space; while the
% bifurcations are characterised in terms of the maximal Lyapunov exponent.
% The key thing to observe is the decrease in entropy of blanket states
% prior to the Hopf bifurcation and implicit exponential divergence of
% trajectories. This is scored by the maximal Lyapunov exponent crossing
% zero. Here, the form of the Lorenz attractor defines the three states as
% active, sensory and hidden. Note that there are no internal states in
% this example and blanket states become the particular states (i.e., the
% states of a particle).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_HB_and_LE.m 7502 2018-12-02 12:28:03Z karl $

% generative model
%==========================================================================                       % switch for demo
spm_figure('GetWin','DEM'); clf

% flow and Jacobian functions
%--------------------------------------------------------------------------
f  = @(x,v,P,G) v(:) + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64;
Df = @(x,v,P,G) [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]/64 + ...
                [0 0 0; -x(3) 0 0; 0 x(1) 0]/64;

% model with initial states
%--------------------------------------------------------------------------
G.x  = [1; 1; 24];
G.f  = f;

% set up
%--------------------------------------------------------------------------
b     = 32;                         % number of bins for density estimation
T     = 2^20;                       % length of trajectory
W     = 32;                         % precision of intrinsic fluctuations
P     = exp((-16:32)*log(32)/32);   % Rayleigh parameter range
for k = 1:length(P)
    
    % integrated timeseries
    %----------------------------------------------------------------------
    U.u   = randn(T + 256,3)/W;
    Pk    = [10; -8/3; P(k)];
    t     = spm_int_L(Pk,G,U);
    
    % remove intial transients
    %----------------------------------------------------------------------
    t     = t(256:end,:);
    t     = t(1:T,:);
    tt{k} = t;
    
    
    % sample density
    %----------------------------------------------------------------------
    for i = 1:3
        t(:,i) = t(:,i) - min(t(:,i));
        t(:,i) = (b - 2)*t(:,i)/max(t(:,i));
    end
    t     = t + 1;
    p     = zeros(b,b,b) + 1;
    S     = 0;
    for i = 1:T
        
        % accumulate in bins
        %------------------------------------------------------------------
        j = round(t(i,:));
        p(j(3),j(2),j(1)) = p(j(3),j(2),j(1)) + 1;
        
        % Lyapunov exponents
        %------------------------------------------------------------------
        dfdx = Df(tt{k}(i,:)',zeros(1,3)',Pk);
        S    = S + sort(real(eig(dfdx,'nobalance')),'descend');
        
    end
    p     = p/sum(p(:));
    pp{k} = p;
    
    % mutual informations of sample density
    %----------------------------------------------------------------------
    [I,Ii,Ie] = spm_self_entropy(p);
    MI(1,k)   = I;
    MI(2,k)   = Ii;
    MI(3,k)   = Ie;
    LE(:,k)   = S/T; 

    
    % plot mutual informations
    %----------------------------------------------------------------------
    subplot(3,2,3)
    semilogx(P(1:k),MI(1,:),'b',P(1:k),MI(2,:),'b-.',P(1:k),MI(3,:),'b:')
    axis square xy
    title('Expected surpise','Fontsize',16)
    xlabel('Control parameter'),ylabel('Entropies (nats)'),drawnow
    
    % plot maximal Lyapunov exponent
    %----------------------------------------------------------------------
    subplot(3,2,4)
    semilogx(P(1:k),LE(1,:),'b')
    axis square xy
    title('Lyapunov exponent','Fontsize',16)
    xlabel('Control parameter'),ylabel('Principal exponent'),drawnow
    
end

% lines and thresholds
%--------------------------------------------------------------------------
j     = find(LE(1,:) > 0,1);
j     = P(j);
subplot(3,2,3)
hold on, plot([j j],[ 0 6],':'), hold off
hold on, plot([1 1],[ 0 6],':'), hold off
legend({'H(B)','H(B|E)','I(B,E)'})
subplot(3,2,4)
hold on, plot([j j],[-1 1]*0.03,':'), hold off
hold on, plot([1 1],[-1 1]*0.03,':'), hold off
hold on, plot(P,zeros(size(P)),'--'), hold off


% illustrate trajectories and ergodic density
%--------------------------------------------------------------------------
j     = [8 38 length(P)];
for i = 1:length(j)
    
    % exemplar trajectory (plot)
    %----------------------------------------------------------------------
    subplot(3,3,i)
    t    = tt{j(i)};
    plot(t(1:1024,2),t(1:1024,3),'k')
    axis square xy
    title('Trajectory','Fontsize',16)
    axis([-30 30 -5 60])
    
    % image format
    %----------------------------------------------------------------------
    subplot(3,3,6 + i)
    p = pp{j(i)};
    imagesc(1-squeeze(sum(p,2))'),axis xy square
    title('Marginal density','Fontsize',16)
    ylabel('State'), xlabel('State')
    
%     subplot(3,3,6 + i)
%     [W,S,X,beta] = spm_power_law(t');
%     plot(W,S,'b.',W,X*beta,'b','LineWidth',1)
%     title(sprintf('alpha = %-2.2f',beta(2)),'FontSize',16)
%     ylabel('Log power'), xlabel('Log frequency')
%     axis square, axis xy, spm_axis tight
    
    
end

return


function [HB,HBH,IBH] = spm_self_entropy(pHxB)
% FORMAT [HB,HBH,IBH] = spm_self_entropy(pHxB)
% Entropies
% HS  = H(B)              % self entropy
% HBH = H(B|H)            % conditional entropy
% IBH = I(B,H)            % mutual information
%
% This subroutine assumes that the first dimension of the joint density
% corresponds to a hidden or external state and the rest are particular
% or blanket states

% evaluate joint density and posterior
%--------------------------------------------------------------------------
pB    = sum(pHxB,1);
pH    = sum(sum(sum(pHxB,2),3),4);

% inline functions
%--------------------------------------------------------------------------
ln    = @(p)log(spm_vec(p) + 1e-16);
H     = @(p)-spm_vec(p)'*ln(p);

% relative entropies
%--------------------------------------------------------------------------
HB    = H(pB);
IBH   = H(pH) + H(pB) - H(pHxB);
HBH   = HB - IBH;

return

function D = spm_Kaplan_Yorke(LE)
% FORMAT Kaplan Yorke estimate of dimensional complexity
%--------------------------------------------------------------------------
L     = real(LE);
for i = 1:size(L,2)
    l = sort(L(:,i),'descend');
    j = find(cumsum(l) > 0,1);
    if isempty(j), j = 0; end
    D(i) = j + sum(l(1:j))/abs(l(j + 1));
end


function [W,S,X,beta] = spm_power_law(x)
% FORMAT spm_power_law(x)

% illustrate power law scaling
%--------------------------------------------------------------------------
N     = floor(log2(size(x,2)));
s     = abs(fft(x(1,:)')).^2;
w     = (1:2^12)';
W     = w;
S     = s(w + 1);

S     = decimate(log(S),N - 4);
W     = log(decimate(W,N - 4));
X     = [ones(size(W)),W];

% plot part of trajectory
%--------------------------------------------------------------------------
[~,i] = max(abs(diff(spm_conv(x(1,:),2^(N - 8)))));
nn    = 2^10;
i     = (-nn:nn) + i;
i     = i(i > 0 & i < size(x,2));

% estimate exponent (alpha)
%--------------------------------------------------------------------------
[~,~,beta] = spm_ancova(X,[],S,[0;1]);
