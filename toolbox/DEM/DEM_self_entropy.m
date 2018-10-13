function DEM_self_entropy
%--------------------------------------------------------------------------
% Routine to produce graphics illustrating self organisation in terms of
% the entropy of blanket states and associated trajectories. A low blanket
% entropy induces anomalous diffusion and itinerancy with power law scaling
% (i.e., self similar dynamics). This example uses a fixed form (quadratic)
% likelihood and optimises the density over over hidden states to minimise
% blanket (i.e., self) entropy explicitly.
%
% In this example, there is just one active and sensory state and one
% hidden state to illustrate noise-phase symmetry breaking as the
% probability density over action reduces the entropy of sensory states
% under a fixed density of hidden or external states.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_self_entropy.m 7446 2018-10-13 15:28:42Z karl $


% set up
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
rng('default')

n     = 128;                            % number of (SxH) bins
m     = 16;                             % numer of (A) active bins
dt    = 1;                              % time step for solution
N     = 18;                             % 2^N solution

% likelihood – mapping from hidden states to sensory states - A = p(S|HxA)
%--------------------------------------------------------------------------
a     = zeros(n,n);
for i = 1:n
    for j = 1:n
        es     = (j - n/2).^2/32 + n/8;
        vs     = (es/8)^2;
        a(i,j) = exp(-(i - es).^2/vs);
    end
end

% reduce entropy for different levels of active states and normalise
%--------------------------------------------------------------------------
A     = zeros(n,n,m);
for k = 1:m
    s        = 1 + 8*(k - 1);
    A(:,:,k) = spm_softmax(log(a)/s);
end
for j = 1:n
    for k = 1:m
        A(:,j,k) = A(:,j,k)/sum(A(:,j,k) + exp(-32));
    end
end

% check marginal likehoods
%--------------------------------------------------------------------------
% for i = 1:min(m,16)
%     subplot(4,4,i)
%     imagesc(1 - A(:,:,i)), axis square, axis xy
% end
% return

% hidden and active states
%--------------------------------------------------------------------------
pH    = exp(-((1:n) - n/2 - n/8).^2/n*2)/2 + exp(-((1:n) - n/2 + n/8).^2/n);
lnpH  = log(pH(:));
lnpA  = ones(m,1);

% progressively optimise mutual information w.r.t. hidden states
%==========================================================================
ni    = [1 3 4];                        % batches of iterations
for g = 1:3
    for ii = 1:ni(g)
        
        % evaluate entropy for current hidden states
        %------------------------------------------------------------------
        [HB,Gi,Ge] = spm_G(A,lnpH,lnpA);
        
        % evaluate joint density and marginals
        %------------------------------------------------------------------
        pH     = spm_softmax(lnpH);
        pA     = spm_softmax(lnpA);
        
        % joint probability
        %--------------------------------------------------------------------------
        for i = 1:n
            for j = 1:m
                pSxHxA(:,i,j) = A(:,i,j)*pH(i)*pA(j);
            end
        end
        pSxH   = sum(pSxHxA,3);
        pS     = sum(pSxH,2);
        
        % mutual informations
        %------------------------------------------------------------------
        Hh(g,1) = HB;
        He(g,1) = Ge;
        Hi(g,1) = Gi;
        
        
        % Optimise self-entropy w.r.t. distribution over hidden states
        %------------------------------------------------------------------
        dp   = spm_diff(@spm_G,A,lnpH,lnpA,3);
        
        % update marginal over hidden states
        %------------------------------------------------------------------
        lnpA = log(spm_softmax(lnpA - dp(:)*8));
        
        % graphics
        %------------------------------------------------------------------
        subplot(3,2,1),     imagesc(1 - sum(A,3))
        title('Likelihood','FontSize',16)
        xlabel('Hidden states'), ylabel('Sensory states')
        axis square, axis xy
        
        subplot(3,2,2),     bar([Hh,He,Hi])
        title('Expected surpise','FontSize',16)
        xlabel('Iteration'),ylabel('Information (nats)')
        axis square, axis xy,
        legend({'H(B)','H(B|E)','I(B,E)'})
        set(gca,'XTickLabel',ni)
        
        subplot(3,3,g + 3), imagesc(1 - pSxH)
        j  = sum(ni(1:(g - 1))) + ii;
        title(sprintf('Iteration %i',j),'FontSize',16)
        xlabel('Hidden states'), ylabel('Sensory states')
        axis square, axis xy
        
        hold on
        plot(pH*n*n/8,'k')
        plot(pS*n*n/8,(1:n),'r')
        hold off
        drawnow
        
    end
    
    % illustrate dynamics
    %======================================================================
    
    % flow
    %----------------------------------------------------------------------
    G       = eye(2,2);                  % amplitude of random fluctuations
    Q       = [0 -1;1 0]/4;              % solenoidal flow
    [gh,gs] = gradient(log(pSxH));       % gradients
    f       = [gh(:),gs(:)]*(G - Q);     % flow
    fh      = spm_unvec(f(:,1),gh);
    fs      = spm_unvec(f(:,2),gs);
    [gi,gj] = meshgrid(1:n,1:n);
    i       = 1:8:n;
    
    subplot(3,2,5), hold off, quiver(gi(i,i),gj(i,i),fh(i,i),fs(i,i),'k')
    title('Flow and trajectories','FontSize',16)
    xlabel('Hidden states'), ylabel('Blanket states')
    axis square, axis xy
    
    
    % solve for a particular trajectory
    %----------------------------------------------------------------------
    [~,k] = max(pSxH(:));
    [p,q] = ind2sub([n,n],k);
    x     = [q;p];
    for t = 1:(2^N)
        x(:,t)     = max(1,min(n,x(:,t)));
        k          = sub2ind([n,n],round(x(2,t)),round(x(1,t)));
        dx         = f(k,:)' + sqrt(G/2)*randn(2,1);
        x(:,t + 1) = x(:,t)  + dx*dt;
    end
    
    % illustrate power law scaling (sensory states)
    %----------------------------------------------------------------------
    s     = abs(fft(x(1,:)')).^2;
    w     = (1:2^12)';
    W     = w;
    S     = s(w + 1);
    
    S     = decimate(log(S),N - 4);
    W     = log(decimate(W,N - 4));
    X     = [ones(size(W)),W];
    
    % plot part of trajectory
    %----------------------------------------------------------------------
    [~,i] = max(abs(diff(spm_conv(x(1,:),2^(N - 8)))));
    nn    = 2^10;
    i     = (-nn:nn) + i;
    i     = i(i > 0 & i < size(x,2));
    hold on, plot(x(1,i),x(2,i),'b'), hold off
    axis([1,n,1,n]),drawnow
    
    % estimate exponent (alpha)
    %----------------------------------------------------------------------
    [~,~,beta] = spm_ancova(X,[],S,[0;1]);
    
    % plot
    %----------------------------------------------------------------------
    if g == 3
        subplot(3,2,6), plot(W,S,'b.',W,X*beta,'b','LineWidth',1)
    elseif g == 1
        subplot(3,2,6), plot(W,S - 4,'c.'), hold on
    else
        subplot(3,2,6), plot(W,S - 2,'m.'), hold on
    end
    title(sprintf('alpha = %-2.2f',beta(2)),'FontSize',16)
    ylabel('Log power'), xlabel('Log frequency')
    axis square, axis xy, spm_axis tight
    
end

return

% subroutines
%==========================================================================

function [HB,HBH,IBH] = spm_G(A,lnpH,lnpA)
% FORMAT [HB,HBH,IBH] = spm_G(A,lnpH,lnpA)
% HB  = H(B)              % self entropy
% HBH = H(B|H)            % conditional entropy
% IBH = I(H,B)            % mutual information


% evaluate marginal over hidden states
%--------------------------------------------------------------------------
n     = size(A,2);
m     = size(A,3);
pH    = spm_softmax(lnpH);
pA    = spm_softmax(lnpA);

% inline functions
%--------------------------------------------------------------------------
ln    = @(p)log(spm_vec(p) + 1e-16);
H     = @(p)-spm_vec(p)'*ln(p);

% joint probability
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:m
        pSxHxA(:,i,j) = A(:,i,j)*pH(i)*pA(j);
    end
end

% marginal blanket
%--------------------------------------------------------------------------
pB    = spm_vec(sum(pSxHxA,2));
HB    = H(pB);
if nargout > 1
    IBH = H(pB) + H(pH) - H(pSxHxA);
    HBH = HB - IBH;
end

return