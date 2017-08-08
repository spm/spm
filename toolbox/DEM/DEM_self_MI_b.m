function DEM_self_MI_b
%--------------------------------------------------------------------------
% Routine to produce graphics illustrating self relative entropy or mutual
% information. A low self mutual information induces anomalous diffusion
% and itinerancy with power law scaling (i.e., self similar dynamics). This
% example reduces the KL divergence between a density with low expected
% self entropy and the density produced by hidden states.
%
% In this example where is just one Markov blanket states and one hidden
% state to illustrate noise phase symmetry breaking as self mutual
% information decreases. The subroutines illustrate the relationship
% between self mutual information, intrinsic mutual information and
% extrinsic cost.

% set up
%--------------------------------------------------------------------------
rng('default'), clf

n     = 64;                             % number of bins
S     = (1:n)'/n;                       % domain of states
dt    = 1;                              % time step for solution
N     = 20;                             % 2^N solution

% likelihood – mapping from hidden states to sensory states - A
%--------------------------------------------------------------------------
v     = @(s) ((1 - s).^2)/32 + 1/512;
m     = @(s) s.^2 + 1/8;
for i = 1:n
    for j = 1:n
        A(:,j) = exp(-(S - m(j/n)).^2/v(j/n));
    end
end
A     = A/diag(sum(A));                 % likelihood

% initial probability density over hidden states
%--------------------------------------------------------------------------
lnpH  = log(hanning(n))*2;

% progressively optimise mutual information w.r.t. hidden states
%==========================================================================
ni    = [1 7 8];                        % number of iterations (batches)
for g = 1:3
    for i = 1:ni(g)
        
        % evaluate self mutual information for current hidden states
        %------------------------------------------------------------------
        [G,Gi,Ge] = spm_G(A,lnpH);
        
        % evaluate joint density and marginals
        %------------------------------------------------------------------
        pH     = spm_softmax(lnpH);
        pS     = A*pH;
        pSxH   = A*diag(pH);
        pSxH   = pSxH/sum(sum(pSxH));
        
        % mutual informations
        %------------------------------------------------------------------
        MIi(g) = pS'*Gi;
        MI2(g) = pS'*Ge;
        MI3(g) = G;
        H(g)   = pS'*(-log(pS));
                
        % Optimise marginal w.r.t. KL divergence
        %------------------------------------------------------------------
        dp     = spm_diff(@spm_KL,A,lnpH,2);
        
        % update marginal over hidden states
        %------------------------------------------------------------------
        lnpH = log(spm_softmax(lnpH + dp(:)*16));
        
        
        % graphics
        %==================================================================
        subplot(3,2,1),     imagesc(1 - A)
        title('Likelihood','FontSize',16)
        xlabel('Hidden states'), ylabel('Blanket states')
        axis square, axis xy
        
        subplot(3,2,2),     bar([MIi;MI2;MI3;H]')
        title('Mutual information','FontSize',16)
        xlabel('Iteration'),ylabel('Mutual information')
        axis square, axis xy, legend({'iMI','MI2','MI3','H'})
        set(gca,'XTickLabel',ni)
        
        subplot(3,3,g + 3), imagesc(1 - pSxH)
        j  = sum(ni(1:(g - 1))) + i;
        title(sprintf('Iteration %i',j),'FontSize',16)
        xlabel('Hidden states'), ylabel('Blanket states')
        axis square, axis xy
        
        hold on
        tS  = spm_softmax((Gi - Ge)*4);
        plot(pH*n*n/8,'k')
        plot(pS*n*n/8,(1:n),'r')
        plot(tS*n*n/8,(1:n),'r:')
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
    
    % illustrate power law scaling
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
        subplot(3,2,6), plot(W,S - 2,'c.'), hold on
    else
        subplot(3,2,6), plot(W,S - 4,'m.'), hold on
    end
    title(sprintf('alpha = %-2.2f',beta(2)),'FontSize',16)
    ylabel('Log power'), xlabel('Log frequency')
    axis square, axis xy, spm_axis tight
    
end

return


function [I3,Gi,Ge] = spm_G(A,lnpH)
% target distribution: intrinsic MI minus KL
% G = Ge - Gi
% E[Gi] = I(H,S'|S)
% E[Ge] = I(H,S)
% E[G]  = I(H,S',S) 

% evaluate marginals and joint density 
%--------------------------------------------------------------------------
n     = size(A,2);
pH    = spm_softmax(lnpH);
pS    = A*pH;
pSxH  = A*diag(pH);
pSxH  = pSxH/sum(sum(pSxH));
pHS   = pSxH'/diag(sum(pSxH,2) + eps);


% entropies and probabilities
%--------------------------------------------------------------------------
for j = 1:n
    ph      = pHS(:,j);
    ps      = A*ph;
    psxh    = A*diag(ph);
    psxh    = psxh/sum(sum(psxh));
    psxh    = psxh(:);
    
    % intrinsic MI minus KL
    %----------------------------------------------------------------------
    Gi(j,1) = psxh'*log(psxh + eps) - ps'*log(ps) - ph'*log(ph);
    Ge(j,1) = ph'*(log(ph) - log(pH));
end
I3   =  pS'*(Ge - Gi);

return


function G  = spm_KL(A,lnpH)

% KL divergence between density over blanket states spm_softmax((Gi - Ge))
%--------------------------------------------------------------------------
pH        = spm_softmax(lnpH);
lnpH      = log(pH);
[G,Gi,Ge] = spm_G(A,lnpH);
pS        = spm_softmax((Gi - Ge));
G         = pS'*(log(pS) - log(A*pH));


