function DEM_self_MI_a
%--------------------------------------------------------------------------
% Routine to produce graphics illustrating self relative entropy or mutual
% information. A low self mutual information induces anomalous diffusion
% and itinerancy with power law scaling (i.e., self similar dynamics). This
% example uses a fixed form (quadratic) likelihood and optimises the density
% over over hidden states to minimise self mutual information explicitly.
%
% In this example where is just one Markov blanket states and one hidden
% state to illustrate noise phase symmetry breaking as self mutual
% information decreases. The subroutines illustrate the relationship
% between self mutual information, intrinsic mutual information and
% extrinsic cost.


% set up
%--------------------------------------------------------------------------
rng('default'), clf

n     = 64;                            % number of bins
S     = (1:n)'/n;                       % domain of states
dt    = 1;                              % time step for solution
N     = 16;                             % 2^N solution

% likelihood – mapping from hidden states to sensory states - A
%--------------------------------------------------------------------------
v     = @(s) ((s*3/2).^2)/512 + 1/512;
m     = @(s) ((s - 1/2)*3/2).^2 + 1/8;
for i = 1:n
    for j = 1:n
        A(:,j) = exp(-(S - m(j/n)).^2/v(j/n));
    end
end
A     = A/diag(sum(A));                % likelihood
lnpH  = log(hanning(n))/2;             % hidden states)

% subroutine to illustrate the decomposition of (relative) self entropy
%==========================================================================
% spm_G_check(A,lnpH);

% progressively optimise mutual information w.r.t. hidden states
%==========================================================================
ni    = [1 7 8];                       % batches of iterations
for g = 1:3
    for i = 1:ni(g)
        
        % evaluate self entropy for current hidden states
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
        MIe(g) = pS'*Ge;
        MI3(g) = pS'*G ;
        H(g)   = pS'*(-log(pS));
        C(g)   = H(g) - MI3(g);
        
        % Optimise marginal w.r.t. second-order mutual information
        %------------------------------------------------------------------
        % dp   = spm_diff(@spm_MI2,A,lnpH,2);
        
        % Optimise marginal w.r.t. third-order mutual information
        %------------------------------------------------------------------
        dp     = spm_diff(@spm_MI3,A,lnpH,2);
        
        % update marginal over hidden states
        %------------------------------------------------------------------
        lnpH = log(spm_softmax(lnpH + dp(:)*8));
        
        % graphics
        %------------------------------------------------------------------
        subplot(3,2,1),     imagesc(1 - A)
        title('Likelihood','FontSize',16)
        xlabel('Hidden states'), ylabel('Blanket states')
        axis square, axis xy
        
        subplot(3,2,2),     bar([MIe;MIi;H;C;MI3]')
        title('Mutual information','FontSize',16)
        xlabel('Iteration'),ylabel('Mutual information')
        axis square, axis xy,
        legend({'extrinsic','intrinsic','entropy','conditional','self-entropy'})
        set(gca,'XTickLabel',ni)
        
        subplot(3,3,g + 3), imagesc(1 - pSxH)
        j  = sum(ni(1:(g - 1))) + i;
        title(sprintf('Iteration %i',j),'FontSize',16)
        xlabel('Hidden states'), ylabel('Blanket states')
        axis square, axis xy
        
        hold on
        tS  = spm_softmax(-G);
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

% subroutines
%==========================================================================

function [G,Gi,Ge] = spm_G(A,lnpH)
% FORMAT [G,Gi,Ge] = spm_G(A,lnpH)
% G = Ge - Gi               % self entropy (extrinsic cost - intrinsic MI)
% E[Gi] = I(H,S'|S)         % conditional MI
% E[Ge] = I(H,S)            % (second-order) MI
% E[G]  = I(H,S',S)         % (third order)  MI

% evaluate marginal over hidden states
%--------------------------------------------------------------------------
n     = size(A,2);
pH    = spm_softmax(lnpH);

% evaluate joint density and posterior
%--------------------------------------------------------------------------
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
    
    % intrinsic MI minus KL (extrinsic cost)
    %----------------------------------------------------------------------
    Gi(j,1) = psxh'*log(psxh + eps) - ps'*log(ps) - ph'*log(ph);
    Ge(j,1) = ph'*(log(ph) - log(pH));
end

% self entropy
%--------------------------------------------------------------------------
G    = Ge - Gi;

return

function MI2 = spm_MI2(A,lnpH)

% evaluate marginal over hidden states
%--------------------------------------------------------------------------
pH    = spm_softmax(lnpH);
pS    = A*pH;

% evaluate joint density and posterior
%--------------------------------------------------------------------------
pSxH  = A*diag(pH);
pSxH  = pSxH/sum(sum(pSxH));

% second-order mutual information
%--------------------------------------------------------------------------
pSxH  = pSxH(:);
MI2   = pSxH'*log(pSxH + eps) - pS'*log(pS) - pH'*log(pH);

return

function MI3 = spm_MI3(A,lnpH)

% evaluate marginal over hidden states
%--------------------------------------------------------------------------
pH    = spm_softmax(lnpH);
pS    = A*pH;

% (negative) third order (self) mutual information 
%--------------------------------------------------------------------------
MI3   = - pS'*spm_G(A,lnpH);


function [G,Gi,Ge] = spm_G_check(A,lnpH)
% FORMAT [G,Gi,Ge] = spm_G(A,lnpH)
% A    - likelihhod p(S|H)
% pnpH - log prior log(p(H))
%
% G = Ge - Gi                 % self entropy (extrinsic cost - intrinsic MI)
% E[Gi] = I(H,S'|S)           % conditional MI
% E[Ge] = I(H,S)              % (second-order) MI
% E[G]  = I(H,S',S) = I(S,S') % (third order)  MI = self entropy

% evaluate marginal over hidden states
%--------------------------------------------------------------------------
n     = size(A,2);
pH    = spm_softmax(lnpH);

% evaluate joint density and posterior
%--------------------------------------------------------------------------
pSxH  = A*diag(pH);
pSxH  = pSxH/sum(sum(pSxH));
pHS   = pSxH'/diag(sum(pSxH,2) + eps);
pS    = sum(pSxH,2);

% entropies and probabilities
%--------------------------------------------------------------------------
for j = 1:n
    ph            = pHS(:,j);
    ps            = A*ph;
    psxh          = A*diag(ph);
    psxh          = psxh/sum(sum(psxh));
    psxhxS(:,:,j) = psxh*pS(j);
    
    % intrinsic MI minus KL (extrinsic cost)
    %----------------------------------------------------------------------
    Gi(j,1) = psxh(:)'*log(psxh(:) + eps) - ps'*log(ps) - ph'*log(ph);
    Ge(j,1) = ph'*(log(ph) - log(pH));
    
    Gp(j,1) = ps'*log(ps) + ph'*log(ph) + ph'*log(ph)...
              - psxh(:)'*log(psxh(:)) - ph'*log(pH);
    Gp(j,1) = ps'*log(ps) + ph'*log(ph) - psxh(:)'*log(pSxH(:));
       
end

% self entropy and its various forms
%--------------------------------------------------------------------------
G      = Ge - Gi;

% intrinsic and extrinsic expected surprises (and prior)
%--------------------------------------------------------------------------
IsHIS  = Gi'*pS;
IHS    = Ge'*pS;
IsHS   =  G'*pS;

% underlying probability distributions
%--------------------------------------------------------------------------
psxh   = squeeze(sum(psxhxS,3));
psxS   = squeeze(sum(psxhxS,2));
pHxS   = squeeze(sum(psxhxS,1));
ps     = sum(psxh,2);

% and entropies
%--------------------------------------------------------------------------
HsHS   = -psxhxS(:)'*log(psxhxS(:));
HS     = -pS'*log(pS);
Hs     = -ps'*log(ps);
HH     = -pH'*log(pH);
HsH    = -psxh(:)'*log(psxh(:));
HsS    = -psxS(:)'*log(psxS(:));
HHS    = -pHxS(:)'*log(pHxS(:));

% equivalent formulation of self entropy in terms of entropy
%--------------------------------------------------------------------------
I3     = HsHS + Hs + HH + HS - HsH - HHS - HsS;
I3     = Hs + HS - HsS;
HsHS   = HsH + HHS - HH;

subplot(3,1,1)
bar([IsHIS IHS IsHS 0 HS (HsS - HS) I3])
str  = {'I(H,s|S)','I(H,S)','I(S,s)',' ','H(s)','H(s|S)','I(S,s)'};
set(gca,'XTickLabel',str)
title('Self entropy and information')

return


