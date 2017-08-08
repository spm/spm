function DEM_I3_and_TS
%--------------------------------------------------------------------------
% This routine is a phenomenological examination of the relationship
% between conditional mutual information (i.e.,expected intrinsic mutual
% information and the exponential divergence of trajectories as the
% Rayleigh parameter of a Lorenz attractoris increased (through a pitchfork
% bifurcation and subsequent (subcritical) Hopf bifurcation. The
% (stochastic) Lorentz system is integrated for different values of the
% relay (control) parameter. The nonequilibrium steady-state density is
% then estimated by bidding into a discrete state space; while the
% bifurcations are characterised in terms of the maximal Lyapunov exponent.
% The key thing to observe is the increase in  conditional mutual
% information following the Hopf bifurcation and implicit exponential
% divergence of trajectories. This is scored by the maximal Lyapunov
% exponent crossing zero.


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
T     = 2^18;                       % length of trajectory
W     = 32;                        % precision of intrinsic fluctuations
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
    p     = zeros(b,b,b) + 1/T;
    S     = 0;
    for i = 1:T
        
        % accumulate in bins
        %------------------------------------------------------------------
        j = round(t(i,:));
        p(j(1),j(2),j(3)) = p(j(1),j(2),j(3)) + 1;
        
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
    semilogx(P(1:k),MI(1,:),'b:',P(1:k),MI(2,:),'b',P(1:k),MI(3,:),'b-.')
    axis square xy
    title('Self mutual information','Fontsize',16)
    xlabel('Control parameter'),ylabel('MI3 (nats)'),drawnow
    
    % plot maximal Lyapunov exponent
    %----------------------------------------------------------------------
    subplot(3,2,4)
    semilogx(P(1:k),LE(1,:),'b')
    axis square xy
    title('Lyapunov exponent','Fontsize',16)
    xlabel('Control parameter'),ylabel('MI3 (nats)'),drawnow
    
end

% lines and thresholds
%--------------------------------------------------------------------------
j     = find(LE(1,:) > 0,1);
j     = P(j);
subplot(3,2,3)
hold on, plot([j j],[ 0 3],':'), hold off
hold on, plot([1 1],[ 0 3],':'), hold off
legend({'Ir','Ii','Ie'})
subplot(3,2,4)
hold on, plot([j j],[-1 1]*0.03,':'), hold off
hold on, plot([1 1],[-1 1]*0.03,':'), hold off
hold on, plot(P,zeros(size(P)),'--'), hold off


% illustrated lectures and ergodic density
%--------------------------------------------------------------------------
j     = [8 38 length(P)];
for i = 1:length(j)
    
    % exemplar trajectory (plot)
    %----------------------------------------------------------------------
    subplot(3,3,i)
    t    = tt{j(i)};
    plot(t(1:1024,2),t(1:1024,3),'k')
    axis square xy
    title('trajectory','Fontsize',16)
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


function [I,Ii,Ie] = spm_self_entropy(pSxH)
% FORMAT [I,Ii,Ie] = spm_self_entropy(pSxH)
% mutual informations
% G     = Ge - Gi
% E[Gi] = I(H,S'|S)
% E[Ge] = I(H,S)
% E[G]  = I(H,S',S) = I

% size of probability distribution array
%--------------------------------------------------------------------------
n     = size(pSxH);

% evaluate joint density and posterior
%--------------------------------------------------------------------------
pH    = squeeze(sum(pSxH,1));
pS    = sum(sum(sum(pSxH,2),3),4);
pHS   = pSxH;
for i = 1:n(1)
    pHS(i,:,:,:) = pSxH(i,:,:,:)/sum(pSxH(i,:));
end
pSH   = spm_norm(pSxH);
pS    = pS(:);
pH    = pH(:);

% entropies and probabilities
%--------------------------------------------------------------------------
for j = 1:n(1)
    ph    = squeeze(pHS(j,:,:,:));
    ps    = pSH(:,:)*ph(:);
    psxh  = spm_unnorm(pSH,ph);
    
    psxh  = psxh(:);
    ps    = ps(:);
    ph    = ph(:);
    
    % intrinsic (mutual information) and extrinsic (cost)
    %----------------------------------------------------------------------
    Gi(j,1) = psxh'*log(psxh + eps) - ps'*log(ps) - ph'*log(ph);
    Ge(j,1) = ph'*(log(ph) - log(pH));
end

% expected values
%--------------------------------------------------------------------------
I    = pS'*(Ge - Gi);
Ii   = pS'*Gi;
Ie   = pS'*Ge;

return

function A = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                S = sum(A(:,i,j,k,l),1);
                if S > 0
                    A(:,i,j,k,l) = A(:,i,j,k,l)/S;
                else
                    A(:,i,j,k,l) = 1/size(A,1);
                end
            end
        end
    end
end

function A = spm_unnorm(A,B)
% a normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                if size(B,1) > 1
                    A(:,i,j,k,l) = A(:,i,j,k,l)*B(i,j,k,l);
                else
                    A(:,i,j,k,l) = A(:,i,j,k,l)*B(1,i,j,k,l);
                end
            end
        end
    end
end

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

% plot
%--------------------------------------------------------------------------
% plot(W,S,'b.',W,X*beta,'b','LineWidth',1)
% title(sprintf('alpha = %-2.2f',beta(2)),'FontSize',16)
% ylabel('Log power'), xlabel('Log frequency')
% axis square, axis xy, spm_axis tight

