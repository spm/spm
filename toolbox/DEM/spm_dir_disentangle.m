function [e,v,c] = spm_dir_disentangle(a,n)
% information geometry of a likelihood mapping
% FORMAT [e,v,c] = spm_dir_disentangle(a,n)
% a{g}    - Dirichlet tensor for modality g
% n       - nummber of (first) latent states to plot
%
% This routine illustrates the information geometry inherent in a
% likelihood mapping. It computes a similarity matrix based upon an
% approximation to the information length between columns of a likelihood
% tensor (as approximated by the KL divergence of the implicit categorical
% distributions). It then plots a certain number of latent states in the
% ensuing eigenspace; along with the eiegnspectrum.
%
% seee also spm_information_distance.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging

% preliminaries
%--------------------------------------------------------------------------
if nargin < 2
    n = size(a{1}(:,:),2);  % number of states to plot
end

% reduce tensor if necessary
%--------------------------------------------------------------------------
Ns    = size(a{1}(:,:,:),2:3);
for g = 1:numel(a)
    for i = 1:Ns(2)
        A{g}(:,:,i) = spm_dir_norm(a{g}(:,1:n,i));
    end
end
a     = A;
Ns    = size(a{g}(:,:,:),2:3);
N     = prod(Ns);

% information geometry – divergence : likelihood mapping
%==========================================================================

% KL divergence squared
%--------------------------------------------------------------------------
% spm_KL = @(q,p) (q'*(spm_log(q) - spm_log(p))).^2;

c     = zeros(N,N);
C     = zeros(N,N);
for g = 1:numel(a)
    q = a{g};
    for i = 1:N
        p = a{g}(:,i);
        c(:,i) = sum(q.*minus(spm_log(q),spm_log(p))).^2;
    end
    C  = c + c';
end

% similarity matrix (assuming normalised vectors on a hypersphere)
%==========================================================================
C     = real(C);
rho   = C/max(C,[],'all');                       % normalise distance
rho   = 1 - 2*rho;                               % correlation matrix

if Ns(2) > 1
    X = kron(eye(Ns(2),Ns(2)),ones(Ns(1),1));   % classification
    R = eye(N,N) - X*pinv(X);                   % residual projector
    c = rho - R*rho*R;                          % classification ssq.
else
    c = rho;
end

[e,v] = svd(spm_cov2corr(c));                   % eigenvectors
v     = diag(real(v));
[v,i] = sort(v,'descend');
e     = real(e(:,i));

% return unless called for plotting
%--------------------------------------------------------------------------
if Ns(2) == 1 || nargout
    return
end

% Display [eigen] space
%==========================================================================
spm_figure('GetWin','Information geometry'); clf;


% eigenspace
%--------------------------------------------------------------------------
subplot(2,2,2), bar(v(1:min(numel(v),16))), 
xlabel('eigenvectors'), ylabel('eigenvalues'), title('Eigenvalues')
axis square

subplot(2,2,1), imagesc(rho);
xlabel('latent states'), ylabel('latent states')
title('Correlation matrix'), axis image

e(:,1) = [];                                    % supress 1st vector
col    = get(gca,'ColorOrder');
col    = [col; spm_softmax(2*randn(3,10))'];
for n  = 1:Ns(2)

    subplot(2,1,2)
    j = (1:Ns(1)) + Ns(1)*(n - 1);
    plot3(e(j,1),e(j,2),e(j,3),'.','Color',col(n,:)), axis square, hold on
    plot3(e(j,1),e(j,2),e(j,3),'o','Color',col(n,:)), axis square, hold on
    xlabel('2nd'), ylabel('3rd'), zlabel('4th'), title('Embedding sapce')
    grid on

    text(mean(e(j,1)),mean(e(j,2)),mean(e(j,3)), ...
        num2str(n - 1),'Color',col(n,:),...
        'FontSize',24)

end
hold off


return



% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% subroutines for notes
%==========================================================================

%  objective function
%--------------------------------------------------------------------------
function F = spm_MIB(b,G,DIM)
% b    - transition matrix (unfactorised) Dirichlet parameters
% G    - grouping operator
% DIM  - tensor size (compressed)
%
% F    - Free energy (MI)
%__________________________________________________________________________

% MI = H(S) - H(S|O)
%--------------------------------------------------------------------------
F     = 0;
P     = G'*b*G;
P     = reshape(P,[DIM,DIM]);
Ndim  = numel(DIM);
for i = 1:Ndim
    for j = 1:Ndim
        if i ~= j
            p = squeeze(spm_sum(P,[i,j + Ndim]));
            F = F - spm_MDP_MI(p);
        end
    end
end
P    = squeeze(spm_sum(P,1:Ndim));
P    = P/sum(P,'all');
p    = spm_vec(P);
P    = spm_kron(spm_marginal(P));            % kron
F    = -P'*(spm_log(P) - spm_log(p));        % H(S)

return

%  objective function
%--------------------------------------------------------------------------
function F = spm_MIX(a,G,DIM)
% a{g} - vectorised  matrix form (uncompressed) Dirichlet  parameters
% G    - grouping operator
% DIM  - tensor size (compressed)
%
% F    - Free energy (MI)
%__________________________________________________________________________


% MI = H(S) - H(S|O)
%--------------------------------------------------------------------------
g    = 1;
a{g} = a{g}(:,:)*G;
a{g} = a{g}/sum(a{g},'all');
N{g} = sum(a{g},1);
P    = N{g};
p    = spm_vec(P);
P    = reshape(P,DIM);                       % grouped
P    = spm_kron(spm_marginal(P));            % kron
F    = - P'*(spm_log(P) - spm_log(p));       % H(S)

return


%% NOTES on disentanglement in Dirichlet matrices
%==========================================================================

% iillustration of  gradient descent using softmax operators on  G
%--------------------------------------------------------------------------
Ns    = [8,2];
MI    = @(a,G,Ns)spm_MIX(a,spm_softmax(G),Ns) + spm_MDP_MI(spm_softmax(G))*0;
G     = eye(size(a{1},2),size(a{1},2));
g     = linspace(-8,8,128);
F     = [];
for i = g
   F(end + 1) = MI(a,G*i,Ns);
end
subplot(2,2,2), plot(g,F), hold on

G     = randn(size(G));
F     = [];
for i = g
   F(end + 1) = MI(a,G*i,Ns);
end
subplot(2,2,2), plot(g,F), hold off


% illustrate basic tensor algebra
%--------------------------------------------------------------------------
A    = @(a) a/sum(a,'all');
H    = @(a) A(a(:))'*spm_log(A(a(:)));
n    = 4;
m    = 3;
a    = A(rand(m,n));
x{1} = A(rand(m,1));
x{2} = A(rand(n,1));
xx   = x{1}*x{2}';

disp(x{1}'*a*x{2})
disp(spm_vec(a)'*kron(x{2},x{1}))
disp(' ')
disp(H(kron(x{2},x{1})))
disp(H(x{2}) + H(x{1}))
disp(H(sum(xx,1)) + H(sum(xx,2)))


%% create a multiset tensor by replicating a Dirichlet tensor
%==========================================================================
n   = 4;
m   = 3;
a0  = eye(n,n);
a   = kron(8*rand(1,m),a0) + 1/32;
a   = a(:,randperm(size(a,2)));

% now try to recover the original set (a0) by maximising the mutual
% information MI = H(O,S) - H(O) - H(S)
%--------------------------------------------------------------------------
MI   = @(a,R)spm_MDP_MI(a*spm_softmax(R')');
R    = randn(size(a,2),size(a,2))/16;

% maximise mutual information with respect to R
%--------------------------------------------------------------------------
for i = 1:128
    [dFdR,f] = spm_diff(MI,a,R,2);
    R        = R + 64*reshape(dFdR(:),size(R));
    F(i)     = f;

    subplot(2,2,2)
    imagesc(a*spm_softmax(R')')
    subplot(2,2,1)
    plot(F), drawnow
end








