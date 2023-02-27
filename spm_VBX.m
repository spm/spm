function [Q,F] = spm_VBX(O,P,A,METHOD)
% vvariational Bayes estimate of categorical posterior over factors
% FORMAT [Q,F] = spm_VBX(O,P,A,[METHOD])
%
% O{g}    -  outcome probabilities over each of G modalities
% P{f}    -  (empirical) prior over each of F factors
% A{g}    -  likelihood tensor for modality g
%
% Q{f}    -  variational posterior for each of F factors
% F       -  (-ve)  variational free energy or ELBO
%
% This routine is a simple implementation of variational Bayes for discrete
% state space models  under a mean field approximation, in which latent
% states are partitioned into factors (and the distribution over outcomes
% is also assumed to be conditionally independent). It takes cell arrays of
% outcome probabilities,  prior probabilities over factors and a likelihood
% tensor parameterising the  likelihood of an outcome for any combination
% of latent states. The optional argument METHOD [default: exact] switches
% among number of approximate schemes:
%
% 'full'    :  a vanilla variational scheme that uses a coordinate descent
% over a small number (four) iterations
%
% 'exact'   :  a non-iterative heuristic but numerically accurate scheme
% that replaces the variational density over hidden factors with the
% marginal over the exact posterior
%
% 'sparse'  :  as for the exact scheme but suitable for sparse tensors
%
% 'marginal':  a heuristic scheme  that uses the log of the marginalised
% likelihood and log prior to estimate the lot posterior
%
% see: spm_MDP_VB_XXX.m (NOTES)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
if nargin < 4, METHOD = 'exact'; end

switch METHOD

    case 'full'

        %  (iterative) variational scheme
        %==================================================================

        % accumulate likelihoods over modalities
        %------------------------------------------------------------------
        L     = 1;
        for g = 1:numel(O)
            L = L.*spm_dot(A{g},O{g});
        end
        L     = spm_log(L);

        % preclude numerical overflow of log likelihood
        %------------------------------------------------------------------
        L     = max(L,max(L,[],'all') - 8);

        % log prior
        %------------------------------------------------------------------
        for f = 1:numel(P)
            LP{f} = spm_vec(spm_log(P{f}));
        end

        % variational iterations
        %------------------------------------------------------------------
        Q     = P;
        Z     = -Inf;
        for v = 1:4
            F     = 0;
            for f = 1:numel(P)

                % log likelihood
                %----------------------------------------------------------
                LL   = spm_vec(spm_dot(L,Q,f));

                % posterior
                %----------------------------------------------------------
                Q{f} = spm_softmax(LL + LP{f});

                % (-ve) free energy (partition coefficient)
                %----------------------------------------------------------
                F    = F + Q{f}'*(LL + LP{f} - spm_log(Q{f}));
            end

            % convergence
            %--------------------------------------------------------------
%             if F > 0
%                 warning('positive ELBO in spm_VBX')
%             end
            dF = F - Z;
            if dF < 1/128
                break
            elseif dF < 0
                warning('ELBO decreasing in spm_VBX')
            else
                Z = F;
            end
        end



    case 'exact'

        % approximation with marginals of exact posterior
        %==================================================================
        Nf    = size(A{1});
        L     = 1;
        for g = 1:numel(O)
            L = L.*spm_dot(A{g},O{g});             % likelihood over modalities
        end
        U     = L.*spm_cross(P);                   % posterior unnormalised
        Z     = sum(U,'all');                      % partition coefficient
        F     = spm_log(Z);                        % negative free energy 
        U     = reshape(U/Z,[Nf(2:end),1]);        % joint posterior
        Q     = spm_marginal(U);                   % marginal  posteriors

        % (-ve) free energy (partition coefficient)
        %------------------------------------------------------------------
        F     = 0;
        for f = 1:numel(P)
            LL = spm_vec(spm_dot(spm_log(L),Q,f));
            LP = spm_vec(spm_log(P{f}));
            F  = F + Q{f}'*(LL + LP - spm_log(Q{f}));
        end

    case 'sparse'

        % approximation with marginals suitable for sparse tensors
        %==================================================================
        Nf    = size(A{1});
        L     = 1;
        for g = 1:numel(O)
            L = L.*(O{g}'*A{g}(:,:));              % likelihood over modalities
        end
        U     = spm_vec(L).*spm_vec(spm_cross(P)); % posterior unnormalised
        Z     = sum(U,'all');                      % partition coefficient
        F     = spm_log(Z);                        % negative free energy 
        U     = reshape(U/Z,[Nf(2:end),1]);        % joint posterior
        Q     = spm_marginal(U);                   % marginal  posteriors

    case 'marginal'

        %  marginalised likelihood
        %==================================================================
        for g = 1:numel(O)
            L{g} = spm_marginal(spm_dot(A{g},O{g}));
        end

        % log prior
        %------------------------------------------------------------------
        for f = 1:numel(P)
            LP{f} = spm_vec(spm_log(P{f}));
        end

        F     = 0;
        for f = 1:numel(P)
            LL    = 0;
            for g = 1:numel(O)
                LL = LL + spm_log(L{g}{f});
            end
            Q{f} = spm_softmax(LL + LP{f});

            % (-ve) free energy (partition coefficient)
            %--------------------------------------------------------------
            F   = F + Q{f}'*(LL - spm_log(Q{f}));
        end


    otherwise
        disp('unknown method')

end

return



%% NOTES on disentanglement in Dirichlet matrices
%==========================================================================
% create a multiset likelihood mapping by replicating a Dirichlet
% likelihood

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

%  maximise mutual information with respect to R
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



%% now repeat but using MI = H(O) - H(O|S) where H(O) is constant under
% the operation of R
%--------------------------------------------------------------------------
clf, clear F

A    = @(a)a/sum(a,'all');
L    = @(a)spm_dir_norm(a);
HO   = @(a)sum(A(a),2)'*spm_log(sum(A(a),2));
HOS  = @(a)sum(L(a).*spm_log(L(a)),1)*sum(A(a),1)';

disp(HO(a*spm_softmax(R')'));
disp(HO(a));
disp(' ');
disp(MI(a,1));
disp(HOS(a) - HO(a));
disp(' ');

H    = @(a,R)HOS(a*spm_softmax(R')');
R    = randn(size(a,2),size(a,2))/16;

%  maximise mutual information with respect to R
%--------------------------------------------------------------------------
for i = 1:128
    [dFdR,f] = spm_diff(H,a,R,2);
    R        = R + 64*reshape(dFdR(:),size(R));
    F(i)     = f;
    subplot(2,2,1)
    plot(F), drawnow
end

% plot distributed Dirichlet parameters
%--------------------------------------------------------------------------
subplot(2,2,2)
imagesc(a*spm_softmax(R')')

% illustrate evaluation for very large matrices
%--------------------------------------------------------------------------
hso   = 0;
s     = sum(A(a),1);
N     = L(a);
for i = 1:size(a,2)
    hso  = hso + s(i)*N(:,i)'*spm_log(N(:,i));
end

disp(HOS(a))
disp(hso)


%% now repeat but using a multimodal/multifactorial distribution
%--------------------------------------------------------------------------
clf, clear

% specify a iikelihood tensor for three binary outcomes and two factors
%--------------------------------------------------------------------------
aa{1}(1,:,:) = [1 0;
                1 0];    % if x1
aa{2}(1,:,:) = [1 1;
                0 0];    % if x2
aa{3}(1,:,:) = [1 0;
                0 0];    % if x1 & x2

% vectorise likelihood tensor to matrix
%--------------------------------------------------------------------------
for g = 1:numel(aa)
    aa{g}(2,:,:) = 1 - aa{g}(1,:,:);
    ag{g}        = aa{g}(:,:);
end

% convert into a likelihood mapping to unique outcomes, for each state
%--------------------------------------------------------------------------
for i = 1:size(ag{1},2)  % for the total number of combinations of states
    d = 1;
    for g = 1:numel(ag)
        d = kron(d,ag{g}(:,i));
    end
    a0(:,i) = d;
end

% simulate Dirichlet counts for three repetitions of each state
%--------------------------------------------------------------------------
m   = 8*rand(1,3);                  % three repetitions
j   = randperm(size(ag{1},2)*3);    % assigned randomly to each combination
for g = 1:numel(ag)
    a{g}   = kron(m,ag{g}) + 1/32;
    a{g}   = a{g}(:,j);
end

% set up objective function H(O|S) accumulated over modalities
%--------------------------------------------------------------------------
A    = @(a)a/sum(a,'all');
L    = @(a)spm_dir_norm(a);
HOS  = @(a)sum(L(a{1}).*spm_log(L(a{1})))*sum(A(a{1}))' + ...
           sum(L(a{2}).*spm_log(L(a{2})))*sum(A(a{2}))' + ...
           sum(L(a{3}).*spm_log(L(a{3})))*sum(A(a{3}))';

% add reallocation (R)
%--------------------------------------------------------------------------
H    = @(a,R)HOS({a{1}*spm_softmax(R')',a{2}*spm_softmax(R')',a{3}*spm_softmax(R')'});
R    = randn(size(a{1},2),size(a{1},2))/32;

% maximise mutual information with respect to R
%--------------------------------------------------------------------------
for i = 1:256
    [dFdR,f] = spm_diff(H,a,R,2);
    R        = R + 32*reshape(dFdR(:),size(R));
    F(i)     = f;
    subplot(2,2,1)
    plot(F), drawnow
end

%  redistribute Dirichlet parameters and plot results
%--------------------------------------------------------------------------
Na    = numel(a);
for g = 1:Na
    subplot(2*Na,2,2*g)
    a{g} = a{g}*spm_softmax(R')'
    imagesc(a{g})
end


%%  now recover factorisation  of the state space
%==========================================================================

%  remove redundant state combinations
%--------------------------------------------------------------------------
n     = 0;
for g = 1:Na
    n = n + sum(a{g});
end
[n,i] = sort(n,'descend');
i     = i(1:4);
for g = 1:Na
    a{g} = a{g}(:,i);
end


%  redistribution operator
%--------------------------------------------------------------------------
R = eye(size(a{1},2),size(a{1},2))/2;
T = @(a,R) reshape(a*spm_softmax(R')',size(a,1),2,2);
G = @(a,R) {T(a{1},R) T(a{2},R) T(a{3},R)};

% ....


%  maximise mutual information with respect to R
%--------------------------------------------------------------------------
for i = 1:256
    [dFdR,f] = spm_diff(H,a,R,2);
    R        = R + 32*reshape(dFdR(:),size(R));
    F(i)     = f;
    subplot(2,2,1)
    plot(F), drawnow
end






