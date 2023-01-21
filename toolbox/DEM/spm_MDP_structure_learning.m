function MDP = spm_MDP_structure_learning(MDP)
% structure learning of factorised Markov decision processes
% FORMAT mdp = spm_MDP_structure_learning(MDP)
% FORMAT mdp = spm_MDP_structure_learning(o)
%
% MDP - generative process or
% o   - cell array of outcomes
%
% mdp - generative model: with mdp.a, mdp.b (and mdp.k)
%
% This routine returns a generative model in the form of an MDP, based upon
% a sequence of outcomes. If the outcomes are not supplied, then they are
% generated automatically from a generative process, specified with an MDP
% structure. The generative model learns from successive epochs of data
% generated under the first level of each factor of the process. By
% exploring different extensions to the model (using Bayesian model
% comparison) successive epochs are assimilated under a model structure
% that accommodates context sensitive dynamics. This routine makes certain
% assumptions about the basic structural form of generative models at any
% given level of a hierarchical model. These are minimal assumptions:
% 
% (i) Dynamics are conditionally independent of outcomes. This means that the
% generative model can be factorised into a likelihood mapping (A) and
% transition probabilities over latent states (B)
% 
% (ii) Latent states can be partitioned into factors, whose dynamics are
% conditionally independent 
% 
% (iii) The dynamics for each factor can be partitioned
% into discrete paths.
% 
% This leads to a generic form for any level of a
% hierarchical (deep) Markov decision process in which the likelihood
% mapping (for each modality) is a tensor whose trailing dimensions
% correspond to the dimensions of each factor. The (transition) priors are
% tensors, with a probability transition matrix for each path. In addition,
% the initial state and path of each factor is specified with D and E. With
% this form, structure learning can simply consider the addition of a
% latent state, a latent path or a new factor. 
% 
% It is assumed that the first path of any factor has no dynamics and
% corresponds to an identity operator. Subsequent paths can have any form.
% Because outcomes are assumed to be generated under the first level of
% each factor, they generate the same outcome. In other words, the
% likelihood mapping is shared by the first state of every factor. In turn,
% this means that adding a factor entails adding a second state to the
% implicit first state of the new factor.
%
% See: spm_MDP_log_evidence.m, spm_MDP_VB_update and spm_MDP_VB_sleep.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_sleep.m 7679 2019-10-24 15:54:07Z spm $
%__________________________________________________________________________


% options for model inversion (and evaluation)
%==========================================================================
OPTIONS.A = 0;                             % suppress explicit action
OPTIONS.B = 0;                             % suppress backward pass
OPTIONS.D = 0;                             % suppress backward pass
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.P = 0;                             % suppress plotting
OPTIONS.G = 1;                             % suppress graphics

% generate sequence of outcomes if nnecessary
%==========================================================================
o  = {};                                   % initialise outcome cell array
if isstruct(MDP)

    % number of outcomes, states, controls and policies
    %----------------------------------------------------------------------
    [Nf,Ns,Nu,Ng,No] = spm_MDP_size(MDP);
    T = max(Ns) + 2;

    % ensure factors with dynamics generate sequences first
    %----------------------------------------------------------------------
    disp('generating outcomes for structure learning'), disp(' ')
    if any(diff(Nu) > 0)
        warning('Please reorder B{i} with descending number of paths: size(B{i},3)')
    end

    % passive generation
    %----------------------------------------------------------------------
    MUP   = MDP;
    MUP.U = zeros(1,Nf);
    MUP.T = T;                              % number of moves

    % cycle over latent factors to generate outcomes
    %----------------------------------------------------------------------
    for i = 1:Nf                            % cycle over factors

        % cycle over paths
        %------------------------------------------------------------------
        for j = 1:Nu(i)

            % cycle over states
            %--------------------------------------------------------------
            for k = 1:Ns(i)

                % generate outcomes from initial states (and paths) of this
                % factor, under the first states (and paths) of others
                %----------------------------------------------------------
                MUP.s = ones(Nf,1); MUP.s(i) = k;
                MUP.u = ones(Nf,1); MUP.u(i) = j;
                MUP.o = [];
                MUP   = spm_MDP_VB_XXX(MUP,OPTIONS);

                % outcomes for this epoch
                %==========================================================
                o{end + 1} = MUP.o;

            end
        end
    end
else

    % number of modalities and length of epochs
    %----------------------------------------------------------------------
    [Ng,T] = size(o{1});

    % number of outcomes for each modality
    %----------------------------------------------------------------------
    O     = spm_cat(o);
    for g = 1:Ng
        No(g) = max(O(g,:));
    end

end

%  structure learning from observations
%==========================================================================

% initialise model (one factor, one state, no dynamics)
%--------------------------------------------------------------------------
mdp.T = T;                                 % number of moves
mdp.p = 1/32;                              % initial Dirichlet counts
mdp.q = 1024;                              % precise Dirichlet counts
for g = 1:Ng
    mdp.a{g} = zeros([No(g),1]) + mdp.p;   % likelihood tensor
end
mdp.b = {eye(1,1)*mdp.q + mdp.p};          % transition probabilities

% cycle over latent factors to generate outcomes and grow the model
%--------------------------------------------------------------------------
N     = 2;                                 % number of iterations
U     = 3;                                 % Occam's razor (nats)
FF    = [];                                % successive ELBO
for t = 1:numel(o)                         % cycle over epochs


    % Bayesian model selection for this epoch
    %======================================================================
    [Nf,Ns,Nu] = spm_MDP_size(mdp);
    mdp.o  = o{t};
    mdp.U  = zeros(1,Nf);

    % no expansion
    %----------------------------------------------------------------------
    hdp    = {mdp};
    hdp{1} = spm_expand(mdp,1,0,0);
    F      = spm_evaluate(hdp{1},N) + U;

    % Bayesian model expansion: add state to each factor
    %----------------------------------------------------------------------
    for f = 1:Nf
        if Nu(f) < 2
            hdp{end + 1} = spm_expand(mdp,f,1,0);      % expand
            F(end + 1)   = spm_evaluate(hdp{end},N);   % evaluate
        end
    end

    % Bayesian model expansion: add path to each factor
    %----------------------------------------------------------------------
    for f = 1:Nf
        if Ns(f) > 1
            hdp{end + 1} = spm_expand(mdp,f,0,1);      % expand
            F(end + 1)   = spm_evaluate(hdp{end},N);   % evaluate
        end
    end

    % Bayesian model expansion: add factor
    %----------------------------------------------------------------------
    if min(Ns) > 1
        if FF(end) - max(F) > mdp.T*U
            hdp{end + 1} = spm_expand(mdp,Nf + 1,0,0);
            F(end + 1)   = spm_evaluate(hdp{end},N);
        end
    end

    % Bayesian model selection
    %----------------------------------------------------------------------
    [m,h]       = max(F);
    FF(end + 1) = m;


    % Parameter learning under selected model
    %----------------------------------------------------------------------
    for i = 1:N
        hdp{h} = spm_MDP_VB_XXX(hdp{h},OPTIONS);
    end
    mdp.a = hdp{h}.a;
    mdp.b = hdp{h}.b;

    % lossless compression
    %----------------------------------------------------------------------
    Nf = spm_MDP_size(mdp);
    for f = 1:Nf
        mdp = spm_reduce(mdp,f);
    end

    % graphics
    %======================================================================
    if OPTIONS.G
        spm_figure('getwin','Structure learning'); clf

        for f = 1:numel(mdp.b)
            for u = 1:size(mdp.b{f},3)
                subplot(6,6,u + (f - 1)*6)
                imagesc(spm_dir_norm(mdp.b{f}(:,:,u)))
                title('Transition priors')
                axis square
            end
            subplot(2,2,3)
            for g = 1:Ng, a{g} = mdp.a{g}(:,:); end
            imagesc(spm_dir_norm(spm_cat(a')))
            axis square, title('Likelihood')
            subplot(4,2,6)
            bar(F - min(F))
            axis square, title('ELBO')
            subplot(4,2,8)
            bar(FF)
            axis square, title('Successive ELBO')
        end
        drawnow
    end
end


% transcribe Dirichlet parameters to output
%--------------------------------------------------------------------------
mdp    = spm_expand(mdp,1,0,0);   % default priors
MDP.a  = mdp.a;                   % Dirichlet likelihood parameters
MDP.b  = mdp.b;                   % Dirichlet transition parameters
MDP.D  = mdp.D;                   % initial states
MDP.E  = mdp.E;                   % initial control
MDP.k  = mdp.U;                   % controllable factors

return


function mdp = spm_expand(mdp,f,s,u)
% Augment with an additional path
% FORMAT mdp = spm_expand(mdp)
% mdp  - MDP structure
% f    - factor to augment or add [if f > Nf] 
% s    - if logical(s) then add state
% u    - if logical(u) then add path
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% size of factor
%--------------------------------------------------------------------------
[Nf,Ns,Nu,Ng,No] = spm_MDP_size(mdp);

%  default priors over states and control
%--------------------------------------------------------------------------
for i = 1:Nf
    mdp.D{i} = ones(Ns(i),1);        % initial state
    mdp.E{i} = ones(Nu(i),1);        % initial control
end

% add factor: with 2 states, because the first is shared by all factors
%--------------------------------------------------------------------------
if f > Nf

    % add prior transition parameters: a precise identity mapping
    % ---------------------------------------------------------------------
    mdp.b{end + 1} = eye(2,2)*mdp.q + mdp.p;
    Nf             = Nf + 1;
    Ns(Nf)         = 2;
    Nu(Nf)         = 1;

    % add likelihood parameters
    % ---------------------------------------------------------------------
    for g = 1:Ng
        a{g}      = zeros([No(g),Ns]) + mdp.p;
        k         = numel(mdp.a{g});
        a{g}(1:k) = mdp.a{g};
        mdp.a{g}  = a{g};
    end

    % priors over initial states and control (first state and control)
    %----------------------------------------------------------------------
    for i = 1:Nf
        mdp.D{i} = sparse(1,1,1,Ns(i),1);        % initial state
        mdp.E{i} = sparse(1,1,1,Nu(i),1);        % initial control
    end

    % starting in the second state of the new factor
    %----------------------------------------------------------------------
    mdp.D{Nf}    = sparse(2,1,1,Ns(Nf),1);       % initial state
    mdp.E{Nf}    = sparse(1,1,1,Nu(Nf),1);       % initial control

    [Nf,Ns,Nu,Ng,No] = spm_MDP_size(mdp);

end

% add latent state
%--------------------------------------------------------------------------
if logical(s)

    % augment priors: the first path is a precise identity mapping  
    % ---------------------------------------------------------------------
    mdp.b{f}(end + 1,:,:) = mdp.p;
    mdp.b{f}(:,end + 1,:) = mdp.p;
    mdp.b{f}(end,end,1)   = mdp.q + mdp.p;
    Ns(f)                 = Ns(f) + 1;

    % and likelihood: and imprecise otherwise
    % ---------------------------------------------------------------------
    for g = 1:Ng
        a{g}      = zeros([No(g),Ns]) + mdp.p;
        k         = numel(mdp.a{g});
        a{g}(1:k) = mdp.a{g};
        mdp.a{g}  = a{g};
    end

    % priors over initial states and control (first state and control)
    %----------------------------------------------------------------------
    for i = 1:Nf
        mdp.D{i} = sparse(1,1,1,Ns(i),1);        % initial state
        mdp.E{i} = sparse(1,1,1,Nu(i),1);        % initial control
    end

    % new state is the initial state
    %----------------------------------------------------------------------
    mdp.D{f}     = sparse(Ns(f),1,1,Ns(f),1);
    mdp.E{f}     = ones(Nu(f),1);

end

% add control (i.e., path)
%--------------------------------------------------------------------------
if logical(u)

    % augment priors: were the new path is imprecise  
    % ---------------------------------------------------------------------
    mdp.b{f}(:,:,end + 1) = mdp.p;
    Nu(f)                 = Nu(f) + 1;

    % priors over initial states and control (first state and control)
    %----------------------------------------------------------------------
    for i = 1:Nf
        mdp.D{i} = sparse(1,1,1,Ns(i),1);        % initial state
        mdp.E{i} = sparse(1,1,1,Nu(i),1);        % initial control
    end

    % new path is the initial path
    %----------------------------------------------------------------------
    mdp.D{f}     = ones(Ns(f),1);
    mdp.E{f}     = sparse(Nu(f),1,1,Nu(f),1);
end

% priors over controllable factors
%--------------------------------------------------------------------------
mdp.U = zeros(1,Nf);                             % actionable controls

%  remove A and B if necessary
%--------------------------------------------------------------------------
if isfield(mdp,'A'), mdp = rmfield(mdp,'A'); end
if isfield(mdp,'B'), mdp = rmfield(mdp,'B'); end

if isfield(mdp,'d'), mdp = rmfield(mdp,'d'); end
if isfield(mdp,'e'), mdp = rmfield(mdp,'e'); end

return


function mdp = spm_reduce(mdp,f)
% FORMAT mdp = spm_reduce(mdp,f)
% mdp  - MDP structure
% f    - factor to reduce
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging


% reduction operator (R): lossless compression
%==========================================================================
Nu       = size(mdp.b{f},3);         % number of hidden controls
mdp.E{f} = ones(Nu,1)/Nu;

if Nu == 1, return, end

% reduce transition matrices
%--------------------------------------------------------------------------
p     = mdp.p;                       % prior counts
b     = mdp.b{f};                    % factor to reduce
E     = mdp.E{f};                    % initial control states
M0    = spm_MDP_MI(b);
for i = 1:Nu
    for j = (i + 1):Nu

        % evaluate reduced likelihood
        %------------------------------------------------------------------
        b(:,:,i) = b(:,:,i) + b(:,:,j) - p;
        b(:,:,j) = p;
        
        E(i) = E(i) + E(j);
        E(j) = 0;

        %  accept reduction if no loss of mutual information
        %------------------------------------------------------------------
        M  = spm_MDP_MI(b);
        if M >= M0
            mdp.b{f} = b;
            mdp.E{f} = E;
            M0       = spm_MDP_MI(b);
        else
            b = mdp.b{f};
            E = mdp.E{f};
        end
    end
end

% remove redundant paths
%--------------------------------------------------------------------------
i        = find(mdp.E{f});
mdp.b{f} = mdp.b{f}(:,:,i);


%  remove A and B if necessary
%--------------------------------------------------------------------------
if isfield(mdp,'A'), mdp = rmfield(mdp,'A'); end
if isfield(mdp,'B'), mdp = rmfield(mdp,'B'); end

if isfield(mdp,'d'), mdp = rmfield(mdp,'d'); end
if isfield(mdp,'e'), mdp = rmfield(mdp,'e'); end

return


function F = spm_evaluate(mdp,N)
% Augment with an additional control state
% FORMAT mdp = spm_evaluate(mdp,N)
% mdp  - MDP structure
% N    - number of iterations
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% suppress expplicit action and backward pass
%--------------------------------------------------------------------------
OPTIONS.A = 0;
OPTIONS.B = 0;
OPTIONS.N = 0;

% accumulate free energy (states and parameters)
%--------------------------------------------------------------------------
F     = 0;
for t = 1:N
    mdp = spm_MDP_VB_XXX(mdp,OPTIONS);
    F   = F + sum(mdp.Fa) + sum(mdp.Fb) + sum(mdp.F);
end


return



% NOTES
%==========================================================================
%  These notes concern reduction of a likelihood tensor by moving
%  posterior probability mass (Dirichlet counts) and assessing  the
%  resulting possterior beliefs (about parameters) in terms of expected
%  free energy which, in this instance, reduces to mutual information
%__________________________________________________________________________


% reduction operator (R): lossless compression
%==========================================================================
Ns    = size(mdp.b{f},1);            % number of hidden states
Nu    = size(mdp.b{f},3);            % number of hidden controls

% mutual information of likelihood
%--------------------------------------------------------------------------
M0    = 0;
for g = 1:numel(mdp.a)
    a{g} = mdp.a{g}(:,:,f);
    M0   = M0 + spm_MDP_MI(a{g});    % mutual information
end

% reduce likelihood matrix
%--------------------------------------------------------------------------
p     = mdp.p;                       % prior counts
F     = zeros(Ns,Ns);
for i = 1:Ns
    for j = (i + 1):Ns

        % evaluate reduced likelihood
        %------------------------------------------------------------------
        M     = 0;
        for g = 1:numel(mdp.a)
            a{g}      = mdp.a{g}(:,:,f);
            a{g}(:,i) = a{g}(:,i) + a{g}(:,j) - p;
            a{g}(:,j) = p;

            % mutual information
            %------------------------------------------------------------------
            M = M + spm_MDP_MI(a{g});

        end

        % score loss of mutual information
        %------------------------------------------------------------------
        F(i,j) = M - M0;
    end
end


% reduce likelihood matrix: an alternative formulation
%--------------------------------------------------------------------------
p     = mdp.p;                      % prior counts
R     = eye(Ns,Ns);
F     = zeros(Ns,Ns);
for i = 1:Ns
    for j = (i + 1):Ns

        % evaluate reduced likelihood
        %------------------------------------------------------------------
        r      = R;
        r(:,i) = r(:,i) + r(:,j);
        r(:,j) = 0;

        M     = 0;
        for g = 1:numel(a)
            M = M + spm_MDP_MI((a{g} - p)*r + p);
        end

        % score loss of mutual information
        %------------------------------------------------------------------
        F(i,j) = M - M0;

    end
end

% lossless compression
%--------------------------------------------------------------------------
for i = 1:Ns
    for j = find(F(i,:) > 0)
        R(:,i) = R(:,i) + R(:,j);
        R(:,j) = 0;
    end
end

% apply reduction operator to likelihood and transition matrices
%--------------------------------------------------------------------------
for g = 1:numel(g)
    mdp.a{g}(:,:,f) = (a{g} - p)*R + p;
end
mdp.b{f}(:,:,Nu) = R'*(mdp.b{f}(:,:,Nu) - p)*R + p;





