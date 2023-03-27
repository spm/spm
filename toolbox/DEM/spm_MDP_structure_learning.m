function MDP = spm_MDP_structure_learning(MDP,mdp)
% structure learning of factorised Markov decision processes
% FORMAT mdp = spm_MDP_structure_learning(MDP,[mdp])
% FORMAT mdp = spm_MDP_structure_learning(o,[mdp])
%
% MDP   - generative process or
% o     - cell array of outcomes
% mdp.p - initial Dirichlet counts [1/16]
% mdp.q - precise Dirichlet counts [512]
%
% mdp   - generative model: with mdp.a, mdp.b (and mdp.k)
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
% (i) Dynamics are conditionally independent of outcomes. This means that
% the generative model can be factorised into a likelihood mapping (A) and
% transition probabilities over latent states (B)
% 
% (ii) Latent states can be partitioned into factors, whose dynamics are
% conditionally independent 
% 
% (iii) The dynamics for each factor can be partitioned into discrete
% paths.
% 
% This leads to a generic form for any level of a hierarchical (deep)
% Markov decision process in which the likelihood mapping (for each
% modality) is a tensor whose trailing dimensions correspond to the
% dimensions of each factor. The (transition) priors are tensors, with a
% probability transition matrix for each path. In addition, the initial
% state and path of each factor is specified with D and E. With this form,
% structure learning can simply consider the addition of a latent state, a
% latent path or a new factor.
% 
% It is assumed that the first path of any factor has no dynamics and
% corresponds to an identity operator. Subsequent paths can have any form.
% Because outcomes are assumed to be generated under the first level of
% each factor, they generate the same outcome. In other words, the
% likelihood mapping is shared by the first state of every factor. In turn,
% this means that adding a factor entails adding a second state to the
% implicit first state of the new factor.
%
% If called with two arguments, the outcomes are assimilated into an
% existing generative model.
%
% See: spm_MDP_log_evidence.m, spm_MDP_VB_update and spm_MDP_VB_sleep.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8439 2023-03-27 18:41:45Z guillaume $
%__________________________________________________________________________


% options for model inversion (and evaluation)
%==========================================================================
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.P = 0;                             % suppress plotting
OPTIONS.B = 1;                             % replay
OPTIONS.G = 1;                             % suppress graphics

% generate sequence of outcomes if necessary
%==========================================================================
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
    MUP.T = T;                              % number of outcomes

    % cycle over latent factors to generate outcomes
    %----------------------------------------------------------------------
    o     = {};                             % initialise outcome cell array
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
                MUP   = spm_MDP_VB_XXX(MUP);

                % outcomes for this epoch
                %==========================================================
                o{end + 1} = MUP.o;

            end
        end
    end

else

    % MDP are outcomes
    %----------------------------------------------------------------------
    o     = MDP;
    if isnumeric(o{1})

        % number of modalities
        %------------------------------------------------------------------
        Ng    = size(o{1},1);

        % number of outcomes for each modality
        %------------------------------------------------------------------
        O     = spm_cat(o);
        for g = 1:Ng
            No(g) = max(max(O(g,:)),2);
        end

    else  

        % probabilistic outcomes
        %------------------------------------------------------------------
        Ng    = size(o{1},1);

        % number of outcomes for each modality
        %------------------------------------------------------------------
        for g = 1:Ng
            No(g) = numel(o{1}{g,1});
        end

    end

end

% initialise model (mdp)
%==========================================================================
try mdp.l; catch, mdp.l  = 1;       end      % structure learning
try mdp.p; catch, mdp.p  = exp(-3); end      % initial Dirichlet counts
try mdp.q; catch, mdp.q  = exp(16); end      % precise Dirichlet counts
try mdp.r; catch, mdp.r  = 1;       end      % repetitions

% check there are sufficient outcomes for each modality
%--------------------------------------------------------------------------
try
    [nf,ns,nu,ng,no] = spm_MDP_size(mdp);
    for g = 1:ng
        do = No(g) - no(g);
        if do > 0
            i = (1:do) + no(g);
            mdp.a{g}(i,:) = mdp.p;
        end
    end

catch

    % initialise model (one factor, one state, no dynamics)
    %----------------------------------------------------------------------
    for g = 1:Ng
        mdp.a{g} = zeros([No(g),1]) + mdp.p;  % likelihood tensor
    end
    mdp.b = {eye(1,1)*mdp.q + mdp.p};         % transition probabilities

end


%  structure learning from observations
%==========================================================================

% cycle over latent factors to generate outcomes and grow the model
%--------------------------------------------------------------------------
U     = 0;                                    % Occam's razor (nats)
FF    = [];                                   % successive ELBO
for t = 1:size(o,2)                           % cycle over epochs


    % Bayesian model selection for this epoch
    %======================================================================
    [Nf,Ns,Nu] = spm_MDP_size(mdp);

    if isnumeric(o{t})
        mdp.T     = size(o{t},2);             % number of outcomes
        mdp.o     = o{t};                     % outcomes
        OPTIONS.O = 0;                        % use o
    else
        mdp.T     = size(o{t},2);             % number of outcomes
        mdp.O     = o{t};                     % outcomes
        OPTIONS.O = 1;                        % use O
    end
    mdp.U  = zeros(1,Nf);                     % supress active sampling

    % no expansion
    %----------------------------------------------------------------------
    hdp    = {mdp};
    if mdp.l

        % learning of the likelihood and transition probabilities
        %------------------------------------------------------------------
        hdp{1} = spm_expand(mdp,1,0,0);

    else
        % learning of the likelihood mapping
        %------------------------------------------------------------------
        hdp{1} = spm_expand(mdp,0,0,0);
    end
    F      = spm_evaluate(hdp{1},OPTIONS);


    % Bayesian model expansion: add state to each factor
    %----------------------------------------------------------------------
    for f = 1:Nf
        if Nu(f) == 1 && mdp.l
            hdp{end + 1} = spm_expand(mdp,f,1,0);              % expand
            F(end + 1)   = spm_evaluate(hdp{end},OPTIONS);     % evaluate
        end
    end

    % Bayesian model expansion: add path to each factor
    %----------------------------------------------------------------------
    for f = 1:Nf
        if Ns(f) > 1 && f == Nf && mdp.l
            hdp{end + 1} = spm_expand(mdp,f,0,1);              % expand
            F(end + 1)   = spm_evaluate(hdp{end},OPTIONS);     % evaluate
        end            
    end

    % Bayesian model expansion: add factor
    %----------------------------------------------------------------------
    if min(Ns) > 1 && Nu(1) > 1 && mdp.l
        hdp{end + 1} = spm_expand(mdp,Nf + 1,0,0);
        F(end + 1)   = spm_evaluate(hdp{end},OPTIONS);
    end

    % Bayesian model selection
    %----------------------------------------------------------------------
    [m,h]       = max(F);
    FF(end + 1) = m;

    % Parameter learning under selected model
    %----------------------------------------------------------------------
    for i = 1:mdp.r

        % evidence accumulation
        %------------------------------------------------------------------
        hdp{h} = spm_MDP_VB_XXX(hdp{h},OPTIONS);

    end
    if mdp.l
        mdp.a = hdp{h}.a;
        mdp.b = hdp{h}.b;
    else
        mdp.a = hdp{h}.a;
    end

    % lossless compression
    %----------------------------------------------------------------------
    [Nf,Ns,Nu] = spm_MDP_size(mdp);
    for f = 1:Nf
        mdp = spm_reduce(mdp,f);
    end

    % graphics
    %======================================================================
    if OPTIONS.G
        
        spm_figure('getwin','Structure learning'); clf

        if mdp.l
            for f = 1:numel(mdp.b)

                % plot transistions and likelihood mapping
                %----------------------------------------------------------
                for u = 1:size(mdp.b{f},3)
                    subplot(6,6,u + (f - 1)*6)
                    imagesc(spm_dir_norm(mdp.b{f}(:,:,u)))
                    title('Transition priors')
                    axis square
                end
                subplot(2,2,3)
                for g = 1:Ng
                    a{g} = spm_dir_norm(mdp.a{g}(:,:));
                end
                A = spm_cat(a');
                if size(A,1) > 256
                    spm_spy(A);
                else
                    imagesc(A)
                end
                axis square, title('Likelihood')

                % free energies (ELBO)
                %----------------------------------------------------------
                subplot(4,2,6)
                bar(F - min(F))
                axis square, title('ELBO')
                subplot(4,2,8)
                bar(FF)
                axis square, title('Successive ELBO')
            end
        else

            % state estimation
            %--------------------------------------------------------------
            spm_MDP_VB_trial(hdp{h});

            for f = 1:numel(mdp.b)
                
                % plot likelihood mapping
                %----------------------------------------------------------
                subplot(2,2,4)
                for g = 1:Ng
                    a{g} = spm_dir_norm(mdp.a{g}(:,:));
                end
                A = spm_cat(a');
                if size(A,1) > 256
                    spm_spy(A);
                else
                    imagesc(A)
                end
                axis square, title('Likelihood')

            end

        end
        drawnow
    end
end


% transcribe Dirichlet parameters to output
%--------------------------------------------------------------------------
mdp    = spm_expand(mdp,1,0,0);   % default priors
try
    MDP.a  = mdp.a;                   % Dirichlet likelihood parameters
    MDP.b  = mdp.b;                   % Dirichlet transition parameters
    MDP.D  = mdp.D;                   % initial states
    MDP.E  = mdp.E;                   % initial control
    MDP.k  = mdp.U;                   % controllable factors
catch
    mdp.k  = mdp.U;
    MDP    = mdp;
end
return

% NOTES: 
%==========================================================================
%     % Bayesian model reduction
%     %--------------------------------------------------------------------
%     for g = 1:Ng
%         mdp.a{g} = spm_MDP_VB_prune(mdp.a{g},a0{g},0,2,0,'SIMPLE');
%     end
%
%     NB: a0 = hdp{h}.a;               % priors for BMR

%     % Contraction of likelihood mapping
%     %--------------------------------------------------------------------
%     [Nf,Ns,Nu] = spm_MDP_size(mdp);
%     if Nu(end) == 1 && Nf < 3
%         mdp = spm_resolve(mdp);
%     end


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

% preclude learning of transitions
%--------------------------------------------------------------------------
if f == 0
    for f = 1:Nf
        mdp.b{f} = (mdp.b{f} - mdp.p)*exp(16);
    end

else

    %  default priors over states and control
    %----------------------------------------------------------------------
    for i = 1:Nf
        mdp.D{i} = ones(Ns(i),1);        % initial state
        mdp.E{i} = ones(Nu(i),1);        % initial control
    end

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

    % likelihood: and imprecise otherwise
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
if isfield(mdp,'k'), mdp = rmfield(mdp,'k'); end

return


function mdp = spm_reduce(mdp,f)
% FORMAT mdp = spm_reduce(mdp,f)
% mdp  - MDP structure
% f    - factor to reduce
%
% This subroutine compresses transition tensors of Dirichlet parameters
% summarising the joint distribution over successive states and paths. This
% compression effectively minimises expected free energy by minimising
% ambiguity; namely, the entropy over successive states, conditioned upon
% paths. Because this minimisation entails moving posterior probability
% mass between paths, the transition entropy remains the same and therefore
% minimising ambiguity of paths corresponds to maximising mutual
% information. This can also be read as minimising information loss.
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
if isfield(mdp,'k'), mdp = rmfield(mdp,'k'); end

return


function mdp = spm_resolve(mdp)
% FORMAT mdp = spm_resolve(mdp)
% mdp  - MDP structure
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
%==========================================================================
% This subroutine concerns the reduction of a likelihood tensor by moving
% posterior probability mass (Dirichlet counts) and assessing  the
% resulting possterior beliefs (about parameters) in terms of expected free
% energy which, in this instance, reduces to mutual information
%
% This subroutine compresses transition tensors of Dirichlet parameters
% summarising the joint distribution over outcomes and states. This
% compression effectively minimises expected free energy by minimising
% ambiguity; namely, the average conditional entropy over outcomes over
% states. Because this minimisation entails moving posterior probability
% mass between states, the entropy over outcomes remains the same and
% therefore minimising ambiguity corresponds to maximising the mutual
% information between states and outcomes. This can also be read as
% minimising information loss.
%__________________________________________________________________________


% reduction operator (R): lossless compression
%==========================================================================
[Nf,Ns] = spm_MDP_size(mdp);

% mutual information of likelihood
%--------------------------------------------------------------------------
M0    = 0;
for g = 1:numel(mdp.a)
    M0 = M0 + spm_MDP_MI(mdp.a{g});
end

% reduce likelihood matrix
%--------------------------------------------------------------------------
p     = mdp.p;                       % prior counts
F     = zeros(1,Ns(end) - 1);
for i = 1:(Ns(end) - 1)

    % evaluate reduced likelihood
    %----------------------------------------------------------------------
    M     = 0;
    for g = 1:numel(mdp.a)
        a{g}  = mdp.a{g};
        if Nf == 1
            a{g}(:,i)   = a{g}(:,i) + a{g}(:,end) - p;
            a{g}(:,end) = p;

        elseif Nf == 2
            a{g}(:,:,i)   = a{g}(:,:,i) + a{g}(:,:,end) - p;
            a{g}(:,:,end) = p;

        end

        % mutual information
        %------------------------------------------------------------------
        M = M + spm_MDP_MI(a{g});

    end

    % score loss of mutual information
    %----------------------------------------------------------------------
    F(i) = M - M0;
end

% ccontract likelihood matrix if there is no loss of information
%--------------------------------------------------------------------------
i = find(F > 0,1);
if isempty(i), return, end

for g = 1:numel(mdp.a)
    if Nf == 1
        mdp.a{g}(:,i)   = mdp.a{g}(:,i) + mdp.a{g}(:,end) - p;
        mdp.a{g}        = mdp.a{g}(:,1:end - 1);

    elseif Nf == 2
        mdp.a{g}(:,:,i) = mdp.a{g}(:,:,i) + mdp.a{g}(:,:,end) - p;
        mdp.a{g}        = mdp.a{g}(:,:,1:end - 1);

    end
end
mdp.b{end}(i,i,:)   = mdp.b{end}(i,i,:) + mdp.b{end}(end,end,:) - p;
mdp.b{end}          = mdp.b{end}(1:end - 1,1:end - 1,:);

% remove b if there is only one state
%--------------------------------------------------------------------------
if numel(mdp.b{end}) == 1
    mdp.b = mdp.b(1:end - 1);
end

return

function F = spm_evaluate(mdp,OPTIONS)
% Augment with an additional control state
% FORMAT F = spm_evaluate(mdp,OPTIONS)
% mdp  - MDP structure
%
% This subroutine evaluates the (negative) free energy (pertaining to
% states) by inverting an epoch of outcomes.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% options
%--------------------------------------------------------------------------
if nargin < 2
    OPTIONS.B = 1;
end

% accumulate free energy (ELBO)
%--------------------------------------------------------------------------
mdp = spm_MDP_VB_XXX(mdp,OPTIONS);

% path integral of ELBO
%--------------------------------------------------------------------------
F   = mdp.F;

% mutual information of likelihood
%--------------------------------------------------------------------------
% for g = 1:numel(mdp.a)
%     E(g) = spm_MDP_MI(mdp.a{g});
% end
% F     = F + sum(E);

return






