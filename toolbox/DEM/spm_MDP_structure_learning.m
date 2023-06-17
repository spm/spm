function mdp = spm_MDP_structure_learning(mdp,o,OPTIONS)
% structure learning of factorised Markov decision processes
% FORMAT mdp = spm_MDP_structure_learning(mdp,O,[OPTIONS])
%
% mdp.p - initial Dirichlet counts [exp( 0)]
% mdp.q - precise Dirichlet counts [exp(16)]
% O     - generative process or exemplars [o{.} or O{{.}}]
%
% OPTIONS.N  [0]   - suppress neuronal responses
% OPTIONS.P  [0]   - suppress plotting
% OPTIONS.B  [1]   - replay
% OPTIONS.G  [1]   - suppress graphics
%
% OPTIONS.NF [4]   - maxmium number of factors
% OPTIONS.NS [32,...] - maxmium number of states
% OPTIONS.NU [16,...] - maxmium number of paths
% OPTIONS.UB [1]      - model prior: bound on Ns
% OPTIONS.UG [0]      - model prior: expected FE
%
%
% mdp   - generative model: with mdp.a, mdp.b (and mdp.k)
%
% This routine returns a generative model (mdp) in the form of an MDP,
% given a sequence of outcomes. If the outcomes are not available, then
% they can be generated automatically from a generative process, specified
% with an MDP structure (see spm_MDP_structure_teaching). The generative
% model learns from successive epochs of data generated under the first
% level of each factor of the process. By exploring different extensions to
% the model (using Bayesian model comparison) successive epochs are
% assimilated under a model structure that accommodates paths. This routine
% makes certain assumptions about the structural form of generative models
% at any given level of a hierarchical model. These are minimal
% assumptions:
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
% state and path of each factor are specified with D and E. With this form,
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
% If outcomes are specified in terms of an MDP structure, the generating
% tensors are saved in the output for subsequent simulations
%
% If no outcomes are specified,  modeel priors (b) 
%
% See: spm_MDP_log_evidence.m, spm_MDP_VB_update and spm_MDP_VB_sleep.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8439 2023-03-27 18:41:45Z guillaume $
%__________________________________________________________________________

% options for model inversion (and evaluation)
%==========================================================================
try OPTIONS.N; catch, OPTIONS.N = 0; end     % suppress neuronal responses
try OPTIONS.P; catch, OPTIONS.P = 0; end     % suppress plotting
try OPTIONS.B; catch, OPTIONS.B = 1; end     % replay
try OPTIONS.G; catch, OPTIONS.G = 1; end     % suppress graphics

% initialise model (mdp)
%==========================================================================
try mdp.p; catch, mdp.p  = 1/32;     end     % initial Dirichlet counts
try mdp.q; catch, mdp.q  = exp(16); end      % precise Dirichlet counts
try mdp.w; catch, mdp.w  = 0;        end     % Occam's window

% upper bounds on states, paths and factors
%==========================================================================
NF     = 4;                                  % maximum number of factors
NS     = [32,16,16,4];                       % maximum number of states
NU     = [8,8,8,8];                          % maximum number of paths
UB     = 1;                                  % model prior: bound on Ns
UG     = 0;                                  % model prior: expected FE

% maxmium number of states and paths for each factor
%--------------------------------------------------------------------------
try NF = OPTIONS.NF; end
try NS = OPTIONS.NS; end
try NU = OPTIONS.NU; end
try UB = OPTIONS.UB; end
try UG = OPTIONS.UG; end

% precision of expected free energy priors
%--------------------------------------------------------------------------
    try, zeta = mdp.zeta; catch, zeta = 1; end

% inversion scheme
%--------------------------------------------------------------------------
spm_evaluate = @spm_MDP_VB_XXX;              % free energy evaluation
spm_G        = @spm_MDP_MI;                  % expected free energy

        % number of outcomes for each modality
%--------------------------------------------------------------------------
Ng = size(o{1},1);
    if ~isfield(mdp,'No')

        if isnumeric(o{1})

            % index outcomes
            %--------------------------------------------------------------
            oo    = spm_cat(o);
        for g = 1:Ng
                mdp.No(g) = max(max(oo(g,:)),2);
        end

    else  

        % probabilistic outcomes
            %--------------------------------------------------------------
        for g = 1:Ng
                mdp.No(g) = numel(o{1}{g,1});
            end
        end
    end

% specify single state prior (b) if necessary
%--------------------------------------------------------------------------
if ~isfield(mdp,'b')
    mdp.b{1} = mdp.q + mdp.p;
    end

% specify likelihoods (a) if necessary
%--------------------------------------------------------------------------
if ~isfield(mdp,'a')
    for g = 1:numel(mdp.No)
        mdp.a{g} = zeros([mdp.No(g),1]) + mdp.p;
    end
end


%  structure learning from observations
%==========================================================================

% cycle over latent factors to generate outcomes and grow the model
%--------------------------------------------------------------------------
FF    = [];                                   % successive ELBO
for t = 1:size(o,2)                           % cycle over epochs

    % Bayesian model selection for this epoch
    %======================================================================
    [Nf,Ns,Nu] = spm_MDP_size(mdp);
    mdp.T      = size(o{t},2);                % number of outcomes
    mdp.U      = zeros(1,Nf);                 % suppress active sampling

    if isnumeric(o{t})
        mdp.o     = o{t};                     % outcomes
        OPTIONS.O = 0;                        % use o
    else
        mdp.O     = o{t};                     % outcomes
        OPTIONS.O = 1;                        % use O
    end

    % no expansion
    %----------------------------------------------------------------------
    hdp    = {spm_expand(mdp,0,0,0,OPTIONS)};
    hdp{1} = spm_evaluate(hdp{1},OPTIONS);
    F      = sum(hdp{1}.F);

    % Bayesian model expansion: add state to last factor
    %----------------------------------------------------------------------
    if Nu(Nf) == 1 && Ns(Nf) <= NS(Nf)
        hdp{end + 1} = spm_expand(mdp,0,1,0,OPTIONS);        % expand
        hdp{end}     = spm_evaluate(hdp{end},OPTIONS);       % evaluate
        F(end + 1)   = sum(hdp{end}.F);                      % record
        
        % model priors: expected free energy
        %------------------------------------------------------------------
        if UG
            G      = spm_G(hdp{end}.a) - spm_G(hdp{1}.a);
            F(end) = F(end) + zeta*G;
    end

        % model priors: (upper bound on) number of states
        %------------------------------------------------------------------
        if UB
            G      = spm_log(1 - Ns(Nf)/NS(Nf));
            F(end) = F(end) + zeta*G;
        end            

    end

    % Bayesian model expansion: add path to last factor
    %----------------------------------------------------------------------
    if Ns(Nf) > 1 && mdp.T > 1 && Nu(Nf) <= NU(Nf)
        hdp{end + 1} = spm_expand(mdp,0,0,1,OPTIONS);
        hdp{end}     = spm_evaluate(hdp{end},OPTIONS);
        F(end + 1)   = sum(hdp{end}.F);

    end

    % Bayesian model expansion: add factor
    %----------------------------------------------------------------------
    if min(Ns) > 1 &&  mdp.T > 1 && Nf <= NF
        hdp{end + 1} = spm_expand(mdp,1,0,0,OPTIONS);
        hdp{end}     = spm_evaluate(hdp{end},OPTIONS);
        F(end + 1)   = sum(hdp{end}.F);
    end

    % Bayesian model selection
    %----------------------------------------------------------------------
    [m,h]       = max(F);
    FF(end + 1) = m;

    % Parameter learning under selected model
    %----------------------------------------------------------------------
        mdp.a = hdp{h}.a;
        mdp.b = hdp{h}.b;

    % lossless compression
    %----------------------------------------------------------------------
    [Nf,Ns,Nu] = spm_MDP_size(mdp);
    for f = 1:Nf
        mdp = spm_reduce(mdp,f);
    end

    % graphics
    %======================================================================
    if OPTIONS.G
        
        % state estimation
        %------------------------------------------------------------------
        if OPTIONS.G == 2
            spm_figure('getwin','Inference'); clf
            spm_MDP_VB_trial(hdp{h});
        end

        spm_figure('getwin','Structure learning'); clf
            for f = 1:numel(mdp.b)

                % plot transistions and likelihood mapping
            %--------------------------------------------------------------
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
            %--------------------------------------------------------------
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

% default priors
%--------------------------------------------------------------------------
[Nf,Ns,Nu] = spm_MDP_size(mdp);
for f = 1:Nf
    mdp.D{f} = ones(Ns(f),1);              % initial state
    mdp.E{f} = ones(Nu(f),1);              % initial control
end
mdp.k   = zeros(1,Nf);                      % suppress active sampling
mdp     = rmfield(mdp,'T');
try
    mdp = rmfield(mdp,'o');
end
try
    mdp = rmfield(mdp,'O');
end


return

% NOTES: 
%==========================================================================
%     % Contraction of likelihood mapping
%     %--------------------------------------------------------------------
%     [Nf,Ns,Nu] = spm_MDP_size(mdp);
%     if Ns(1) > 4 && Nu(1) == 1 && Nf == 1 && h == 2
%         mdp = spm_resolve(mdp);
%     end

%     % Bayesian model reduction
%     %--------------------------------------------------------------------
%     for g = 1:Ng
%         mdp.a{g} = spm_MDP_VB_prune(mdp.a{g},a0{g},0,0,0,'MI');
%     end
%
%     NB: a0 = hdp{h}.a;               % priors for BMR


function mdp = spm_expand(mdp,n,s,u,OPTIONS)
% Augment with an additional path
% FORMAT mdp = spm_expand(mdp)
% mdp  - MDP structure
% n    - if logical(f) then add factor 
% s    - if logical(s) then add state
% u    - if logical(u) then add path
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

OPTIONS.N = 0; % suppress neuronal responses
OPTIONS.P = 0; % suppress plotting
OPTIONS.G = 0; % suppress graphics
OPTIONS.B = 1; % replay

% size of factors
%--------------------------------------------------------------------------
[Nf,Ns,Nu,Ng,No] = spm_MDP_size(mdp);

    % priors over initial states and control (first state and control)
%--------------------------------------------------------------------------
for f = 1:Nf
        try
        mdp.D{f};                                % pre-specifed
        catch
        mdp.D{f} = sparse(1,1,1,Ns(f),1);        % initial state
        mdp.E{f} = sparse(1,1,1,Nu(f),1);        % initial control
        end
    end

% ELBO under most likely state of last factor
%--------------------------------------------------------------------------
if ~logical(n) && ~logical(s) && ~logical(u)

    pdp       = spm_MDP_VB_XXX(mdp,OPTIONS);
    [d,i]     = max(pdp.X{Nf}(:,1));

    mdp.D{Nf} = sparse(i(1),1,1,Ns(Nf),1);
    mdp.E{Nf} = ones(Nu(f),1);

end

% add factor: with 2 states, because the first is shared by all factors
%--------------------------------------------------------------------------
if logical(n)

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
    mdp.b{Nf}(end + 1,:,:) = mdp.p;
    mdp.b{Nf}(:,end + 1,:) = mdp.p;
    mdp.b{Nf}(end,end,1)   = mdp.q + mdp.p;
    Ns(Nf)                 = Ns(Nf) + 1;

    % likelihood: and imprecise otherwise
    % ---------------------------------------------------------------------
    for g = 1:Ng
        a{g}      = zeros([No(g),Ns]) + mdp.p;
        k         = numel(mdp.a{g});
        a{g}(1:k) = mdp.a{g};
        mdp.a{g}  = a{g};
    end

    % new state is the initial state
    %----------------------------------------------------------------------
    mdp.D{Nf}     = sparse(Ns(Nf),1,1,Ns(Nf),1);
    mdp.E{Nf}     = ones(Nu(Nf),1);

end

% add control (i.e., path)
%--------------------------------------------------------------------------
if logical(u)

    % augment priors: were the new path is imprecise  
    % ---------------------------------------------------------------------
    mdp.b{Nf}(:,:,end + 1) = mdp.p;
    Nu(Nf)                 = Nu(Nf) + 1;

    % new path is the initial path
    %----------------------------------------------------------------------
    mdp.D{Nf}     = ones(Ns(Nf),1);
    mdp.E{Nf}     = sparse(Nu(Nf),1,1,Nu(Nf),1);
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

if Nu == 1, return, end

% reduce transition matrices
%--------------------------------------------------------------------------
p     = mdp.p;                       % prior counts
b     = mdp.b{f};                    % factor to reduce
EE    = ones(Nu,1);                  % initial control states
E     = EE;
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
            EE       = E;
            M0       = spm_MDP_MI(b);
        else
            b = mdp.b{f};
            E = EE;
        end
    end
end

% remove redundant paths
%--------------------------------------------------------------------------
i        = find(EE);
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
for i = 1:(Ns(end) - 1)

    % evaluate reduced likelihood
    %----------------------------------------------------------------------
    M     = 0;
    for g = 1:numel(mdp.a)
        
        a{g}  = mdp.a{g};
            a{g}(:,i)   = a{g}(:,i) + a{g}(:,end) - p;
            a{g}(:,end) = p;

        % mutual information
        %------------------------------------------------------------------
        M = M + spm_MDP_MI(a{g});

    end

    % score loss of mutual information
    %----------------------------------------------------------------------
    G(i) = M - M0;

end

% for comparison: reduce likelihood matrix using BMR
%--------------------------------------------------------------------------
% F     = zeros((Ns(end) - 1),1);
% for i = 1:(Ns(end) - 1)
% 
%     % evaluate reduced likelihood
%     %----------------------------------------------------------------------
%     for g = 1:numel(mdp.a)
% 
%         qa = mdp.a{g}(:,end);
%         pa = qa - qa + p;
%         ra = mdp.a{g}(:,i);
% 
%         % score loss of mutual information
%         %----------------------------------------------------------------------
%         F(i) = F(i) - spm_MDP_log_evidence(qa,pa,ra);
% 
%     end
% 
% end


% contract likelihood matrix if there is no loss of information
%--------------------------------------------------------------------------
[G,i] = max(G); i = i(1);

if G < 0, return, end

for g = 1:numel(mdp.a)
        mdp.a{g}(:,i)   = mdp.a{g}(:,i) + mdp.a{g}(:,end) - p;
        mdp.a{g}        = mdp.a{g}(:,1:end - 1);
end
mdp.b{end}(i,i,:)   = mdp.b{end}(i,i,:) + mdp.b{end}(end,end,:) - p;
mdp.b{end}          = mdp.b{end}(1:end - 1,1:end - 1,:);

return

