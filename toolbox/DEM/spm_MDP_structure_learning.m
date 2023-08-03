function mdp = spm_MDP_structure_learning(mdp,O,OPTIONS)
% structure learning of factorised Markov decision processes
% FORMAT mdp = spm_MDP_structure_learning(mdp,O,[OPTIONS])
%
% mdp.p - initial Dirichlet counts [1/32]
% mdp.q - precise Dirichlet counts [512]
% O     - probabilitic exemplars of paths
%
% OPTIONS.N  [0]   - suppress neuronal responses
% OPTIONS.P  [0]   - suppress plotting
% OPTIONS.G  [1]   - suppress graphics
% OPTIONS.B  [1]      - replay
%
% OPTIONS.NF [4]   - maxmium number of factors
% OPTIONS.NS [32,...] - maxmium number of states
% OPTIONS.NU [16,...] - maxmium number of paths
% OPTIONS.UB [1]      - model prior: bound on Ns
%
% mdp   - generative model: with mdp.a, mdp.b
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
% See: spm_MDP_log_evidence.m, spm_MDP_VB_update and spm_MDP_VB_sleep.m
%      spm_MDP_structure_teaching.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8439 2023-03-27 18:41:45Z guillaume $
%__________________________________________________________________________

% options for model inversion (and evaluation)
%==========================================================================
try OPTIONS.N; catch, OPTIONS.N = 0; end       % neuronal responses
try OPTIONS.P; catch, OPTIONS.P = 0; end     % suppress plotting
try OPTIONS.B; catch, OPTIONS.B = 1; end     % replay
try OPTIONS.G; catch, OPTIONS.G = 1; end       % graphics

% initialise model (mdp)
%==========================================================================
try mdp.p; catch, mdp.p  = 1/32;     end     % initial Dirichlet counts
try mdp.q; catch, mdp.q  = 512;      end       % precise Dirichlet counts

% upper bounds on states, paths and factors
%==========================================================================
NF     = 4;                                  % maximum number of factors
NS     = [32,16,16,4];                       % maximum number of states
NU     = [8,8,8,8];                          % maximum number of paths
UB     = 0;                                    % model prior: bound on Ns
UG     = 0;                                  % model prior: expected FE

% maxmium number of states and paths for each factor
%--------------------------------------------------------------------------
try NF = OPTIONS.NF; end
try NS = OPTIONS.NS; end
try NU = OPTIONS.NU; end
try UB = OPTIONS.UB; end

% inversion scheme
%--------------------------------------------------------------------------
mdp.beta     = 0;                              % supress active learning
mdp.eta      = exp(16);                        % supress forgetting
spm_evaluate = @spm_MDP_VB_XXX;              % free energy evaluation
spm_G        = @spm_MDP_MI;                  % expected free energy
spm_F        = @(mdp)sum(mdp.F);               % variational free energy


        % number of outcomes for each modality
%--------------------------------------------------------------------------
Ng    = size(O{1},1);
        for g = 1:Ng
    mdp.No(g) = numel(O{1}{g,1});
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
for t = 1:size(O,2)                            % cycle over epochs

    % Bayesian model selection for this epoch
    %======================================================================
    [Nf,Ns,Nu] = spm_MDP_size(mdp);
    mdp.T      = size(O{t},2);                 % number of outcomes
    mdp.U      = zeros(1,Nf);                 % suppress active sampling
    mdp.O      = O{t};                         % outcome probabilities

    % no expansion
    %----------------------------------------------------------------------
    hdp    = {spm_expand(mdp,0,0,0,OPTIONS)};
    hdp{1} = spm_evaluate(hdp{1},OPTIONS);
    F      = spm_F(hdp{1});

    % Bayesian model expansion: add state to last factor
    %----------------------------------------------------------------------
    if Nu(Nf) == 1 && Ns(Nf) < NS(Nf)

        hdp{end + 1} = spm_expand(mdp,0,1,0,OPTIONS);        % expand
        hdp{end}     = spm_evaluate(hdp{end},OPTIONS);       % evaluate
        F(end + 1)   = spm_F(hdp{end});                      % record
        

        % lossless compression (based on expected free energy)
        %------------------------------------------------------------------
        [d,s]     = max(hdp{1}.D{Nf});
        if F(end) == max(F) && Nf == 1 && Ns(Nf) > 1
            hdp{end} = spm_resolve(hdp{end},s);
    end

        % model priors: number of states
        %------------------------------------------------------------------
        if UB

            % Species discovery
            %--------------------------------------------------------------
            % p      = 1 - (1 - 1/Ns(Nf)/2).^sum(mdp.a{1},'all');
            % p      = 1 - exp(-8*Ns(Nf)/NS(Nf));
            p        = Ns(Nf)/NS(Nf);
            G        = spm_log((1 - p)) - spm_log(p);
            F(end)   = F(end) + G;

    end
    end

    % Bayesian model expansion: add path to last factor
    %----------------------------------------------------------------------
    if Ns(Nf) > 1 && mdp.T > 1 && Nu(Nf) < NU(Nf)
        hdp{end + 1} = spm_expand(mdp,0,0,1,OPTIONS);

        % priors over intial states
        %------------------------------------------------------------------
        hdp{end}.D{Nf} = hdp{1}.D{Nf};
        hdp{end}     = spm_evaluate(hdp{end},OPTIONS);
        F(end + 1)     = spm_F(hdp{end});

        % lossless compression (based on expected free energy)
        %------------------------------------------------------------------
        [d,u]     = max(hdp{1}.E{Nf});
        if F(end) == max(F) && Nu(Nf) > 1
            hdp{end} = spm_reduce(hdp{end},u);
        end

    end

    % Bayesian model expansion: add factor
    %----------------------------------------------------------------------
    if min(Ns) > 1 &&  mdp.T > 1 && Nf < NF
        hdp{end + 1} = spm_expand(mdp,1,0,0,OPTIONS);
        hdp{end}     = spm_evaluate(hdp{end},OPTIONS);
        F(end + 1)   = spm_F(hdp{end});
    end

    % Bayesian model selection
    %----------------------------------------------------------------------
    [m,h]       = max(F);
    FF(end + 1) = m;

    % selected parameters
    %----------------------------------------------------------------------
        mdp.a = hdp{h}.a;
        mdp.b = hdp{h}.b;



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
        spm_MDP_params(mdp)

                % free energies (ELBO)
            %--------------------------------------------------------------
                subplot(4,2,6)
                bar(F - min(F))
                axis square, title('ELBO')
                subplot(4,2,8)
                bar(FF)
                axis square, title('Successive ELBO')

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
mdp     = rmfield(mdp,'T');
try
    mdp = rmfield(mdp,'o');
end
try
    mdp = rmfield(mdp,'O');
end
try
    mdp = rmfield(mdp,'k');
end
try
    mdp = rmfield(mdp,'U');
end


return

% NOTES: 
%==========================================================================
%     % Bayesian model reduction
%     %--------------------------------------------------------------------
%     for g = 1:Ng
%         mdp.a{g} = spm_MDP_VB_prune(mdp.a{g},mdp.p,0,0,0,'SIMPLE');
%     end

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
        mdp.E{f} = sparse(1,1,1,Nu(f),1);       % initial path
        end
    end

if logical(n)

% add factor: with 2 states, because the first is shared by all factors
    %======================================================================

    % add prior transition parameters: a precise identity mapping
    % ---------------------------------------------------------------------
    mdp.b{end + 1} = eye(2,2)*mdp.q + mdp.p;
    Nf             = Nf + 1;
    Ns(Nf)         = 2;
    Nu(Nf)         = 1;

    % add likelihood parameters
    % ---------------------------------------------------------------------
    for g = 1:Ng
        a        = zeros([No(g),Ns]) + mdp.p;
        k         = numel(mdp.a{g});
        a(1:k)   = mdp.a{g};
        mdp.a{g} = a;
    end

    % starting in the second state of the new factor
    %----------------------------------------------------------------------
    mdp.D{Nf}    = sparse(2,1,1,Ns(Nf),1);       % initial state
    mdp.E{Nf}    = sparse(1,1,1,Nu(Nf),1);       % initial path


elseif logical(s)

% add latent state
    %======================================================================

    % augment priors: the first path is a precise identity mapping  
    % ---------------------------------------------------------------------
    mdp.b{Nf}(end + 1,:,:) = mdp.p;
    mdp.b{Nf}(:,end + 1,:) = mdp.p;
    mdp.b{Nf}(end,end,1)   = mdp.q + mdp.p;
    Ns(Nf)                 = Ns(Nf) + 1;

    % likelihood: and imprecise otherwise
    % ---------------------------------------------------------------------
    for g = 1:Ng
        a        = zeros([No(g),Ns]) + mdp.p;
        k         = numel(mdp.a{g});
        a(1:k)   = mdp.a{g};
        mdp.a{g} = a;
    end

    % new state is the initial state
    %----------------------------------------------------------------------
    mdp.D{Nf}     = sparse(Ns(Nf),1,1,Ns(Nf),1);
    mdp.E{Nf}     = ones(Nu(Nf),1);

elseif logical(u)

% add control (i.e., path)
    %======================================================================

    % augment priors: were the new path is imprecise  
    % ---------------------------------------------------------------------
    mdp.b{Nf}(:,:,end + 1) = mdp.p;
    Nu(Nf)                 = Nu(Nf) + 1;

    % new path is the initial path
    %----------------------------------------------------------------------
    mdp.D{Nf}     = ones(Ns(Nf),1);
    mdp.E{Nf}     = sparse(Nu(Nf),1,1,Nu(Nf),1);


else % original model

    % ELBO under most likely state of last factor
    %======================================================================
    mdp.D{Nf} = ones(Ns(Nf),1);
    mdp.E{Nf} = ones(Nu(Nf),1);

    pdp       = spm_MDP_VB_XXX(mdp,OPTIONS);
    [d,i]     = max(pdp.X{Nf}(:,1));
    [d,j]     = max(pdp.P{Nf}(:,1));

    mdp.D{Nf} = sparse(i,1,1,Ns(Nf),1);
    mdp.E{Nf} = sparse(j,1,1,Nu(Nf),1);

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


function mdp = spm_reduce(mdp,u)
% FORMAT mdp = spm_reduce(mdp,u)
% mdp  - MDP structure
% u    - paths to evaluate: default [1:Nu]
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
Nu = size(mdp.b{end},3);                        % number of hidden controls

if Nu == 1, return, end

% factors and initial Dirichlet counts
%--------------------------------------------------------------------------
try iu = [s,Nu]; catch, iu = 1:Nu; end
try p = mdp.p; catch, p  = 1/32; end  

% reduce transition matrices
%--------------------------------------------------------------------------
spm_G = @spm_MDP_MI;                            % expected free energy
G     = zeros(Nu,1);
db    = mdp.b{end}(:,:,Nu) - p;
for i = iu

    % evaluate mutual information
    %----------------------------------------------------------------------
    b         = mdp.b{end};
        
    b(:,:,Nu) = b(:,:,Nu) - db;
    b(:,:,i)  = b(:,:,i)  + db;

    % mutual information (EFE)
    %----------------------------------------------------------------------
    G(i) = spm_G(b);

end

% remove redundant paths
%--------------------------------------------------------------------------
[g,i] = max(G);
if i < Nu
    b          = mdp.b{end};
    b(:,:,i) = b(:,:,i) + b(:,:,Nu) - p;
    mdp.b{end} = b(:,:,1:(Nu - 1));
end

%  remove A and B if necessary
%--------------------------------------------------------------------------
if isfield(mdp,'A'), mdp = rmfield(mdp,'A'); end
if isfield(mdp,'B'), mdp = rmfield(mdp,'B'); end
if isfield(mdp,'d'), mdp = rmfield(mdp,'d'); end
if isfield(mdp,'e'), mdp = rmfield(mdp,'e'); end
if isfield(mdp,'k'), mdp = rmfield(mdp,'k'); end


return


function mdp = spm_resolve(mdp,s)
% FORMAT mdp = spm_resolve(mdp)
% mdp  - MDP structure
% s    - states to evaluate: default [1:Ns]
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
Ns    = size(mdp.a{1}(:,:),2);
if Ns == 1, return, end


% states and initial Dirichlet counts
%--------------------------------------------------------------------------
try is = [s,Ns]; catch, is = 1:Ns; end
try p = mdp.p; catch, p  = 1/32; end    

% reduce likelihood matrix
%--------------------------------------------------------------------------
spm_G = @spm_MDP_MI;                            % expected free energy
G     = zeros(Ns,1);
for i = is

    % evaluate reduced likelihood
    %----------------------------------------------------------------------
    for g = 1:numel(mdp.a)
        
        % evaluate mutual information
        %------------------------------------------------------------------
        da      = mdp.a{g}(:,Ns) - p;
        a       = mdp.a{g};
        a(:,Ns) = a(:,Ns) - da;
        a(:,i)  = a(:,i) + da;

        % mutual information (EFE)
        %------------------------------------------------------------------
        G(i) = G(i) + spm_G(a);

    end
end

% contract likelihood matrix if there is no loss of information
%--------------------------------------------------------------------------
[g,i] = max(G);
if i < Ns
    for g = 1:numel(mdp.a)
        mdp.a{g}(:,i) = mdp.a{g}(:,i) + mdp.a{g}(:,Ns) - p;
        mdp.a{g}      = mdp.a{g}(:,1:Ns - 1);
    end
    mdp.b{1}          = mdp.b{1}(1:end - 1,1:end - 1,:);
end

%  remove A and B if necessary
%--------------------------------------------------------------------------
if isfield(mdp,'A'), mdp = rmfield(mdp,'A'); end
if isfield(mdp,'B'), mdp = rmfield(mdp,'B'); end
if isfield(mdp,'d'), mdp = rmfield(mdp,'d'); end
if isfield(mdp,'e'), mdp = rmfield(mdp,'e'); end
if isfield(mdp,'k'), mdp = rmfield(mdp,'k'); end


return

% for comparison: reduce likelihood matrix using BMR
%--------------------------------------------------------------------------
% for i = 1:Ns
% 
%     % evaluate reduced likelihood
%     %--------------------------------------------------------------------
%     for g = 1:numel(mdp.a)
% 
%         qa = mdp.a{g}(:,end);
%         pa = qa - qa + p;
%         ra = mdp.a{g}(:,i);
% 
%         % BMR
%         %----------------------------------------------------------------
%         G(i) = G(i) - spm_MDP_log_evidence(qa,pa,ra);
% 
%     end
% 
% end
