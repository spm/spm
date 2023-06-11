function mdp = DEM_syntax
% Demo of active inference and structure learning (language)
%__________________________________________________________________________
%
% This routine uses a Markov decision process formulation of active
% inference to illustrate structure learning. Structure learning here is
% read as optimising the structure of a generative model that, crucially,
% includes dynamics. This foregrounds the sequential order and temporal
% scheduling of various updates. In this example, we start with a simple
% problem in which one or more objects can be removed around in a
% two-dimensional space. The kind of structure learning considered here can
% be likened to nonparametric Bayes; namely, the addition of a model
% component if licensed in terms of model evidence or marginal likelihood.
% Procedurally, this means that with each new stimulus (sequence of
% observations) various models are compared that entail the addition of a
% new latent state, path or factor. If the ELBO (i.e., negative variational
% free energy) increases the addition is accepted but not otherwise. The
% training sequences are carefully constructed to reflect the ordinal
% structure of observations. In other words, structure learning is
% predicated on both the content and dynamics of the generative process.
% 
% This demonstration calls a belief propagation scheme with factorisation
% of latent states into factors. Furthermore, the likelihood mapping is
% factorised into conditionally independent outcome modalities. This means
% that the size of any requisite tensor for belief updating is upper
% bounded by the factorisation or mean field approximation). This mitigates
% the van Neumann bottleneck; leading to increased efficiency at all three
% levels of optimisation (inference, learning and model selection).
% 
% A key aspect of this demonstration routine is that it deals with discrete
% state space generative models (and observations). This means that belief
% propagation and updating can be implemented using linear (tensor)
% operators; without worrying about nonlinearities of the sort found in
% continuous state space models.
%
%_________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_surveillance.m 8313 2022-09-30 18:33:43Z karl $


%% set up and preliminaries
%==========================================================================
learn  = @(o,mdp)spm_MDP_structure_learning(spm_word2ind(o),mdp);

% syntax
%--------------------------------------------------------------------------
syntax  = {'subject','verb','object','.'};
subject = {'the boy','the girl','the cat'};


% Illustrate structure learning via assimilation of epochs of observations
%==========================================================================
mdp.r = 16;
for i = 1:numel(syntax)
    mdp = learn(syntax(i),mdp);
end

mdp = learn(syntax,mdp);

for i = 1:numel(subject)
    mdp = learn(subject(i),mdp);
end



learn({'the boy','verb','object','.'},mdp);

return





[Nf,Ns,Nu] = spm_MDP_size(mdp);

for i = 1:Nu(1)
    for j = 1:Ns(2)
        o        = syntax{i};
        k        = find(ismember(o,'object'));
        o(k)     = object(j);

        mdp.D{1} = full(sparse(1,1,1,Ns(1),1));
        mdp.D{2} = full(sparse(j,1,1,Ns(2),1));
        mdp.E{1} = full(sparse(i,1,1,Nu(1),1));
        mdp.E{2} = full(sparse(1,1,1,Nu(2),1));

        mdp      = learn(o,mdp);

    end
end
mdp.l = 1;

learn({'subject','and','object','.'},mdp);

return

learn({{'subject','and','the boy','.'},{'subject','and','the girl','.'}},mdp)

verb  = {'verb','sat','fell','slipped'};
for i = 1:numel(verb)
    mdp = learn(verb(i),mdp);
end


return

mdp   = learn({'subject','and','the boy','.'},mdp);
mdp   = learn({'the kite'},mdp);
mdp   = learn({'the ball'},mdp);

MDP.T = 4;
MDP.A = mdp.a;
MDP.B = mdp.b;
MDP.D = mdp.D;
MDP.E = mdp.E;
MDP.U = mdp.k;
MDP.s = [1,2]';
MDP.u = [3,1]';

for i = 1:8
    PDP = spm_MDP_VB_XXX(MDP);epoch(PDP)
    disp(spm_ind2word(PDP.o))
end

return



% specify a trial, with initial conditions, s
%--------------------------------------------------------------------------
MDP   = mdp;
MDP.s = [1,1,2];
MDP   = spm_MDP_VB_XXX(MDP);


% illustrate behavioural responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP(1));

% illustrate physiological responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP(1),[],3);

% illustrate action and perceptual inference
%--------------------------------------------------------------------------
spm_figure('GetWin','optimal behaviour'); clf
spm_report(MDP)



function [o,O] = spm_word2ind(word)
% FORMAT [o,O] = spm_word2ind(word)
% word  -   word string
% r     -   number of repetitions
%
% o = matrix
% O - cell [' ' are treated as impecise outcomes]
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging


% create array of observations: each column is a word (i.e., letter sequence)
%--------------------------------------------------------------------------
Nc    = 8;                              % number of characters per word
Nw    = numel(word);                    % number of words per phrase
o     = repmat(double(' ') - 31,Nc,Nw);
for i = 1:Nw
    o(1:numel(word{i}),i) = double(word{i}) - 31;
end

% return a cell for this word
%--------------------------------------------------------------------------
o = {repmat(o,1,2)};

if nargout < 2, return, end

% convert to probability cell array
%--------------------------------------------------------------------------
O     = cell(Nc,Nw);
for i = 1:Nc
    for j = 1:Nw
        if o(i,j) > 1
            O{i,j} = sparse(o(i,j),1,1,96,1);
        else
            O{i,j} = spm_dir_norm(ones(96,1));
        end
    end
end

% return a cell for this word
%--------------------------------------------------------------------------
O = {repmat(O,1,2)};

return


function phrase = spm_ind2word(o)
% FORMAT phrase = spm_ind2word(o)
% word  -   word string
% Nc    -   number of characters (with trailing spces)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

%--------------------------------------------------------------------------
phrase = [];
for  i = 1:size(o,2)
    phrase = [phrase,deblank(char(o(:,i)' + 31)),' '];
end

return

