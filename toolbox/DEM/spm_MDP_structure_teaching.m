function [O,o] = spm_MDP_structure_teaching(MDP,OPTIONS)
% Generates probabilistic training sequences for structure learning
% FORMAT [O,o] = spm_MDP_structure_teaching(MDP,[OPTIONS])
%
% MDP     - generative process or exemplars
%
% O{{.}}  - probabilistic exemplars or training sequence
% o{[.]}  - indices of exemplars or training sequence
%
% OPTIONS.N  [0]   - suppress neuronal responses
% OPTIONS.P  [0]   - suppress plotting
% OPTIONS.B  [0]   - replay
% OPTIONS.G  [0]   - suppress graphics
%
% This routine generates a sequence of (probabilistic) outcomes from a
% specified POMDP structure. This sequence is appropriate for structure
% learning. It comprises a sequence of epochs, where each epoch is
% generated in a specific order: starting from the first factor, the
% outcomes associated with each hidden state are generated under the first
% path. By construction, the first path is stationary. After all hidden
% states have generated the outcomes, successive paths are generated
% starting from each hidden state. After all paths have been generated, the
% process is repeated for subsequent factors; under the first state and
% path of previous factors (noting, that the first path is always
% stationary; i.e., an identity transition mapping). Unless otherwise
% specified, outcomes comprise two observations.
%
% This routine is used in conjunction with spm_MDP_structure_learning.m
% that will, in principle, recover the factorial structure of MDP given,
% this sequence of observations.
%
% See: spm_MDP_structure_learning
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8446 2023-06-11 15:20:29Z karl $
%__________________________________________________________________________

% options for model inversion (and evaluation)
%==========================================================================
try OPTIONS.N; catch, OPTIONS.N = 0; end     % suppress neuronal responses
try OPTIONS.P; catch, OPTIONS.P = 0; end     % suppress plotting
try OPTIONS.B; catch, OPTIONS.B = 0; end     % replay
try OPTIONS.G; catch, OPTIONS.G = 0; end     % suppress graphics

% inversion scheme
%--------------------------------------------------------------------------
spm_evaluate = @spm_MDP_VB_XXX;

% generate outcomes
%==========================================================================

% number of outcomes, states, controls and policies
%--------------------------------------------------------------------------
[Nf,Ns,Nu] = spm_MDP_size(MDP);

% ensure factors with dynamics generate sequences first
%--------------------------------------------------------------------------
if any(diff(Nu) > 0)
    warning('Please reorder B{i} with descending number of paths: size(B{i},3)')
end

% number of outcomes (T)
%--------------------------------------------------------------------------
MDP.T = 2;                              % for state transitions
MDP.U = zeros(1,Nf);                    % remove any control
MDP.o = [];                             % remove any outcomes
MDP.O = {};                             % remove any outcomes

% cycle over latent factors to generate outcomes
%--------------------------------------------------------------------------
o     = {};                             % initialise outcome cell array
O     = {};                             % initialise outcome cell array
for i = 1:Nf                            % cycle over factors

    % cycle over paths
    %----------------------------------------------------------------------
    for j = 1:Nu(i)

        % cycle over states
        %------------------------------------------------------------------
        for k = 1:Ns(i)

            % generate outcomes from initial states (and paths) of this
            % factor, under the first states (and paths) of others
            %--------------------------------------------------------------
            MDP.s = ones(Nf,1); MDP.s(i) = k;
            MDP.u = ones(Nf,1); MDP.u(i) = j;
            PDP   = spm_evaluate(MDP,OPTIONS);

            % outcomes for this epoch
            %==============================================================
            o{end + 1} = PDP.o;
            O{end + 1} = PDP.O;
        end
    end
end

return

