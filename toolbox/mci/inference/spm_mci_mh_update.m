function [next,accepted,bayes_fb,dL] = spm_mci_mh_update (curr,prop,verbose)
% Update parameters using Metropolis-Hastings
% FORMAT [next,accepted,bayes_fb,dL] = spm_mci_mh_update (curr,prop,verbose)
%
% curr      quantities re current state
% prop      quantities re proposed state
% verbose   1 for text output
%
% next      next state
% accepted  1 for accepted proposal
% bayes_fb  Log Bayes factor for forward versus backward transition
% dL        Proposed difference in log joint
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_mh_update.m 6548 2015-09-11 12:39:47Z will $

accepted=1;

if isempty(curr)
    % must be first move - accept proposal
    next = prop;
    bayes_fb = 0;
    dL = 0;
else
    % Evaluate transition probabilities forward and backward
    ef = prop.pos-curr.mu;
    Lforward = -0.5*curr.logdetCp - 0.5*ef'*curr.iCp*ef;
    eb = curr.pos-prop.mu;
    Lback = -0.5*prop.logdetCp - 0.5*eb'*prop.iCp*eb;
    bayes_fb = Lforward - Lback;
    dL = prop.L - curr.L;
    
    Ltot = prop.L - curr.L + Lback - Lforward;
    r = exp(Ltot);
    alpha = min(1,r);
    test_prob = rand(1);
    
    % Hastings ratio
    if alpha > test_prob
        if verbose
            display(['*********** sample accepted *****************']);
        end
        % accept and update parameters
        next = prop;
    else
        % reject move
        next = curr;
        accepted = 0;
    end
end