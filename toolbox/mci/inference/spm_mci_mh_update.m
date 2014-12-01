function [x,E,accepted] = mh_update (x,E,pos,verbose)
% Update parameters using Metropolis-Hastings
% FORMAT [x,E,accepted] = mh_update (x,E,pos,verbose)
%
% x         parameters
% E         energy
% pos       proposal
% verbose   1 for text output
%
% x         parameters
% E         energy
% accepted  1 for accepted proposal
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_mh_update.m 6275 2014-12-01 08:41:18Z will $

i=length(E);
new_energy=E(i);
accepted=1;

if i > 1
    old_energy = E(i-1);
    compratio = exp(old_energy - new_energy);
    alpha = min(1,compratio);
    test_prob = rand(1);
    
    % Hastings ratio
    if alpha > test_prob
        if verbose
            display(['*********** sample accepted *****************']);
        end
        % accept and update parameters
        x(i,:) = pos;
    else
        % reject move
        x(i,:) = x(i-1,:);
        E(i) = E(i-1);
        accepted = 0;
    end
else
    % Always accept first move
    x(i,:) = pos;
end