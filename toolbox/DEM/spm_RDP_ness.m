function [MDP] = spm_RDP_ness(MDP,S,chi)
% BMR using eigen-decomposition of allowable transitions 
% FORMAT [MDP] = spm_RDP_ness(MDP,S,chi)
% MDP - hierarchal RDP
% S    - list of modality streams;    e.g., [2, 3]
% chi  - log prior for each modality; e.g., [8,-8] 
%
% This routine uses the eigen-decomposition of allowable transitions (i.e.,
% flow among states) to identify steady-state distributions (i.e., the
% principal eigenvector with a real eigenvalue of unity). These
% nonequilibrium steady-state (NESS) distributions correspond to the
% attracting set; i.e., a pullback attractor.
%
% The derivative of expected cost is then evaluated with respect to
% reducing the transition probabilities to each state in turn (provided
% it's removal does not render another state childless). If reducing the
% probability of visiting a generalised state increases expected goals
% (i.e., expected occupancy of goal states under the nonequilibrium
% steady-state density), it is removed.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


% prune absorbing states
%--------------------------------------------------------------------------
MDP   = spm_RDP_prune(MDP);

% get constraints
%--------------------------------------------------------------------------
MDP   = spm_set_goals(MDP,S,chi);
h     = MDP{end}.id.hid;

% eigenvalue decomposition of flows
%--------------------------------------------------------------------------
B     = spm_dir_norm(sum(MDP{end}.b{1},3) > 0);
Ns    = size(B,1);
[e,v] = eig(B,'nobalance');
[~,j] = max(real(diag(v)));
e     = spm_dir_norm(abs(e(:,j)));

% Find states that can be removed without creating absorbing states
%--------------------------------------------------------------------------
k     = any(B == 1,2) | ~any(B,2);
k     = find(~k)';

% And pre-compute matrices for evaluating partial derivatives
%--------------------------------------------------------------------------
iB    = pinv(eye(Ns,Ns) - B);
iB    = max(chi)*sum(iB(h,:),1);
dr    = 1/128;
dGdr  = zeros(1,Ns);
for i = k
    b       = B;
    j       = logical(b(i,:));
    b(i,j)  = b(i,j)*exp(dr);
    b(:,j)  = spm_dir_norm(b(:,j));
    dBdr    = (b(:,j) - B(:,j))/dr;
    dGdr(i) = iB*(dBdr*e(j));
end

% Retain states that increase expected goal occupancy (E) at NESS
%--------------------------------------------------------------------------
r  = dGdr > -1/128;

% restriction matrix (R) and reduction of MDP
%--------------------------------------------------------------------------
R   = speye(Ns,Ns);
MDP = spm_RDP_compress(MDP,R(:,r));

% prune absorbing states
%--------------------------------------------------------------------------
MDP = spm_RDP_prune(MDP);

% and set constraints
%--------------------------------------------------------------------------
MDP = spm_set_goals(MDP,S,chi);

return


% NOTES
%--------------------------------------------------------------------------
% [e,v] = eig(B,'nobalance');
% [v,j] = max(real(diag(v)));
% e     = spm_dir_norm(abs(e(:,j)));
% R     = sum(e(h));
% r     = true(1,Ns);
% for i = 1:Ns
%     b      = B;
%     b(i,:) = 0;
%     if all(any(b,1))
%         b      = spm_dir_norm(b);
%         [e,v]  = eig(b,'nobalance');
%         [v,j]  = max(real(diag(v)));
%         e      = spm_dir_norm(abs(e(:,j)));
%         E      = sum(e(h));
%         if E > (R + exp(-16))
%             R = E;
%             B = b;
%             r(i) = false;
%         end
%     end
% end
