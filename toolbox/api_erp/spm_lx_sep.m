function [L] = spm_lx_sep(P,M)
% observer matrix for a neural mass model of erps: y = G*x
% FORMAT [G] = spm_lx_sep(P,M)
% x      - state vector
%   x(:,1) - voltage (spiny stellate cells)
%   x(:,2) - voltage (pyramidal cells) +ve
%   x(:,3) - voltage (pyramidal cells) -ve
%   x(:,4) - current (spiny stellate cells)    depolarizing
%   x(:,5) - current (pyramidal cells)         depolarizing
%   x(:,6) - current (pyramidal cells)         hyperpolarizing
%   x(:,7) - voltage (inhibitory interneurons)
%   x(:,8) - current (inhibitory interneurons) depolarizing
%   x(:,9) - voltage (pyramidal cells)
%
% G        - y = G*x
%
% This is like spm_erp_lx but allows each subpopulation within a source to
% have its own lead field
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lx_sep.m 1174 2008-02-27 20:22:30Z karl $
 
% parameterised lead field ECD for each population
%--------------------------------------------------------------------------
J      = [1 7 9];                            % contributing states
nc     = M.l;                                % number of channels
ns     = M.n/9;                              % number of sources
L      = cell(1,9);
[L{:}] = deal(sparse(nc,ns));
for  i = 1:length(P.L)
    Q       = P;
    Q.L     = P.L{i};
    L{J(i)} = spm_erp_L(Q,M);
end
 
% concatenate
%--------------------------------------------------------------------------
L = spm_cat(L);
