function [L] = spm_lx_erp(P,dipfit)
% observer matrix for a neural mass model: y = G*x
% FORMAT [G] = spm_lx_erp(P,dipfit)
% FORMAT [G] = spm_lx_erp(P,M)
%
% M.dipfit - spatial model specification
%
% G        - where y = L*x; G = dy/dx
% x        - state vector
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lx_erp.m 6720 2016-02-15 21:06:55Z karl $

% extract dipfit from model if necessary
%--------------------------------------------------------------------------
if  isfield(dipfit,'dipfit'), dipfit = dipfit.dipfit; end
if ~isfield(dipfit,'type'),   dipfit = 'LFP';         end

% parameterised lead field times source contribution to ECD
%--------------------------------------------------------------------------
L     = spm_erp_L(P,dipfit);               % lead field per source
if isnumeric(P.J)
    L = kron(P.J,L);                       % lead-field per state
else
    
    % construct lead field for each source
    %----------------------------------------------------------------------
    for i = 1:numel(P.J)
        G{i} = L(:,i)*P.J{i};
    end
    L = spm_cat(G);
end
