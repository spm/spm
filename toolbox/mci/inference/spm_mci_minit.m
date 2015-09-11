function [M] = spm_mci_minit (M)
% Check and initialise model strucuture
% FORMAT [M] = spm_mci_minit (M)
%
% eg. Pre-compute quantities for computing log-joint
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_minit.m 6548 2015-09-11 12:39:47Z will $

if isstruct(M.pE)
    M.vpE=spm_vec(M.pE);
else
    M.vpE=M.pE;
end
if isfield(M,'Ce')
    M.logdet_Ce=spm_logdet(M.Ce);
    try
        M.iCe=M.iCe;
    catch
        M.iCe = inv(M.Ce);
    end
end

M = spm_mci_priors (M);