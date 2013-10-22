function [P] = spm_dcm_fmri_graph_gen(x,v,P)
% Generates adjacency graph for DCM for CSD (fMRI)
% FORMAT [g] = spm_dcm_fmri_graph_gen(x,v,P)
%
% This routine computes the connecivity graph for DCM
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_graph_gen.m 5708 2013-10-22 09:20:59Z karl $


% compute bias for log connectivity using functional space
%==========================================================================

% Distance-based bias on (empirical) prior mean of log connectivity
%--------------------------------------------------------------------------
m     = size(v.x,2);
for i = 1:m
    for j = (i + 1):m
        
       P.A(i,j) = v.a - sum((v.x(:,i) - v.x(:,j)).^2)/2;
       P.A(j,i) = P.A(i,j);
       
    end
end







