function [g] = spm_dcm_fmri_graph_gen(x,v,P)
% Generates adjacency graph for DCM for CSD (fMRI)
% FORMAT [g] = spm_dcm_fmri_graph_gen(x,v,P)
%
% This routine computes the connecivity graph for DCM
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_graph_gen.m 5696 2013-10-15 19:10:26Z karl $


% global DCM and evaluate orginal generative model
%==========================================================================
global DCM

% compute bias for log connectivity using functional space
%==========================================================================

% Distance matrix D (and radial distance R)
%--------------------------------------------------------------------------
m     = size(v,2);
for i = 1:m
    for j = (i + 1):m
       
       % Distance - d (and radial -r)
       %-------------------------------------------------------------------
       d  = sqrt(sum((v(:,i) - v(:,j)).^2));
       r  = sqrt(sum(v(:,j).^2)) - sqrt(sum(v(:,i).^2));

       % Forward and backward bias
       %-------------------------------------------------------------------
       B(i,j,1) = -2*log(1 + exp(-r)) - d/2;
       B(i,j,2) = -2*log(1 + exp( r)) - d/2;
       B(j,i,1) = B(i,j,2);
       B(j,i,2) = B(i,j,1);
       
    end
end

% add bias to (empirical) prior mean of log connectivity
%--------------------------------------------------------------------------
g   = DCM.M.pE;
g.A = g.A + B;




