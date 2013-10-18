function [g] = spm_dcm_fmri_graph_gen(x,v,P)
% Generates adjacency graph for DCM for CSD (fMRI)
% FORMAT [g] = spm_dcm_fmri_graph_gen(x,v,P)
%
% This routine computes the connecivity graph for DCM
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_graph_gen.m 5702 2013-10-18 11:10:06Z karl $


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
       d  = d/2;
       r  = r*2;

       % Forward and backward bias
       %-------------------------------------------------------------------
       A(i,j,1) = -log(1 + exp(-r)) - d;
       A(i,j,2) = -log(1 + exp( r)) - d;
       A(j,i,1) = A(i,j,2);
       A(j,i,2) = A(i,j,1);
       
    end
end

% add bias to (empirical) prior mean of log connectivity
%--------------------------------------------------------------------------
g   = DCM.M.pE;
g.A = g.A + A;

return

% Notes to illustrate the connection strengths dependencies on distance
%--------------------------------------------------------------------------
r      = -4:1/32:4;
d      = abs(r);
d      = d/2;
r      = r*2;

A(:,1) = -log(1 + exp(-r*k)) - d;
A(:,2) = -log(1 + exp( r*k)) - d;
AA     = exp(A(:,1)) - exp(A(:,2));

subplot(2,2,1)
plot(r,A)
xlabel('hierarchical distance')
xlabel('log connection strength')
title('Log connectivity','FontSize',16)
legend({' forward',' backward'})

subplot(2,2,2)
plot(r,AA)
xlabel('hierarchical distance')
xlabel('Log connection strength')
title('Connectivity','FontSize',16)







