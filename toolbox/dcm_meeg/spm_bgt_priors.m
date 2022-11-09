function [pE,pC] = spm_bgt_priors
% Prior moments for a basal ganglia circuit
% FORMAT [pE,pC] = spm_bgt_priors
% only contains priors for intrinsic parameters
% priors for extrinsic parameters are defined in spm_cmc_priors
%__________________________________________________________________________
 
% Bernadette van Wijk
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging

% synaptic parameters
%--------------------------------------------------------------------------

E.T  = sparse(1,5);   V.T  = [1/4 1/4 1/4 1/4 1/4];

E.G=[ 0          0       0       0       0       0       0       0       0];
V.G=[ 1/2    1/2    1/2    1/2    1/2   1/2    1/2    1/2    1/2]; 

E.S  = 0;             V.S  = 1/16;                % slope of sigmoid  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pE     = E;
pC     = V;