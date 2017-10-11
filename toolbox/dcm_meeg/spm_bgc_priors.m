function [pE,pC] = spm_bgc_priors
% Prior moments for a basal ganglia circuit
% FORMAT [pE,pC] = spm_bgc_priors
% only contains priors for intrinsic parameters
% priors for extrinsic parameters are defined in spm_cmc_priors
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging
 
% Bernadette van Wijk
% $Id: spm_bgc_priors.m 7185 2017-10-11 10:10:04Z spm $


% synaptic parameters
%--------------------------------------------------------------------------

E.T  = sparse(1,5);   V.T  = [1/4 1/4 1/4 1/4 1/4];

E.G=[ 0          0       0       0       0       0       0       0       0];
V.G=[ 1/2    1/2    1/2    1/2    1/2   1/2    1/2    1/2    1/2]; 

E.S  = 0;             V.S  = 1/16;                % slope of sigmoid  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pE     = E;
pC     = V;