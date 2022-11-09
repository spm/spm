function [u,P] = spm_fx_tfm_P(u,P)
% Exogenous input and input dependent parameters
% FORMAT [u,P] = spm_fx_tfm_P(u,P)
%
% arguments:
% u  - inputs
% P  - parameters
%
% returns:
% u  - exogenous (conductance) inputs driving states
% P  - input dependent parameters
%
% This is a help routine for the microcircuit models equations of motion -
% it simply separates inputs into those affecting (driving) his neuronal
% states and those modulating parameters. It returns the exogenous
% (conductance) inputs and input dependent parameters.
%___________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging
 

% input dependent (intrinsic connection) parameters
%==========================================================================
j     = [4 3 2 1];
for i = 2:size(P.C,2)
    P.G(:,j(i - 1)) = P.G(:,j(i - 1)) + P.C(:,i) + u(i);
end

% exogenous inputs
%--------------------------------------------------------------------------
u     = exp(P.C(:,1))*u(1);
