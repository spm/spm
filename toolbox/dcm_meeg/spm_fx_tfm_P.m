function [u,P] = spm_fx_tfm_P(u,P)
% returns exogenous input and input dependent parameters
% FORMAT [u,P] = spm_fx_tfm_P(u,P)
%
% arguments:
% u  – inputs
% P  – parameters
%
% returns:
% u  – exogenous (conductance) inputs driving states
% P  – input dependent parameters
%
% tthis is a help are routine for the, microcircuit models equations of
% motion – it simply separates inputs into those affecting (driving) his
% neuronal states and those modulating parameters. It returns the exogenous
% (conductance) inputs and input dependent parameters.
%___________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_tfm_P.m 4714 2012-04-10 13:30:44Z karl $
 

% input dependent (intrinsic connection) parameters
%==========================================================================
j     = [4 3 2 1];
for i = 2:size(P.C,2)
    P.G(:,j(i - 1)) = P.G(:,j(i - 1)) + P.C(:,i) + u(i);
end

% exogenous inputs
%--------------------------------------------------------------------------
u     = exp(P.C(:,1))*u(1);


