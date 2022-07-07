function [y] = spm_gx_state_fmri(x,u,P,M)
% Simulated BOLD response and copied state vector
% FORMAT [y] = spm_gx_state_fmri(x,u,P,M)
% y          - BOLD response and copied state vector
%
% x          - state vector     (see spm_fx_fmri)
% P          - Parameter vector (see spm_fx_fmri)
% M          - model specification structure (see spm_nlsi)
%
% The `copied state vector' passes the first hidden variable in each region
% to the output variable y, so that 'neural activities' can be plotted 
% by spm_dcm_generate.m
%
% See spm_fx_fmri.m and spm_dcm_generate.m
%__________________________________________________________________________
 
% Will Penny
% Copyright (C) 2011-2022 Wellcome Centre for Human Neuroimaging


y = spm_gx_fmri(x,u,P,M);

% Copy first hidden state (neural activity) from each region
nregions = size(x,1);
for i=1:nregions
    y = [y;x(i,1)];
end
y = full(y);
