function [x] = spm_dcm_neural_x(P,M)
% Returns the the fixed point or steady-state of a neural mass DCM
% FORMAT [x,f] = spm_dcm_x_neural(P,'model')
%
%  P      - parameter structure
% 'model'   - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%
% x   - initial states
% f   - state euquation
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_neural_x.m 4718 2012-04-19 15:34:45Z karl $
 
 
% solve for fixed point (with 64 ms burn in) - if no exogenous input
%--------------------------------------------------------------------------
if strcmp(M.f,'spm_fx_mfm')
    M.g  = {};
    U.u  = sparse(8,M.m);
    U.dt = 8/1000;
    x    = spm_int_L(P,M,U);
    x    = spm_unvec(x(end,:),M.x);
else
    x    = M.x;
end