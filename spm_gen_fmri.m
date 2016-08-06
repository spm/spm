function [y] = spm_gen_fmri(P,M,U)
% Generates a prediction of trial-specific source activity
% FORMAT [y] = spm_gen_erp(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-effects
%   U.u  - inputs
%   U.dt - time bins for within-trial effects
%
% y - BOLD predictions
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_fmri.m 6855 2016-08-06 10:06:35Z karl $

% integrate neural system - when imput changes 
%--------------------------------------------------------------------------

N.f = M.fn;
N.x = M.x.x;
N.m = M.m;
u   = spm_int_L(P,N,U);

% neurovascular imput in terms of intrinsic and extrinsic afferetns
%--------------------------------------------------------------------------


% integrate dynamics and generate prediction
%--------------------------------------------------------------------------
N.f = @spm_fx_hdm;
N.g = @spm_gx_hdm;
M.x = M.x.h;
N.m = M.m;
y   = spm_int_L(P,N,u);