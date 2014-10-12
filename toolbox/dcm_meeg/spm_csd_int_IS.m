function [y] = spm_csd_int_IS(P,M,U)
% wrapper for erp and csd response of a neural mass model
% FORMAT [y] = spm_csd_int_IS(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - time-dependent input
%
% y{1}  - erp
% y{2}  - csd
%__________________________________________________________________________
%
% This integration routine evaluates the responses of a neural mass model
% to exogenous input - in terms of neuronal states. These are then used as
% expansion point to generate complex cross spectral responses due to
% random neuronal fluctuations. The ensuing spectral (induced) response is
% then convolved (in time) with a window that corresponds to the window of
% a standard wavelet transform. In other words, this routine generates
% predictions of data features based upon a wavelet transform
% characterisation of induced responses.
%__________________________________________________________________________
% Copyright (C) 2012-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_int_IS.m 6234 2014-10-12 09:59:10Z karl $


% check input - default: one trial (no between-trial effects)
%--------------------------------------------------------------------------
[erp,csd] = spm_csd_int(P,M,U);

% concatenate erp and csd into cell array
%--------------------------------------------------------------------------
y = {};
for i = 1:numel(erp)
    y{end + 1} = erp{i};
    y{end + 1} = csd{i};
end