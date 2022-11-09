function [y] = spm_csd_int_IS(P,M,U)
% Wrapper for erp and csd response of a neural mass model
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

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


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
