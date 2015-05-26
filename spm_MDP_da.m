function spm_MDP_da(MDP)
% Simulated histograms of dopamine firing
% FORMAT spm_MDP_da(MDP)
%
% See also: spm_MDP_game, which generalises this scheme and replaces prior
% beliefs about KL control with minimisation of expected free energy.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_da.m 6451 2015-05-26 09:26:03Z karl $

% deconvolve to simulate dopamine responses
%--------------------------------------------------------------------------
wn   = MDP.d;
da   = pinv( tril(toeplitz(exp(-((1:length(wn)) - 1)'/8))) )*wn;
da   = da(4:end);

% peristimulus time histogram
%--------------------------------------------------------------------------
pst  = (1:length(da))*256;
r    = 256;
psth = r*da + randn(size(da)).*sqrt(r*da);
bar(pst,psth,1)
title('Simulated dopamine responses','FontSize',16)
xlabel('Peristimulus time (ms)','FontSize',12)
ylabel('Spikes per bin','FontSize',12)
axis([pst(1) pst(end) 0 max(psth)*(1 + 1/4)])
axis square
