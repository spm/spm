function [mh] = spm_mci_popdef (scale,tune,samp)
% Set default parameters for population MCMC
% FORMAT [mh] = spm_mci_popdef (scale,tune,samp)
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

mh.nscale=scale;
mh.ntune=tune;
mh.nsamp=samp;
mh.ind_samp=[scale+tune+1:scale+tune+samp];
mh.J=1; % Number of temperatures
mh.gprob=0;
mh.remove_burn_in=0;
mh.verbose=0;
