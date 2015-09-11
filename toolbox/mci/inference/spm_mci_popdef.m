function [mh] = spm_mci_popdef (scale,tune,samp)
% Set default parameters for population MCMC
% FORMAT [mh] = spm_mci_popdef (scale,tune,samp)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_popdef.m 6548 2015-09-11 12:39:47Z will $

mh.nscale=scale;
mh.ntune=tune;
mh.nsamp=samp;
mh.ind_samp=[scale+tune+1:scale+tune+samp];
mh.J=1; % Number of temperatures
mh.gprob=0;
mh.remove_burn_in=0;
mh.verbose=0;
