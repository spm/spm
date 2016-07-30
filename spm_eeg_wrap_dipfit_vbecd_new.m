function [y,leads] = spm_eeg_wrap_dipfit_vbecd_new(P,M,U)
% A cost function/wrapper to sit between non-linear optimisation spm_nlsi_gn.m
% and dipole fit routine spm__eeg_inv_vbecd.m
% sens and vol structures should be passed in M, where
%   sens=M.Setup.forward.sens;
%   vol=M.Setup.forward.vol;
% P contains a list of the free parameters (assuming all position
%   parameters come first (in triplets) followed by all moment paameters
%   (also in triplets)
% At the momnent reduces the rank of the MEG leadfield 2 dimensions.
% leads are the lead fields of the dipoles fit
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging
%__________________________________________________________________________

% $Id: spm_eeg_wrap_dipfit_vbecd_new.m 6847 2016-07-30 10:35:32Z karl $


% forward moral
%==========================================================================
sens      = M.Setup.forward.sens;
vol       = M.Setup.forward.vol;
siunits   = M.Setup.forward.siunits;
chanunits = M.Setup.forward.chanunits;

Nchan     = numel(chanunits);
Nmods     = size(P.q,3);
Ndips     = size(P.p,1);

% restricting rank of MEG data, could change this in future
%--------------------------------------------------------------------------
if ft_senstype(sens, 'meg')
    RANK = 2;
else
    RANK = 3;
end
y     = zeros(Nchan,Nmods);
leads = zeros(Ndips,3,Nchan);
for i = 1:Ndips,
    if siunits
        gmn = ft_compute_leadfield(1e-3*P.p(i,:), sens, vol, 'reducerank',RANK,  'dipoleunit', 'nA*m', 'chanunit', chanunits);
    else 
        gmn = ft_compute_leadfield(P.p(i,:), sens, vol, 'reducerank',RANK);
    end
    leads(i,:,:) = gmn';
    gmn   = gmn*M.Lscal^2;
    for j = 1:Nmods
        y(:,j) = y(:,j) + gmn*P.q(i,:,j)';
    end
    
end; % for i






