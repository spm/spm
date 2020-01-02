function [y,outside,leads] = spm_eeg_wrap_momfit_vbecd(P,M,U)
% A cost function/wrapper to sit between non-linear optimisation spm_nlsi_gn.m
% and dipole fit routine spm_cfg_eeg_momentfit.m
% FORMAT [y,outside,leads] = spm_eeg_wrap_momfit_vbecd(P,M,U)
% sens and vol structures should be passed in M, where
%   sens = M.Setup.forward.sens;
%   vol  = M.Setup.forward.vol;
% P contains a list of the free parameters (assuming all position
%   parameters come first (in triplets) followed by all moment paameters
%   (also in triplets)
% U is unused
% At the momnent reduces the rank of the MEG leadfield 2 dimensions.
% leads are the lead fields of the dipoles fit
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_eeg_wrap_momfit_vbecd.m 7764 2020-01-02 15:34:03Z spm $


x = U.u; % input , unused

sens = M.Setup.forward.sens;
vol  = M.Setup.forward.vol;

siunits = M.Setup.forward.siunits;
if ~siunits
    warning('Data not in SI units, scaling (and therefore priors) maybe wrong');
end
chanunits = M.Setup.forward.chanunits;

mompars = P; % one moment parameter per source
Ndips   = length(mompars);
allpos  = M.pos; % source positions
allori  = M.ori; % source orientations

if (size(allpos,1)~=size(allori,1))|| (size(allpos,1)~=length(mompars))
    error('There must be equal number of moments, positions and orientations');
end

Ndips = length(mompars);

allmom = allori.*repmat(spm_vec(mompars),1,3);

if ft_senstype(sens, 'meg')
    RANK = 2; % restricting rank of MEG data, could change this in future
else
    RANK = 3;
end

y = 0;
outside = 0;
leads = zeros(Ndips,3,numel(sens.label));
for i=1:Ndips
    
    pos = allpos(i,:);
    
    mom = allmom(i,:); %% in nAm
    
    if siunits
        [tmp] = ft_compute_leadfield(1e-3*pos, sens, vol, 'reducerank',RANK, 'dipoleunit', 'nA*m', 'chanunit', chanunits);
        
    else
        [tmp] = ft_compute_leadfield(pos, sens, vol, 'reducerank',RANK);
        
    end
    
    gmn = tmp;
    leads(i,:,:) = gmn';
    
    y = y + gmn * mom';
    
end % for i
