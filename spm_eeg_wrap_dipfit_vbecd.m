function [y,outside,leads]=spm_eeg_wrap_dipfit_vbecd(P,M,U)
% A cost function/wrapper to sit between non-linear optimisation spm_nlsi_gn.m
% and dipole fit routine spm__eeg_inv_vbecd.m
% sens and vol structures should be passed in M, where
%   sens=M.Setup.forward.sens;
%   vol=M.Setup.forward.vol;
% P contains a list of the free parameters (assuming all position
%   parameters come first (in triplets) followed by all moment paameters
%   (also in triplets)
% U is unused
% At the moment this removes the mean level from EEG data
% and reduces the rank of the MEG leadfield 2 dimensions.
%% leads are the lead fields of the dipoles fit
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging
%
% $Id: spm_eeg_wrap_dipfit_vbecd.m 5777 2013-12-04 16:18:12Z vladimir $

x=U.u; %% input , unused


sens=M.Setup.forward.sens;
vol=M.Setup.forward.vol;

siunits = M.Setup.forward.siunits;

chanunits = M.Setup.forward.chanunits;

posandmom=P;

Pospars=3;
Mompars=3;
Ndippars=Pospars+Mompars;
Ndips=length(posandmom)/Ndippars;
if Ndips~=round(Ndips),
    error('only works for 6 params per dipole');
end; % if
allpos=reshape(P(1:Ndips*Pospars),Pospars,Ndips)';
allmom=reshape(P(Ndips*Pospars+1:Ndips*Ndippars),Mompars,Ndips)';

if ft_senstype(sens, 'meg')
    RANK=2; %% restricting rank of MEG data, could change this in future
else
    RANK = 3;
end

y=0;
outside=0;
leads=zeros(Ndips,3,numel(sens.label));
for i=1:Ndips,
    
    pos=allpos(i,:);
    %%mom=allmom(i,:)./1000; %% SCALE BACK FROM SIMILAR UNITS TO LOCATION;
    mom=allmom(i,:); %% in nAm
    % mean correction of LF, only for EEG data.
    
    if siunits
        [tmp] = ft_compute_leadfield(1e-3*pos, sens, vol, 'reducerank',RANK,  'dipoleunit', 'nA*m', 'chanunit', chanunits);
        outside = outside+ ~ft_inside_vol(1e-3*pos,vol);
    else
        [tmp] = ft_compute_leadfield(pos, sens, vol, 'reducerank',RANK);
        outside = outside+ ~ft_inside_vol(pos,vol);
    end
    
    gmn=tmp;
    leads(i,:,:)=gmn';
    if siunits
        rescale = 1;
    else
        rescale=1e3*1e9; %%%% NEED TO DO THIS PROPERLY -GRB MAY 2013
    end
    y=y+gmn*mom'.*rescale;
    
end; % for i


y=y*M.sc_y; %% scale data appropriately
if outside
    y=y.^2;
end;  % penalise sources outside head




