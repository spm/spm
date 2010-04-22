function [gmn, gm, dgm] = spm_eeg_inv_vbecd_getLF(s, sens, vol, step) 
% Estimation of the leadfield matrix and its spatial derivative if required 
% for a set of dipoles used in the VB-ECD solution
%% scales non-eeg data up by a fixed factor (1e8) for compatibility of
%% units
%
% FORMAT [gmn, gm, dgm] = spm_eeg_inv_vbecd_getLF(s, sens, vol, step)
% 
% s      - location vector
% sens   - sensor locations (MNI [mm])
% vol    - volume structure needed by fieldtrip
% step   - stepsize to compute numerical derivatives
%
% gmn    - leadfields (three vectors for each dipole)
% gm     - vectorized leadfields
% dgm    - vectorized partials wrt locations
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips & Stefan Kiebel
% $Id: spm_eeg_inv_vbecd_getLF.m 3833 2010-04-22 14:49:48Z vladimir $


%%% now does rank reduction (to 2) for non eeg data

MEGRANK=2;
MEGSCALE=1e10;  %% should be the same as in spm_eeg_inv_vbecd.m

if nargin<4
     step = 0;
end

gm = [];
for i = 1:length(s)/3
  
    
    % mean correction of LF, only for EEG data.
    if ft_senstype(sens, 'eeg')
       [tmp] = ft_compute_leadfield(s(1+(i-1)*3:i*3)', sens, vol);
        tmp = tmp - repmat(mean(tmp), size(tmp,1), 1);
    else %% reduce rank of leadfield for MEG- assume one direction (radial) is silent
       [tmp] = ft_compute_leadfield(s(1+(i-1)*3:i*3)', sens, vol,'reducerank',MEGRANK);
        tmp=tmp.*MEGSCALE;
    end
    gm = [gm tmp];
end

gmn = gm; % leadfield

[Nc, Np] = size(gmn);
if nargout >= 2
    gm = gmn(:); % vectorized leadfield
end


if all(step > 0) && nargout == 3,
%     if isempty(step),
%         step=randn(size(s));
%     end; % if isempty
    dgm = [];
    for j = 1:length(s)
        ds = s;
        ds(j) = s(j) + step(j);
        dtmp = [];
        for i = 1:length(s)/3
            if ceil(j/3) == i 

                % mean correction of LF, only for EEG data.
                if ft_senstype(sens, 'eeg')
                    [tmp] = ft_compute_leadfield(ds(1+(i-1)*3:i*3)', sens, vol);
                    tmp = tmp - repmat(mean(tmp), size(tmp,1), 1);
                else  % MEG
                   [tmp] = ft_compute_leadfield(ds(1+(i-1)*3:i*3)', sens, vol,'reducerank',MEGRANK);
                    tmp=tmp.*MEGSCALE;
                end
                dtmp = [dtmp tmp];
            else
                dtmp = [dtmp gmn(:, 1+(i-1)*3:i*3)];
            end
        end
        dtmp = dtmp(:);
        dgm = [dgm (-gm + dtmp)./step(j)];
    end
    
    % correct order
    ind = reshape(1:Np^2 , Np, Np)';
    
    dgm = reshape(dgm, size(gmn, 1), Np^2);
    dgm = dgm(:, ind(:));
    dgm = reshape(dgm, Np*Nc, Np);    
end

