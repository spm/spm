function [gmn, gm, dgm] = spm_eeg_inv_vbecd_getLF(s, sens, vol, step) 
% Estimation of the leadfield matrix and its spatial derivative if required 
% for a set of dipoles used in the VB-ECD solution
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
% $Id: spm_eeg_inv_vbecd_getLF.m 3034 2009-04-01 15:12:55Z jean $


if nargin<4
     step = 0;
end

gm = [];
for i = 1:length(s)/3
    [tmp] = forwinv_compute_leadfield(s(1+(i-1)*3:i*3)', sens, vol);
    % mean correction of LF, only for EEG data.
    if forwinv_senstype(sens, 'eeg')
        tmp = tmp - repmat(mean(tmp), size(tmp,1), 1);
    end
    gm = [gm tmp];
end

gmn = gm; % leadfield

[Nc, Np] = size(gmn);
if nargout >= 2
    gm = gmn(:); % vectorized leadfield
end


if all(step > 0) && nargout == 3
    dgm = [];
    for j = 1:length(s)
        ds = s;
        ds(j) = s(j) + step(j);
        dtmp = [];
        for i = 1:length(s)/3
            if ceil(j/3) == i 
                [tmp] = forwinv_compute_leadfield(ds(1+(i-1)*3:i*3)', sens, vol);
                tmp = tmp - repmat(mean(tmp), size(tmp,1), 1);
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

