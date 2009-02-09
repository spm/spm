function [gmn, gm, dgm] = spm_eeg_inv_vbecd_getLF(s, sens, vol, step) 
% FORMAT [gmn, gm, dgm] = spm_eeg_inv_vbecd_getLF(s, sens, vol, channels, step, Bad)
%
% Estimation of the leadfield matrix and is spatial derivative if required 
% for a set of dipoles used in the VB-ECD solution
%
% inputs:
%   s    - location vector
%   sens - sensor locations (MNI [mm])
%   vol  - volume structure needed by fieldtrip
%   step - stepsize to compute numerical derivatives
%
% outputs:
%   gmn  - leadfields (three vectors for each dipole)
%   gm   - vectorized leadfields
%   dgm  - vectorized partials wrt locations (if 3rd output argument specified)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips & Stefan Kiebel
% $Id: spm_eeg_inv_vbecd_getLF.m 2720 2009-02-09 19:50:46Z vladimir $
 
gm = [];
for i = 1:length(s)/3
    [tmp] = forwinv_compute_leadfield(s(1+(i-1)*3:i*3)', sens, vol);
    % mean correction of LF, only for EEG data.
    if forwinv_senstype(sens, 'eeg')
        tmp = tmp - repmat(mean(tmp), size(tmp,1), 1);
    end
%     tmp(Bad, :) = [];
    gm = [gm tmp];
end

gmn = gm; % leadfield

[Nc, Np] = size(gmn);
if nargout >= 2
    gm = gmn(:); % vectorized leadfield
end

if step > 0
    dgm = [];
    for j = 1:length(s)
        ds = s;
        ds(j) = s(j) + step(j);
        dtmp = [];
        for i = 1:length(s)/3
            if ceil(j/3) == i 
                [tmp] = forwinv_compute_leadfield(ds(1+(i-1)*3:i*3)', sens, vol);
                tmp = tmp - repmat(mean(tmp), size(tmp,1), 1);
%                 tmp(Bad, :) = [];
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
