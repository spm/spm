function Q  = spm_dcm_csd_Q(csd)
% Precision of cross spectral density
% FORMAT Q  = spm_dcm_csd_Q(csd)
% 
% csd{i}   - [cell] Array of complex cross spectra
% Q        - normalised precision
%--------------------------------------------------------------------------
%  This routine returns the precision of complex cross spectra based upon
%  the asymptotic results described in: 
%  Camba-Mendez, G., & Kapetanios, G. (2005). Estimating the Rank of the
%  Spectral Density Matrix. Journal of Time Series Analysis, 26(1), 37-48.
%  doi: 10.1111/j.1467-9892.2005.00389.x
%
% NB:  Although a very simple routine, it can be very slow for large
% matrices – and may benefit from coding in C.
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_csd_Q.m 7275 2018-03-07 22:36:34Z karl $

% check for cell arrays
%--------------------------------------------------------------------------
if iscell(csd)
    CSD   = spm_zeros(csd{1});
    n     = numel(csd);
    for i = 1:n
        CSD = CSD + csd{i};
    end
    Q  = spm_dcm_csd_Q(CSD/n);
    Q  = kron(eye(n,n),Q);
    return
end

% get precision
%--------------------------------------------------------------------------
SIZ    = size(csd);
Qn     = spm_length(csd);
Q      = sparse(Qn,Qn);
for Qi = 1:Qn
    for Qj = 1:Qn
        [wi,i,j] = ind2sub(SIZ,Qi);
        [wj,u,v] = ind2sub(SIZ,Qj);
        if wi == wj
            Q(Qi,Qj) = csd(wi,i,u)*csd(wi,j,v);
        end
    end
end
Q      = inv(Q + norm(Q,1)*speye(size(Q))/32);