function Q = spm_dcm_csd_Q(csd)
% Precision of cross spectral density
% FORMAT Q = spm_dcm_csd_Q(csd)
% 
% csd{i}   - [cell] Array of complex cross spectra
% Q        - normalised precision
%--------------------------------------------------------------------------
% This routine returns the precision of complex cross spectra based upon
% the asymptotic results described in: 
% Camba-Mendez, G., & Kapetanios, G. (2005). Estimating the Rank of the
% Spectral Density Matrix. Journal of Time Series Analysis, 26(1), 37-48.
% doi: 10.1111/j.1467-9892.2005.00389.x
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_csd_Q.m 7481 2018-11-09 15:36:57Z peter $


%-Check for cell arrays
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

%-Get precision
%--------------------------------------------------------------------------
SIZ     = size(csd);
Qn      = spm_length(csd);
[w,i,j] = ind2sub(SIZ,1:Qn);
u       = zeros(SIZ(1)*prod(SIZ(2:end))^2,1);
v       = zeros(size(u));
s       = zeros(size(u));
t       = 1;
for Qi  = 1:Qn
    for Qj = 1:Qn
        if w(Qi) == w(Qj)
            u(t) = Qi;
            v(t) = Qj;
            s(t) = csd(w(Qi),i(Qi),i(Qj))*csd(w(Qi),j(Qi),j(Qj));
            t    = t + 1;
        end
    end
end
Q       = sparse(u,v,s,Qn,Qn);
Q       = inv(Q + norm(Q,1)*speye(size(Q))/32);
