function Q = spm_dcm_csd_Q(csd)
% Precision of cross spectral density
% FORMAT Q = spm_dcm_csd_Q(csd)
% 
% csd{i}   - [cell] Array of complex cross spectra
% Q        - normalised precision
%--------------------------------------------------------------------------
% This routine returns the precision of complex cross spectra based upon
% the asymptotic results described in Camba-Mendez & Kapetanios (2005):
% In particular, the scaled difference between the sample spectral density
% (g) and the predicted density (G);
%  
% e = vec(g - G)
% 
% is asymptotically complex normal, where the covariance between e(i,j) and
% e(u,v) is given by Q/h and:
% 
% Q = G(i,u)*G(j,u):  h = 2*m + 1
%  
% Here m represent the number of averages from a very long time series.The
% inverse of the covariance is thus a scaled precision, where the
% hyperparameter (h) plays the role of the degrees of freedom (e.g., the
% number of averages comprising the estimate). In this routine, we use the
% sample spectral density to create a frequency specific precision matrix
% for the vectorised spectral densities - under the assumption that the
% former of this sample spectral density resembles the predicted spectral
% density (which will become increasingly plausible with convergence).
%
% Camba-Mendez, G., & Kapetanios, G. (2005). Estimating the Rank of the
% Spectral Density Matrix. Journal of Time Series Analysis, 26(1), 37-48.
% doi: 10.1111/j.1467-9892.2005.00389.x
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_csd_Q.m 7751 2019-12-06 11:59:09Z peter $


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
Q       = spalloc(Qn,Qn,SIZ(1)*prod(SIZ(2:end))^2);
[w,i,j] = ind2sub(SIZ,1:Qn);
for Qi  = 1:Qn
    for Qj = 1:Qn
        if w(Qi) == w(Qj)
            Q(Qi,Qj) = csd(w(Qi),i(Qi),i(Qj))*csd(w(Qi),j(Qi),j(Qj));
        end
    end
end
Q       = inv(Q + norm(Q,1)*speye(size(Q))/32);

