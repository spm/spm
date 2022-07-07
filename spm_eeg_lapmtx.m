function [T] = spm_eeg_lapmtx(pst)
% Laplace transform basis set for ERPs
% FORMAT [T] = spm_eeg_lapmtx(pst)
%
% pst - perstimulus time in ms
%
% T   - Laplace transform basis set
% 
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% assume a single sample if not specified
%--------------------------------------------------------------------------
pst   = pst(:)/pst(end);
pst   = pst.*(pst > 0);
w     = [1:32]*pi;
k     = [1:4];
T     = [];
S     = sparse(2,0);
for i = 1:length(w)
    for j = 1:length(k)
        W            = w(i);
        K            = k(j)*W;
        s            = sqrt(-1)*w(i) + k(j)*i;
        T(:,end + 1) = pst.^3.*exp(-s*pst);
        S(:,end + 1) = [W; K];
    end
end

% select major modes
%--------------------------------------------------------------------------
T    = spm_svd(imag(T));
