function c = spm_eeg_contrast(SPM, xCon)
% SPM = spm_eeg_contrast(SPM)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_contrast.m 1143 2008-02-07 19:33:33Z spm $

c = 1; 
for j = 1:SPM.eeg.Nfactors            
    c = kron(c, xCon.eeg.Con{1,j});
end

% pad with zero, if multiple design partitions
if SPM.eeg.Ncomp_d > 1
    c = [c; zeros(size(SPM.xX.X, 2) - size(c,1), 1)];
end

% also remove columns!
