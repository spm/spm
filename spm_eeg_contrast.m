function c = spm_eeg_contrast(SPM, xCon)
% SPM = spm_eeg_contrast(SPM)
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id$

c = 1; 
for j = 1:SPM.eeg.Nfactors            
    c = kron(c, xCon.eeg.Con{1,j});
end

% pad with zero, if multiple design partitions
if SPM.eeg.Ncomp_d > 1
    c = [c; zeros(size(SPM.xX.X, 2) - size(c,1), 1)];
end

% also remove columns!
