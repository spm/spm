function [pC] = spm_dcm_symm(pV,pE)
% locks ECD orientations by introducing prior correlations
% FORMAT [pC] = spm_dcm_symm(pV,pE)
%__________________________________________________________________________
%
% pE   - prior expectation
% pV   - prior variance
% pC   - prior covariance
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_symm.m 4095 2010-10-22 19:37:51Z karl $

% Distance between homolgous sources (16mm)
%--------------------------------------------------------------------------
V     = 16;

% symmetry constraints (based on Euclidean distance from mirror image)
%==========================================================================

% diagonalise feilds
%--------------------------------------------------------------------------
feilds = fieldnames(pV);
for  i = 1:length(feilds)
    pF = getfield(pV,feilds{i});    
    pV = setfield(pV,feilds{i},spm_diag(spm_vec(pF)));
end

% impose correlations between orientations (L)
%--------------------------------------------------------------------------
n     = size(pE.Lpos,2);
Mpos  = [-pE.Lpos(1,:); pE.Lpos(2,:); pE.Lpos(3,:)];
D     = Inf*ones(n);
for i = 1:n
    for j = 1:n
        if sign(pE.Lpos(1,i)) == sign(Mpos(1, j))
            D(i,j) = sqrt(sum(pE.Lpos(:,i) - Mpos(:,j)).^2);
        end
    end
end
D     = (D + D')/2;
DD    = zeros(n);
for i = 1:n
    [M, I] = min(D(i,:));
    if M < V
        DD(i,I) = 1;
    end
end
try
    pV.L = pV.L - kron(DD,diag([1 0 0])) + kron(DD,diag([0 1 1]));
end

% and concatenate
%--------------------------------------------------------------------------
for  i = 1:length(feilds)
    pF      = getfield(pV,feilds{i});
    if ~isempty(pF)
        pC{i,i} = pF;  
    end
end
pC    = spm_cat(pC);


