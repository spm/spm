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
% $Id: spm_dcm_symm.m 5376 2013-04-02 09:59:01Z karl $

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
%==========================================================================
n         = size(pE.Lpos,2);
Rpos      = pE.Lpos;
Rpos(1,:) = -Rpos(1,:);
D         = 128*ones(n);

% find symmetrical sources in each hemisphere
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        if sign(pE.Lpos(1,i)) == sign(Rpos(1,j))
            D(i,j) = sqrt(sum(pE.Lpos(:,i) - Rpos(:,j)).^2);
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

% reduce rank of prior covariance matrix of positions
%--------------------------------------------------------------------------
try
    pV.L = pV.L + kron(DD,diag(pV.L(1)*[-1 1 1]));
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


