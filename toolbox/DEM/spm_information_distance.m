function [D,C] = spm_information_distance(a)
% information geometry of a likelihood mapping
% FORMAT [D,C] = spm_information_distance(a)
% a{g}    - Dirichlet tensor for modality g
%
% D       - information distance (nats);
% C       - correlation (similarity) matrix (normalised)
%
% This routine uses the information geometry inherent in a likelihood
% mapping. It computes a similarity matrix (C) based upon the information
% distance (D) between columns of a likelihood tensor (i.e., implicit
% categorical distributions), or vice versa.
%__________________________________________________________________________
% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging

% correlation distance : likelihood mapping
% -------------------------------------------------------------------------
Ng  = size(a,1);
C   = spm_cat(a);
C   = spm_cov2corr(C'*C);                     % correlation matrix
D   = 2*sqrt(2)*Ng*(1 - C);                   % distance matrix 

return

% complementary derivation
%==========================================================================

% information distance : likelihood mapping
% -------------------------------------------------------------------------
Ng  = size(a,1);
Ns  = size(spm_cat(a(1,:)),2);
D   = zeros(Ns,Ns);
for g = 1:Ng
    q  = sqrt(spm_cat(a(g,:)));
    for i = 1:Ns
        for j = 1:Ns
            D(i,j) = D(i,j) + 2*sqrt(sum((q(:,i) - q(:,j)).^2));
        end
    end
end

% similarity matrix (assuming normalised vectors on a hypersphere)
% -------------------------------------------------------------------------
C   = 1 - D/(2*sqrt(2)*Ng);                   % correlation matrix

return

