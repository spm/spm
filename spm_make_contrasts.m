function Con = spm_make_contrasts(k)
% Make contrasts for one, two or three-way ANOVAs
% FORMAT Con = spm_make_contrasts(k)
%
% k        - vector where the ith entry is the number of levels of factor i
%
% Con      - struct array with fields:
% Con(c).c    - Contrast matrix
%       .name - Name
%
% This function computes contrasts for a generic k(1)-by-k(2)-by-k(3)
% design. It is assumed that the levels of the first factor change slowest.
%
% For one-way ANOVAs set k=L, where L is the number of
% levels, for two-way ANOVAs set k=[L1 L2], for three way set k=[L1 L2 L3]
%
% This function generates (transposed) contrast matrices to test
% average effect, main effect of each factor and interactions.
%__________________________________________________________________________
%
% Reference:
%
% For details of Kronecker operations, see section 5 of
%     http://www.fil.ion.ucl.ac.uk/~wpenny/publications/rik_anova.pdf
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_make_contrasts.m 6157 2014-09-05 18:17:54Z guillaume $


%-Number of factors
%--------------------------------------------------------------------------
nf = numel(k);
k  = [k(:)' 1 1 1];

if nf > 3
    fprintf('spm_make_contrasts not written for %d-way ANOVAs.\n',nf);
    Con = [];
    return
end

%-Get common and differential vectors
%--------------------------------------------------------------------------
C1 = ones(k(1),1);
C2 = ones(k(2),1);
C3 = ones(k(3),1);
C4 = ones(k(4),1);
D1 = -diff(eye(k(1)))';
D2 = -diff(eye(k(2)))';
D3 = -diff(eye(k(3)))';
D4 = -diff(eye(k(4)))';

%-Average effect, main effects and interactions contrast vectors
%--------------------------------------------------------------------------
Con(1).c = kron(kron(kron(C1,C2),C3),C4)';
Con(2).c = kron(kron(kron(D1,C2),C3),C4)';
Con(1).name = 'Average effect of condition';
Con(2).name = 'Main effect of factor 1';

if nf > 1
    Con(3).c = kron(kron(kron(C1,D2),C3),C4)';
    Con(4).c = kron(kron(kron(D1,D2),C3),C4)';
    Con(3).name = 'Main effect of factor 2';
    Con(4).name = 'Interaction, factor 1 x 2';
end

if nf > 2
    Con(5).c = kron(kron(kron(C1,C2),D3),C4)';
    Con(6).c = kron(kron(kron(D1,C2),D3),C4)';
    Con(7).c = kron(kron(kron(C1,D2),D3),C4)';
    Con(8).c = kron(kron(kron(D1,D2),D3),C4)';
    Con(5).name = 'Main effect of factor 3';
    Con(6).name = 'Interaction, factor 1 x 3';
    Con(7).name = 'Interaction, factor 2 x 3';
    Con(8).name = 'Interaction, factor 1 x 2 x 3';
end

if nf > 3
    Con(9).c  = kron(kron(kron(C1,C2),C3),D4)';
    Con(10).c = kron(kron(kron(D1,C2),C3),D4)';
    Con(11).c = kron(kron(kron(C1,D2),C3),D4)';
    Con(12).c = kron(kron(kron(C1,C2),D3),D4)';
    Con(13).c = kron(kron(kron(D1,D2),C3),D4)';
    Con(14).c = kron(kron(kron(D1,C2),D3),D4)';
    Con(15).c = kron(kron(kron(C1,D2),D3),D4)';
    Con(16).c = kron(kron(kron(D1,D2),D3),D4)';
    Con(9).name  = 'Main effect of factor 4';
    Con(10).name = 'Interaction, factor 1 x 4';
    Con(11).name = 'Interaction, factor 2 x 4';
    Con(12).name = 'Interaction, factor 3 x 4';
    Con(13).name = 'Interaction, factor 1 x 2 x 4';
    Con(14).name = 'Interaction, factor 1 x 3 x 4';
    Con(15).name = 'Interaction, factor 2 x 3 x 4';
    Con(16).name = 'Interaction, factor 1 x 2 x 3 x 4';
end
