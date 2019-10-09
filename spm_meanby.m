function [M,Mi,i] = spm_meanby(Y,I)
% Means of data in columns by group
% FORMAT [M,Mi,i] = spm_meanby(Y,I)
% Y  - Data matrix, data in columns. (Row vector Y also accepted.)
% I  - Column of indicator vectors, indicating group membership of rows of Y
%    - Multi-column I are treated as multiple factors to be interacted, and
%      means are computed within each unique combination of the factor levels
% M  - Matrix of same size as Y, with observations replaced by the
%      appropriate group mean
% Mi - Mean for observations in each group, one column for each column of Y,
%      one row for each group (or unique factor level combination)
% i  - Group indicator values corresponding to rows of Mi
%__________________________________________________________________________
%
% spm_meanby computes means for grouped data presented as columns of data
% with a vector of group indicators.
%__________________________________________________________________________
% Copyright (C) 2008-2019 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_meanby.m 7675 2019-10-09 10:37:43Z guillaume $


%-Check arguments
%--------------------------------------------------------------------------
if size(Y,1) == 1, bT = true; Y = Y'; else, bT = false; end
if nargin < 2, I = ones(size(Y,1),1); end
if size(I,1) ~= size(Y,1)
    if size(Y,1) == size(I,2)
        I = I';
    else
        str = sprintf('\nData: [%d %d] and Groups: [%d %d]',size(Y),size(I));
        error(['Row dimension mismatch between inputs.' str]);
    end
end


%-Computation - using spm_DesMtx
%--------------------------------------------------------------------------
[X,null,i,idx,jdx] = spm_DesMtx(I);

i  = i';
Mi = (X'*Y)./repmat(sum(X,1)',1,size(Y,2));
M  = Mi(jdx,:);

if bT, M = M'; end
