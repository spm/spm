function [M,Mi,i] = spm_meanby(Y,I)
% Means of data in columns by group
% FORMAT [M,Mi,i] = spm_meanby(Y,I)
% Y  - Data matrix, data in columns. (Row vector Y also accepted.)
% I  - Column of indicator vectors, indicating group membership of rows of Y
%    - Multi-column I are treated as multiple factors to be interacted, 
%      and means are computed within each unique combination of the factor levels
% M  - Matrix of same size as Y, with observations replaced by the
%     appropriate group mean
% Mi - Mean for observations in each group, one column for each column of Y,
%      one row for each group (or unique factor level combination)
% i  - Group indicator values corresponding to rows of Mi
%_______________________________________________________________________
%
% spm_meanby computes means for grouped data presented as columns of
% data with a vector of group indicators
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_meanby.m 1143 2008-02-07 19:33:33Z spm $



%-Check arguments
%=======================================================================
if nargin<0, error('insufficient arguments'), end
if  size(Y,1)==1, bT=1; Y=Y'; else, bT=0; end
if nargin<2, I=ones(size(Y,1),1); end
if  size(I,1)~=size(Y,1), I=I'; end
if  size(I,1)~=size(Y,1), error('row dimension mismatch between inputs'), end


%-Computation - using spm_DesMtx
%=======================================================================
[X,null,i,idx,jdx] = spm_DesMtx(I);

i  = i';
Mi = (X'*Y)./repmat(sum(X,1)',1,size(Y,2));
M  = Mi(jdx,:);

if bT, M=M'; end
