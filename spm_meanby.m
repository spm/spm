function [M,Mi,iVals]=spm_meanby(X,i)
% Mean of data in columns by group
% FORMAT [M,Mi,iVals]=spm_meanby(X,i)
% X	- Data matrix, with data in columns (row vector X is transposed)
% i	- Column of indicatir vectors, indicating group membership of
%	  observations in rows of X
% Mi	- Mean for observations in each group, for each column
% iVals	- Group indicator values corresponding to rows of Mi
% M	- Matrix of same size as X, with observations replaced by the
%	  appropriate group mean
%_______________________________________________________________________
% spm_meanby computes means for grouped data presented as columns of
% data with a vector of group indicators
%_______________________________________________________________________
% %W% Andrew Holmes %E%

if nargin<2, i=ones(size(X,1),1); end
if size(X,1)==1, X=X'; end
if min(size(i))~=1, error('vector i'), end
i=i(:);

if ~all(floor(i)==ceil(i)), error('Non-integer indicator vector'), end
if size(X,1)~=length(i), error('X must have row dimension the length of i'), end

%-Sort out unique indicator levels
tmp    = sort(i);
iVals  = tmp([1;diff(tmp)>0]);

%-Compute indicator matrix
I       = zeros(size(X,1),length(iVals));
for p_i = 1:length(iVals)
	I(:,p_i)=i==iVals(p_i);
end

%-Work out means by index
Mi = I'*X ./ sum(I)'*ones(1,size(X,2));

%-Combine into matrix of group means for each observation
M = Mi(i,:);
