function [M,Mi,V,Vi,iVals]=spm_meanby(X,i)
% Mean (& variances) of data in columns by group
% FORMAT [M,Mi,V,Vi,iVals]=spm_meanby(X,i)
% X	- Data matrix, with data in columns (row vector X is transposed)
% i	- Column of indicator vectors, indicating group membership of
%	  observations in rows of X
% Mi	- Mean for observations in each group, for each column
% Vi    - (Sample) variance for observations in each group
% iVals	- Group indicator values corresponding to rows of Mi
% M	- Matrix of same size as X, with observations replaced by the
%	  appropriate group mean
% V     - Matrix of same size as X, with observations replaced by the
%	  appropriate group sample variance
%_______________________________________________________________________
% spm_meanby computes means for grouped data presented as columns of
% data with a vector of group indicators
%_______________________________________________________________________
% %E% Andrew Holmes %W%

if nargin<2, i=ones(size(X,1),1); end
M = zeros(size(X));
V = zeros(size(X));
if size(X,1)==1, X=X'; end
if min(size(i))~=1, error('vector i'), end
i=i(:);

if ~all(floor(i)==ceil(i)), error('Non-integer indicator vector'), end
if size(X,1)~=length(i), error('X must have row dimension the length of i'), end

%-Computation
%=======================================================================

%-Sort out unique indicator levels
tmp    = sort(i);
iVals  = tmp([1;diff(tmp)>0]);

%-Compute indicator matrix
I       = zeros(size(X,1),length(iVals));
for p_i = 1:length(iVals)
	I(:,p_i)=i==iVals(p_i);
end

%-Work out number of elements by index
Ni = sum(I)'*ones(1,size(X,2));

%-Work out sums by index
Si = I'*X;

%-Work out means by index
Mi = Si ./ Ni;

%-Combine into matrix same size as X:
% (Effectively replacing each observation by its group mean)
M(:) = Mi(i,:);

%-Compute group variances if requested
%-----------------------------------------------------------------------
if nargout<3, return, end

%-Work out sum of squares by index
SSi = I'*X.^2;

%-Work out variances by index
Vi = ( SSi - Si.^2./Ni ) ./ (Ni -1);

%-Combine into matrix same size as X:
V(:) = Vi(i,:);
