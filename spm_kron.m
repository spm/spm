function K = spm_kron(A,B)
% Kronecker tensor product with sparse outputs 
% FORMAT K = spm_kron(A,B)
%        K = spm_kron(A)
%
%   KRON(X,Y) is the Kronecker tensor product of X and Y.
%   The result is a large matrix formed by taking all possible
%   products between the elements of X and those of Y.   For
%   example, if X is 2 by 3, then KRON(X,Y) is
%
%      [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
%        X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
%
%   Class support for inputs X,Y:
%      float: double, single
% 
% When called with a single cell array input, the tensor product
% is formed recursively 

%   Previous versions by Paul Fackler, North Carolina State,
%   and Jordan Rosenthal, Georgia Tech.
%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 7760 $ $Date: 2004/06/25 18:52:18 $
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_kron.m 7760 2019-12-29 17:45:58Z karl $


% Deal with cell arrays
%--------------------------------------------------------------------------
if iscell(A)
    K = 1;
    for i = 1:numel(A)
        K = spm_kron(A{i},K);
    end
    return
end

% Kronecker tensor product
%--------------------------------------------------------------------------
[ma,na] = size(A);
[mb,nb] = size(B);

[ia,ja,sa] = find(A);
[ib,jb,sb] = find(B);
ia = ia(:); ja = ja(:); sa = sa(:);
ib = ib(:); jb = jb(:); sb = sb(:);
ka = ones(size(sa));
kb = ones(size(sb));
t  = mb*(ia-1)';
ik = t(kb,:)+ib(:,ka);
t  = nb*(ja-1)';
jk = t(kb,:)+jb(:,ka);
K  = sparse(ik,jk,sb*sa.',ma*mb,na*nb);


