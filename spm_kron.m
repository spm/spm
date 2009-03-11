function K = spm_kron(A,B)
% as for kron, but forces sparse outputs
%KRON   Kronecker tensor product.
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

%   Previous versions by Paul Fackler, North Carolina State,
%   and Jordan Rosenthal, Georgia Tech.
%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 2863 $ $Date: 2004/06/25 18:52:18 $
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_kron.m 2863 2009-03-11 20:25:33Z guillaume $

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


