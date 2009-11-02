function [C] = spm_Ce(v,a)
% return error covariance constraints (for serially correlated data)
% FORMAT [C] = spm_Ce(v,a)
% v  - (1 x l) v(i) = number of observations for i-th block
% a  - AR coefficient expansion point  (default a = [])
% 
% a  = [] (default) - block diagonal identity matrices specified by v:
%
%   C{i}  = blkdiag( zeros(v(1),v(1)),...,AR(0),...,zeros(v(end),v(end)))
%   AR(0) = eye(v(i),v(i))
%
% otherwise:
%
%   C{i}     = AR(a) - a*dAR(a)/da;
%   C{i + 1} = AR(a) + a*dAR(a)/da;
%
% See also: spm_Q.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_Ce.m 3527 2009-11-02 20:27:13Z karl $
 
 
% defaults
%--------------------------------------------------------------------------
if nargin == 1, a = []; end
 
% create block diagonal components
%--------------------------------------------------------------------------
C    = {};
l    = length(v);
n    = sum(v);
k    = 0;
if l > 1
    for i = 1:l
        dCda  = spm_Ce(v(i),a);
        for j = 1:length(dCda)
            [x y q]    = find(dCda{j});
            x          = x    + k;
            y          = y    + k;
            C{end + 1} = sparse(x,y,q,n,n);
        end
        k     = v(i) + k;
    end
else
    
    % dCda
    %----------------------------------------------------------------------
    if ~isempty(a)
        Q    = spm_Q(a,v);
        dQda = spm_diff('spm_Q',a,v,1);
        C{1} = Q - dQda{1}*a;
        C{2} = Q + dQda{1}*a;
    else
        C{1} = speye(v,v);
    end
 
end
