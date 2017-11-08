function [C] = spm_Ce(t,v,a)
% Error covariance constraints (for serially correlated data)
% FORMAT [C] = spm_Ce(v,a)
% FORMAT [C] = spm_Ce('ar',v,a)
% v  - (1 x n) v(i) = number of observations for i-th block
% a  - AR coefficient expansion point  [Default: a = []]
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
% FORMAT [C] = spm_Ce('fast',v,tr)
% v  - (1 x n) v(i) = number of observations for i-th block
% tr - repetition time
%
% See also: spm_Q.m
%__________________________________________________________________________
% Copyright (C) 2000-2017 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_Ce.m 7203 2017-11-08 12:49:15Z guillaume $
 
 
%-Defaults (and backward compatibility with spm_Ce(v,a) == spm_Ce('ar',v,a))
%--------------------------------------------------------------------------
if ~ischar(t)
    if nargin > 1, a = v; else a = []; end
    v = t;
    t = 'ar';
else
    if nargin == 2, a = []; end
end

%-Error covariance constraints
%--------------------------------------------------------------------------
switch lower(t)
    
    case 'ar'
        
        %-Create block diagonal components
        %------------------------------------------------------------------
        C    = {};
        l    = length(v);
        n    = sum(v);
        k    = 0;
        if l > 1
            for i = 1:l
                dCda  = spm_Ce(v(i),a);
                for j = 1:length(dCda)
                    [x,y,q]    = find(dCda{j});
                    x          = x    + k;
                    y          = y    + k;
                    C{end + 1} = sparse(x,y,q,n,n);
                end
                k     = v(i) + k;
            end
        else
            
            %-dCda
            %--------------------------------------------------------------
            if ~isempty(a)
                Q    = spm_Q(a,v);
                dQda = spm_diff('spm_Q',a,v,1);
                C{1} = Q - dQda{1}*a;
                C{2} = Q + dQda{1}*a;
            else
                C{1} = speye(v,v);
            end
            
        end
        
    case 'fast'
        
        dt = a;
        C  = {};
        n  = sum(v);
        k  = 0;
        for m=1:length(v)
            T     = (0:(v(m) - 1))*dt;
            d     = 2.^(floor(log2(dt/4)):log2(64));
            for i = 1:min(6,length(d))
                for j = 0:2
                    QQ = toeplitz((T.^j).*exp(-T/d(i)));
                    [x,y,q] = find(QQ);
                    x = x + k;
                    y = y + k;
                    C{end + 1} = sparse(x,y,q,n,n);
                end
            end
            k = k + v(m);
        end
        
    otherwise
        error('Unknown error covariance constraints.');
end
