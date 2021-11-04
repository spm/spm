function p = spm_Q_perm(Q)
% returns a cell of permutation indices for separating matrices
% FOTMAT p = spm_Q_perm(Q);
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% John Ashburner & Karl Friston
% $Id: spm_Q_perm.m 8183 2021-11-04 15:25:19Z guillaume $
 
 
% add cell arrays of components
%--------------------------------------------------------------------------
A   = 0;
if iscell(Q)
    for i = 1:length(Q)
        A = A + ~~Q{i};
    end
else
    A = ~~Q;
end
 
% find partition
%--------------------------------------------------------------------------
A = A + A';
p = {};
i = find(any(A,1),1);
while ~isempty(i)
    k     = [];
    while length(k) < length(i)
        k = i;
        j = find(any(A(:,i),2));
        i = j;
    end
    p{end + 1} = i;
    A(i,i)     = 0;
    i          = find(any(A,1),1);
end
