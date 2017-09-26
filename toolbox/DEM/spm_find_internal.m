function nj = spm_find_internal(z,J)
% FORMAT nj = spm_find_internal(z,J)
% finds indices of internal states (that do not contribute to slow modes)
%__________________________________________________________________________
[u,s] = eig(full(J),'nobalance');
[s,j] = sort(real(diag(s)),'descend');
u     = u(:,j);
v     = pinv(u);
nu    = sum(s > (max(s) - 4));
u     = u(:,1:nu);
v     = v(1:nu,:);
s     = real(s(1:nu));
nz    = numel(z);
for i = 1:nz
        dJdj            = spm_zeros(J);
        dJdj(:   ,z{i}) = 1;
        dJdj(z{i},z{i}) = 0;
        nj(i)           = exp(s')*diag(abs(v*dJdj*u));
end
