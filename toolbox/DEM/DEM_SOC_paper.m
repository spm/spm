
%--------------------------------------------------------------------------
clear all
n  = 4;
m  = 8;
dv = 1e-4;
for i = 1:128
    C       = spm_orth(randn(n*m,n*m),1);
    v       = sort(randn(n*m,1).^2,1)*512;
    v       = v - min(v) + exp(-6);
    D       = kron(spm_speye(n,n,1),eye(m,m));
    
    s(i)  = max(real(eig(full(D - C*diag(v)*C'))));
    e(i)  = v(1);
    v(1)  = v(1) + dv;
    ds(i) = (s(i) - max(real(eig(full(D - C*diag(v)*C')))))/dv;
    
end

subplot(2,1,1)
plot(s,ds,'.')

subplot(2,1,2)
bar(-v)