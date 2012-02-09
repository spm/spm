function M = spm_karcher(A)
% Compute Karcher mean
N = size(A,3);
M = eye(size(A,1),size(A,2));

for iter = 1:1024,
    S = zeros(size(M));
    for i=1:N,
        L = real(logm(M\A(:,:,i)));
        S = S + L;
    end
    S = S/N;
    M = M*expm(S);
    %imagesc(M); drawnow
    %fprintf('%d\t%g\n', iter,sum(S(:).^2));
    if sum(S(:).^2)<1e-20,
        break;
    end
end

