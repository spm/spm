
N=100;

% True source dimension is 3
a=[9 7 1;2 4 7;1 3 8; 3 2 1; 1 1 0.5; 7 2 8];
[d,M]=size(a);
s=randn(M,N);
obs_noise=0.1;
X=a*s+sqrt(obs_noise)*randn(d,N);

[p_opt,log_ev,lambda]=spm_pca_order(X);
subplot(2,1,1);
plot(lambda);
title('Eigenspectrum');
subplot(2,1,2);
plot(log_ev);
title('Log Evidence');
disp(sprintf('Estimated number of sources is %d',p_opt));

