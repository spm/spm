
close all
clear all

% Generate data from PCA model y_n = W x_n + mu + e_n
% where W is d*q - first example in Bishop's VPCA paper

N=100;
d=10;
q=d-1;

% Generate orthogonal directions:
W=randn(d,q);
W=orth(W);

% Generate sources
x=randn(q,N);
sd_x=diag([5,4,3,2,1,1,1,1,1]);
x=sd_x(1:q,1:q)*x;

% Generate sensor data
e=randn(d,N);

mu=ones(d,1)*ones(1,N);

%e=zeros(d,N);
t=W*x+mu+e;

% Get pca solution
pca=spm_vpca(t);

figure; imagesc(pca.M_w); colormap gray; title('Bayes estimate');
figure; imagesc(pca.ml.W); colormap gray; title('ML estimate');

figure
plot(pca.Fm_evol);
xlabel('Iterations');
ylabel('Neg. Free Energy');

figure
plot(pca.ml.lambda);
title('Eigenspectrum');

figure
plot(pca.mean_alpha);
title('Prior precision of factors');

disp(' ');
disp('To recover eg 4 hidden sources')
disp('project data, t, onto first 4 columns of factor matrix:');
disp('W=pca.M_w(:,1:4);xhat=W''*t;');
W=pca.M_w(:,1:4);
xhat=W'*t;

disp(' ');
disp('A Bayesian estimate of the data covariance matrix');
disp('is given by:');
disp('obs_noise_var=(1/pca.mean_tau); CBayes=W*W''+obs_noise_var*eye(pca.d)');
obs_noise_var=(1/pca.mean_tau);
CBayes=W*W'+obs_noise_var*eye(pca.d);

C=cov(t');
figure
subplot(2,2,1);
imagesc(C);
colormap gray
colorbar
title('Data Covariance using Cov');
subplot(2,2,3);
imagesc(CBayes);
colormap gray
colorbar
title('Data Covariance using Bayes VPCA');
subplot(2,2,2);
imagesc(C-CBayes);
colormap gray
colorbar
title('Difference');

