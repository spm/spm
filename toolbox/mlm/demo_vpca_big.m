
% Load Alan
x=imread('Alan.jpg','jpg');

% Crop Alan
startx=100;
starty=300;

N=256;
xg=double(x(startx:startx+N-1,starty:starty+N-1,2));
xg=xg-mean(mean(xg));

% Maximum latent space dimension
q=100;
pca=spm_vpca(xg,q);

figure; imagesc(pca.M_w); colormap gray; title('Bayes estimate');
figure; imagesc(pca.ml.W(:,1:q)); colormap gray; title('ML estimate');

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
