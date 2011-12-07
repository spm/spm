
% Show basic use of spm_mlm_bayes

clear all
close all

N=100;
d=2;
p=3;

y=randn(N,d);
x=randn(N,p);
verbose=1;

mlm = spm_mlm_bayes (y,x,'input',verbose);

figure
imagesc(mlm.wmean);
colorbar
ylabel('Inputs');
xlabel('Outputs');
colormap(gray);


disp(sprintf('Model evidence = %1.2f', mlm.fm));