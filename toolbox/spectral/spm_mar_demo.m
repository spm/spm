% Generate MAR(2) and fit MAR model

disp('Generating data from known MAR(2) model');
d=2;
p=2;
T=100;
w=[0;0];

% Coeffs at lag 1
A1 = [ 0.4   1.2;   0.3   0.7 ];
% Coeffs at lag 2
A2 = [ 0.35 -0.3;  -0.4  -0.5 ];
A = [ A1 A2 ];

C = [ 1.00  0.50;   0.50  1.50 ];
lambda_true=inv(C);

%  Generate observations
x = spm_mar_gen (w, A, C, T);

logev=[];
for m=1:5,
    disp(sprintf('Fitting MAR model with %d components',m));
    mar=spm_mar(x,m);
    logev=[logev; mar.fm];
end
logev=logev-min(logev);

figure
subplot(2,1,1);
plot(x);
title('Bivariate time series from MAR(2) model');
subplot(2,1,2);
bar(logev);
xlabel('Number of time lags');
ylabel('Log Evidence');


% Specify prior - this is optional. 
% spm_mar.m runs without the prior being set.
prior=spm_mar_prior(d,p,'global');
[mar,y,y_pred]=spm_mar(x,2,prior);

disp(' ');
disp('Estimates from fitting MAR(2) model');
disp('Lag 1');
disp('True coefficients');
disp(A1);
disp('Estimated coefficients');
disp(-mar.lag(1).a)

disp('Lag 2');
disp('True coefficients');
disp(A2);
disp('Estimated coefficients');
disp(-mar.lag(2).a)

disp('Noise covariance');
disp('True:');
disp(C);
disp('Estimated:');
disp(mar.noise_cov);




