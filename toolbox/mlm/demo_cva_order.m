
clear all
close all

d1=3;
d2=5;
N=100;

disp('Model order selection');

%m=min([d1,d2])
m=3;

% Generate true factor matrices
W1=10*randn(d1,m);
W2=10*randn(d2,m);

% Observation noise covariance
sig=1;
E1=sig*randn(d1,N);
E2=sig*randn(d2,N);

if m==0
    X1=E1;
    X2=E2;
else
    Z=randn(m,N);
    X1=W1*Z+E1;
    X2=W2*Z+E2;
end

for i=1:4,
    CVA = spm_cva_prob (X1,X2,i-1);
    L(i)=CVA.L;
    bic(i)=CVA.bic;
    aic(i)=CVA.aic;
end
figure
plot([0:3],L);
hold on
plot([0:3],bic,'r');
plot([0:3],aic,'g');
xlabel('Number of Canonical Vectors');
legend('LogLike','BIC','AIC');


