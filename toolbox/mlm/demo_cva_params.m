
clear all
close all

d1=3;
d2=5;
N=30;

%m=min([d1,d2])
m=1;

% Generate true factor matrices
W1=10*randn(d1,m);
W2=10*randn(d2,m);

% Observation noise covariance
sig=0.01;

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

CVA = spm_cva_prob (X1,X2);

disp('True');
abs(W1)
disp('Estimated');
abs(CVA.W1)

disp('True');
abs(W2)
disp('Estimated');
abs(CVA.W2)
