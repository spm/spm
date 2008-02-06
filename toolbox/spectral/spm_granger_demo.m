% Granger causality demo based on spm_mar

disp('Generating six sinewaves in two blocks of 3');
disp('The signals are dependent within each block.');
disp(' ');

secs=1;
ns=250;
t=[1/ns:1/ns:secs]';
d=6;
f1=10;
clear x
dev=1*ones(1,6);
y=sin(2*pi*f1*t);
y2=sin(2*pi*12.5*t);
x(:,1)=y+dev(1)*randn(size(t));
for i=2:3,
  x(:,i)=y+dev(i)*randn(size(t));
end
for i=4:6,
  x(:,i)=y2+dev(i)*randn(size(t));
end
for i=1:6,
    x(:,i)=x(:,i)/std(x(:,i));
    x(:,i)=x(:,i)-mean(x(:,i));
end

disp('Estimating order of MAR model');
logev=[];
for m=1:5,
    disp(sprintf('Fitting MAR model with %d components',m));
    mar=spm_mar(x,m);
    logev=[logev; mar.fm];
end
logev=logev-min(logev);
figure
bar(logev);
xlabel('Number of time lags');
ylabel('Log Evidence');

[tmp, p_sel]=max(logev);

disp(sprintf('Using MAR(%d) model ..',p_sel));
[mar,y,y_pred]=spm_mar(x,p_sel);

[G,Psig] = spm_granger (mar);

disp(' ');
disp('True causality matrix');
[ones(3,3),zeros(3,3);zeros(3,3),ones(3,3)]

disp('Granger probability matrix:');
Peffect=ones(6,6)-Psig;
Peffect
disp('where ijth entry is our belief that time series i Granger causes j');
figure
imagesc(Peffect);
title('Granger probability matrix');
colormap gray
colorbar

disp(' ');
disp('Inferred Granger causality matrix:');
disp('This is Granger Prob matrix thresholded at 0.95');
G





