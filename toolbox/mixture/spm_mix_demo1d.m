
% Demonstrate use of spm_mix on 1D data

% Number of data points
N=100;

% Single cluster data
disp('Data set 1: single cluster data');
y=randn(N,1);

logev=[];
for m=1:5,
    disp(sprintf('Fitting mixture model with %d components',m));
    vbmix=spm_mix(y,m,0);
    logev=[logev; vbmix.fm];
end
logev=logev-min(logev);

figure
subplot(2,1,1);
hist(y);
title('Data set 1: Data histogram');
subplot(2,1,2);
bar(logev);
xlabel('Number of mixture components');
ylabel('Log Evidence');
drawnow;

% Three cluster data
disp('Data set 2: three cluster data');
mix.m=3;
mix.state(1).prior=0.3;
mix.state(1).m=-5;
mix.state(1).C=1;

mix.state(2).prior=0.3;
mix.state(2).m=0;
mix.state(2).C=1;

mix.state(3).prior=0.4;
mix.state(3).m=5;
mix.state(3).C=1;

y=spm_samp_mix(mix,N);
logev=[];
for m=1:5,
    disp(sprintf('Fitting mixture model with %d components',m));
    vbmix=spm_mix(y,m,0);
    logev=[logev; vbmix.fm];
end
logev=logev-min(logev);

figure
subplot(2,1,1);
hist(y);
title('Data set 2: Data histogram');
subplot(2,1,2);
bar(logev);
xlabel('Number of mixture components');
ylabel('Log Evidence');



