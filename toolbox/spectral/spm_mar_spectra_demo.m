
% Demo for MAR spectral estimation

f=10;
secs=2;
ns=100;
t=[1/ns:1/ns:secs]';
N=length(t);
noise_dev=0.1;;
noise1=noise_dev*randn(N,1);
noise2=noise_dev*randn(N,1);

delay=50; % ms delay
delay_in_samples=ns*delay/1000;

x1=sin(2*pi*f*t);
y1=x1+noise1;
dy1=[y1(delay_in_samples:end);noise_dev*randn(delay_in_samples-1,1)];
y2=dy1+noise2;
y=[y1,y2];

h=figure;
set(h,'name','Data');
plot(t,y1);
hold on
plot(t,y2+3);
xlabel('Seconds');

p=10; % order of MAR model
freqs=[0.5:0.5:32];
mar = spm_mar(y,p);
mar = spm_mar_spectra (mar,freqs,ns,1);


