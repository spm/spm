
% Demo for MAR spectral estimation

% This simulation is similar to Cassidy and Brown. Spectral Phase
% Estimates in the setting of multidirectional coupling
% J Neurosci Methods. 2003 Aug 1-15;37(3):299.

close all
noise_dev=0.01;

Nsines=100;
f=sqrt(0.2)*randn(Nsines,1)+20;
secs=50;
ns=100;
t=[1/ns:1/ns:secs]';
N=length(t);

y=zeros(N,1);
for n=1:Nsines,
    y=y+sin(2*pi*f(n)*t);
end
y=y/std(y); % Rescale to unit variance

delay=50; % ms delay
delay_in_samples=ns*delay/1000;

y1=y+noise_dev*randn(N,1);
y2=[y1(delay_in_samples:end);zeros(delay_in_samples-1,1)];
y2=y2+noise_dev*randn(N,1);
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


