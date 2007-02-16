
N=200;
fs=100;
t=[1:1:N]/fs;
freqs=[1:45];

% Single sunusoid
x=sin(2*pi*10*t);
p = spm_wavspec (x,freqs,fs,1);
figure
subplot(2,1,1);
plot(t,x);
title('One sinusoid');
subplot(2,1,2);
imagesc(p);
xlabel('Time/Samples');
ylabel('Frequency');

% Two sinusoids
x=x+sin(2*pi*38*t);
p = spm_wavspec (x,freqs,fs);
figure
subplot(2,1,1);
plot(t,x);
title('Two sinusoids');
subplot(2,1,2);
imagesc(p);
xlabel('Time/Samples');
ylabel('Frequency');

% Chirp
load chirp
p = spm_wavspec (x,freqs,fs);
figure
subplot(2,1,1);
plot(t,x);
title('Chirp');
subplot(2,1,2);
imagesc(p);
xlabel('Time/Samples');
ylabel('Frequency');