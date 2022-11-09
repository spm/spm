function [p] = spm_ar_freq (ar, freq, fs)
% Compute spectra from AR coefficients
% FORMAT [p] = spm_ar_freq (ar, freq, fs)
%
% ar    AR model data structure (see spm_ar.m)
% freq  [Nf x 1] vector containing list of frequencies
% fs    sample rate
%
% p     [Nf x 1] vector containing power estimates
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


Nf=length(freq);
Np=length(ar.a_mean);
T=1/fs;
freq=freq(:);
    
if Np < Nf
    % Usual case - fewer AR coeffs than freqs to evaluate at
    exponent=-i*2*pi*freq*T;
    
    denom=ones(Nf,1);
    for j=1:Np,
        denom=denom+ar.a_mean(j)*exp(j*exponent);
    end
    denom=abs(denom).^2;
    p=1./denom;
else
    p_order=[1:Np]';
    for f=1:Nf,
        denom=abs(1+sum(ar.a_mean.*exp(-i*2*pi*freq(f)*T*p_order)));
        p(f)=1/abs(denom)^2;
    end
end

noise_dev=sqrt(1/ar.mean_beta);
p=p*noise_dev*T;
    