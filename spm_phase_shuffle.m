function [y] = spm_phase_shuffle(x,n)
% phase-shuffling of a vector
% FORMAT [y] = spm_phase_shuffle(x,[n])
% x   - data matrix (time-series in columns)
% n   - optional window length for windowed shuffling
%__________________________________________________________________________
% Copyright (C) 2007-2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_phase_shuffle.m 6654 2015-12-22 12:55:36Z spm $
 
 
try
    
    % randomise phase - WFT
    %----------------------------------------------------------------------
    k     = 1:fix(n/2);
    for i = 1:size(x,2);
        C      = spm_wft(x(:,i),k,n);
        W      = abs(C).*exp(1i*angle(C(randperm(size(C,1)),:)));
        y(:,i) = spm_iwft(W,k,n)';
    end
    
catch
    
    % randomise phase - FT
    %----------------------------------------------------------------------
    n              = size(x,1);
    s              = fft(x);
    i              = 2:ceil(n/2);
    r              = rand(length(i),size(x,2))*2*pi - pi;
    p              = zeros(n,size(x,2));
    p(i,:)         =  r;
    p(n - i + 2,:) = -r;
    s              = abs(s).*exp(1i*p);
    y              = real(ifft(s));
 
end
