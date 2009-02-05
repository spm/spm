function M = spm_eeg_morlet(Rtf, ST, f, ff)
% Generate Morlet wavelets
% FORMAT M = spm_eeg_morlet(Rtf, ST, f, ff)
% 
% Rtf - 'wavelet factor', see [1]
% ST  - sample time [ms]
% f   - vector of frequencies [Hz]
% ff  - frequency to fix Gaussian envelope (sigma = Rtf/(2*pi*ff))
%       Default is ff = f, ie.e, a Morlet transform
%       NB: FWHM = sqrt(8*log(2))*sigma_t;
%
% M   - cell vector, where each element contains the filter for each
%       frequency in f
%__________________________________________________________________________
% 
% spm_eeg_morlet generates morlet wavelets for specified frequencies f with
% a specified ratio Rtf, see [1], for sample time ST (ms). One obtains the
% wavelet coefficients by convolution of a data vector with the kernels in
% M. See spm_eeg_tf how one obtains instantaneous power and phase estimates
% from the wavelet coefficients.
%
% [1] C. Tallon-Baudry, O. Bertrand, F. Peronnet and J. Pernier, 1998.
% Induced \gamma-Band Activity during the Delay of a Visual Short-term
% memory Task in Humans. The Journal of Neuroscience (18): 4244-4254.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_morlet.m 2696 2009-02-05 20:29:48Z guillaume $


M = {};
scaling = [];
ctb = [];
for f0 = f
    
    % fixed or scale-dependent window
    %----------------------------------------------------------------------
    if nargin < 4
        sigma_t = Rtf/(2*pi*f0);
    else
        sigma_t = Rtf/(2*pi*ff);
    end
    
    % this scaling factor is proportional to (Tallon-Baudry 98): 
    % (sigma_t*sqrt(pi))^(-1/2);
    %----------------------------------------------------------------------
    t = [0:ST*0.001:5*sigma_t];
    t = [-t(end:-1:2) t];
    M{end+1} = exp(-t.^2/(2*sigma_t^2)) .* exp(2 * 1i * pi * f0 *t);    
    M{end}   = M{end}./(sqrt(0.5*sum(real(M{end}).^2 + imag(M{end}).^2)));
end
