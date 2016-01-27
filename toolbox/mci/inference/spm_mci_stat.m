function [pstat,mu,nse,batch] = spm_mci_stat (post,nbatch,method)
% Test MCMC samples for stationarity
% FORMAT [pstat,mu,nse,batch] = spm_mci_stat (post,nbatch,method)
%
% post      posterior distribution
% nbatch    number of batches (last batch contains last half of samples)
% method    'geweke', 'ar1' (default) or 'geyer'
%
% pstat     p-value for batch mean being different to final batch mean
% mu        batch mean (of energy)
% nse       batch numeric standard error (of energy)
% batch     (n).ind, (n).N indices and number of samples in nth batch
%
% This routine is based on Geweke 1992. But we also allow estimates 
% of the Numeric Standard Error (NSE) to be estimated using an AR1 model 
% or Geyer's method
% 
% J. Geweke (1992) Evaluating the accuracy of sampling-base approaches to 
% the calculation of posterior moments. Bayesian Statistics 4, OUP.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_stat.m 6697 2016-01-27 14:57:28Z spm $

try, meth=method; catch, meth='ar1'; end
try, Nb=nbatch; catch, Nb=6; end
Ns=size(post.P,2);

% Define batches
bl=floor(0.5*Ns/(Nb-1));
boffset=0;
for b=1:Nb-1,
    batch(b).ind=[1:bl]+boffset;
    boffset=boffset+bl;
end
batch(Nb).ind=boffset+1:Ns;

% Get Means and NSE for each batch
for b=1:Nb,
    ind=batch(b).ind;
    batch(b).N=length(ind);
    E=post.E(ind);
    mu(b)=mean(E);
    N(b)=length(E);
    switch meth
        
        case 'geweke',
            [spec,f] = spectrum(E-mean(E));
            s0 = spec(1);
            nse(b)=sqrt(s0);
        
        case 'geyer',
            % Geyer's method
            ess=spm_mci_ess(E);
            nse(b)=std(E)*sqrt(N(b))/sqrt(ess);
            
        case 'ar1',
            % AR(1) method
            ar=spm_ar(E,1);
            nse(b)=sqrt(1/(ar.mean_beta));
    end
end

% Are first Nb-1 batch means different to final batch mean ?
for b=1:Nb-1,
    % Geweke's test statistic is a z-score
    d = mu(b)-mu(Nb);
    s2 = nse(b)^2+nse(Nb)^2;
    pbigger = 1-spm_Ncdf(abs(d),0,s2);
    psmaller = spm_Ncdf(-abs(d),0,s2);
    pstat(b)=pbigger+psmaller;
end

end


% ------------------------------------------------------------

function [y,f]=spectrum(x,nfft,nw)
%SPECTRUM Power spectral density using Hanning window
%  [y,f]=spectrum(x,nfft,nw) 

if nargin < 2 || isempty(nfft)
  nfft = min(length(x),256);
end
if nargin < 3 || isempty(nw)
  nw = fix(nfft/4);
end
noverlap = fix(nw/2);

% Hanning window
w = .5*(1 - cos(2*pi*(1:nw)'/(nw+1)));

n = length(x);
if n < nw
    x(nw)=0;  n=nw;
end
x = x(:);

k = fix((n-noverlap)/(nw-noverlap)); % no of windows
index = 1:nw;
kmu = k*norm(w)^2; % Normalizing scale factor
y = zeros(nfft,1);
for i=1:k
  xw = w.*x(index);
  index = index + (nw - noverlap);
  Xx = abs(fft(xw,nfft)).^2;
  y = y + Xx;
end

y = y*(1/kmu); % normalize

n2 = floor(nfft/2);
y  = y(1:n2);
f  = 1./n*(0:(n2-1));

end
