function [mar] = spm_mar_spectra (mar,freqs,ns,show)
% Get spectral estimates from MAR model
% FORMAT [mar] = spm_mar_spectra (mar,freqs,ns,show)
%
% mar   - MAR data structure (see spm_mar.m)
% freqs - [Nf x 1] vector of frequencies to evaluate spectra at
% ns    - samples per second (default: ns = 2*freqs(end))
% show  - 1 if you wish to plot estimates (default is 0)
%
% The returned mar will have the following fields specified:
%
% .P     [Nf x d x d] Power Spectral Density matrix
% .H     [Nf x d x d] Transfer Function matrix
% .C     [Nf x d x d] Coherence matrix
% .dtf   [Nf x d x d] Kaminski's Directed Transfer Function matrix
% .pve   [Nf x d x d] Geweke's proportion of variance explained
% .gew   [Nf x d x d] Geweke's frequency domain Granger causality
% .pdc   [Nf x d x d] Baccala's Partial Directed Coherence
% .L     [Nf x d x d] Phase matrix
% .f     [Nf x 1] Frequency vector
% .ns    Sample rate
%
% dtf(f,i,j) is the DTF at frequency f from signal j to signal i
% pdc(f,i,j) is the PDC at frequency f from signal j to signal i
% pve(f,i,j) is the proportion of power in signal i at frequency f that can
%            be predicted by signal j. 
% gew(f,i,j) is the Granger casuality from signal j to signal i at frequency f.
%            gew=-log(1-pev)
%
% For DTF and PDC see L. Baccala and K. Sameshima (2001) Biol Cyb 84, 463-474.
% For PVE and GEW see A. Brovelli et al. (2004) PNAS 101(26) 9849-9854.
%
% In addition to the definition of PDC in the above paper, in this
% implementation PDC is also scaled by the observation noise variance
% (Baccala, personal communication).
%
% Also note that PVE and GEW are only valid for d=2 time series
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mar_spectra.m 6198 2014-09-25 10:38:48Z karl $

% options
%--------------------------------------------------------------------------
if nargin < 4  || isempty(show), show = 0; end

% Nyquist
%--------------------------------------------------------------------------
if nargin < 3, ns = 2*freqs(end); end


% dimensions (from mar.lag)
%--------------------------------------------------------------------------
p     = length(mar.lag);
d     = length(mar.lag(1).a);

Nf    = length(freqs);
mar.f = freqs;
w     = 2*pi*freqs/ns;

% precision of innovation
%--------------------------------------------------------------------------
if isfield(mar,'noise_cov');
    noise_cov = mar.noise_cov;
    prec      = diag(diag(inv(noise_cov)));
else
    noise_cov = eye(d,d);
    prec      = eye(d,d);
end

% Get Power Spectral Density matrix and DTF
%==========================================================================
for ff = 1:Nf,
    
  % transfer function (H) and CSD (P)
  %------------------------------------------------------------------------
  af_tmp = eye(d,d);
  for k = 1:p
    af_tmp = af_tmp + mar.lag(k).a*exp(-1i*w(ff)*k);
  end
  iaf_tmp       = inv(af_tmp);
  mar.H(ff,:,:) = iaf_tmp;
  mar.P(ff,:,:) = iaf_tmp * noise_cov * iaf_tmp';
  
  % Get DTF and PDC
  %------------------------------------------------------------------------
  for ii = 1:d
      for j = 1:d
          
          % DTF uses iaf_tmp and normalises wrt rows (to = sink)
          %----------------------------------------------------------------
          mar.dtf(ff,ii,j) = abs(iaf_tmp(ii,j))/sqrt(iaf_tmp(ii,:)*iaf_tmp(ii,:)');
          
          % PDC uses af_tmp and normalises wrt columns (from = source)
          %----------------------------------------------------------------
          mar.pdc(ff,ii,j) = abs(af_tmp(ii,j))/sqrt(abs(af_tmp(:,j)'*prec*af_tmp(:,j)));
      end
  end
end

% Get Coherence and Phase
%--------------------------------------------------------------------------
for k = 1:d
    for j = 1:d
        rkj          = mar.P(:,k,j)./(sqrt(mar.P(:,k,k)).*sqrt(mar.P(:,j,j)));
        mar.C(:,k,j) = abs(rkj);
        mar.L(:,k,j) = atan(imag(rkj)./real(rkj));
    end
end

% Get Geweke's formulation of Granger Causality in the Frequency domain
%--------------------------------------------------------------------------
C     = noise_cov;
for j = 1:d
    for k = 1:d
        rkj = C(j,j)-(C(j,k)^2)/C(k,k);
        sk  = abs(mar.P(:,k,k));
        hkj = abs(mar.H(:,k,j)).^2;
        mar.pve(:,k,j) = rkj*hkj./sk;
        mar.gew(:,k,j) = -log(1 - mar.pve(:,k,j));
    end
end

% Normalise cross spectral density 
%--------------------------------------------------------------------------
mar.P = 2*mar.P/ns;


% plot results if requested
%==========================================================================
if show
    % Plot spectral estimates 
    h=figure;
    set(h,'name','Log Power Spectral Density');
    for k = 1:d
        for j = 1:d
            index=(k-1)*d+j;
            subplot(d,d,index);
            psd=real(mar.P(:,k,j)).^2;
            plot(mar.f,log(psd));
        end
    end
    
    h=figure;
    set(h,'name','Coherence');
    for k = 1:d
        for j = 1:d
            if ~(k==j)
                index=(k-1)*d+j;
                subplot(d,d,index);
                coh=real(mar.C(:,k,j)).^2;
                plot(mar.f,coh);
                axis([min(mar.f) max(mar.f) 0 1]);
            end
        end
    end
    
%     h=figure;
%     set(h,'name','DTF');
%     for k=1:d,
%         for j=1:d,
%             if ~(k==j)
%                 index=(k-1)*d+j;
%                 subplot(d,d,index);
%                 dtf=mar.dtf(:,k,j);
%                 plot(mar.f,dtf);
%             end
%         end
%     end
    
    h=figure;
    set(h,'name','PDC');
    for k = 1:d
        for j = 1:d
            if ~(k==j)
                index=(k-1)*d+j;
                subplot(d,d,index);
                pdc=mar.pdc(:,k,j);
                plot(mar.f,pdc);
                axis([min(mar.f) max(mar.f) 0 1]);
            end
        end
    end
    
%     h=figure;
%     set(h,'name','Phase');
%     for k=1:d,
%         for j=1:d,
%             if ~(k==j)
%                 index=(k-1)*d+j;
%                 subplot(d,d,index);
%                 ph=mar.L(:,k,j);
%                 plot(mar.f,ph);
%             end
%         end
%     end
%     
%     h=figure;
%     set(h,'name','Delay/ms');
%     for k=1:d,
%         for j=1:d,
%             if ~(k==j)
%                 index=(k-1)*d+j;
%                 subplot(d,d,index);
%                 ph=mar.L(:,k,j);
%                 plot(mar.f,1000*ph'./(2*pi*mar.f));
%             end
%         end
%     end
    
end