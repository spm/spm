function [mar] = spm_mar_spectra (mar,freqs,ns,show)
% Get spectral estimates from MAR model
% FORMAT [mar] = spm_mar_spectra (mar,freqs,ns,show)
%
% mar   MAR data structure (see spm_mar.m)
% freqs [Nf x 1] vector of frequencies to evaluate spectra at
% ns    samples per second
% show  1 if you wish to plot estimates (default is 0)
%
% The returned mar will have the following fields specified:
%
% .P     [Nf x d x d] Power Spectral Density matrix
% .C     [Nf x d x d] Coherence matrix
% .dtf   [Nf x d x d] Directed Transfer Function matrix
% .gew   [Nf x d x d] Geweke's frequency domain Granger causality
% .pdc   [Nf x d x d] Partial Directed Coherence
% .L     [Nf x d x d] Phase matrix
% .f     [Nf x 1] Frequency vector
% .ns    Sample rate
%
% dtf(f,i,j) is the DTF at frequency f from signal j to signal i
% pev(f,i,j) is the proportion of power in signal i at frequency f that can
%            be predicted by signal j. 
% gew(f,i,j) is the Granger casuality from signal j to signal i at frequency f.
%            gew=-log(1-pev)
%___________________________________________________________________________
% Copyright (C) 2007 Wellcome Department of Imaging Neuroscience

% Will Penny 
% $Id$

if nargin < 4  | isempty(show)
    show=0;
end

p=mar.p;
d=mar.d;

Nf=length(freqs);
mar.f=freqs;
w=2*pi*freqs/ns;

prec=diag(1./diag(mar.noise_cov));
% Get Power Spectral Density matrix and DTF
for ff=1:Nf,
  af_tmp=eye(d);
  for k=1:p,
    af_tmp=af_tmp+mar.lag(k).a*exp(-i*w(ff)*k);
  end
  iaf_tmp=inv(af_tmp);
  H(ff,:,:)=iaf_tmp;  % transfer function
  mar.P(ff,:,:) = iaf_tmp * mar.noise_cov * iaf_tmp';
  
  % Get DTF and PDC
  for ii=1:d,
      prec_ii=1/sqrt(mar.noise_cov(ii,ii));
      for j=1:d,
          
          % DTF uses iaf_tmp, PDC uses af_tmp
          % DTF normalises wrt rows (to), PDC wrt columns (from)
          mar.dtf(ff,ii,j)=iaf_tmp(ii,j)/sqrt(iaf_tmp(ii,:)*iaf_tmp(ii,:)');
          
          a_j=af_tmp(:,j);
          mar.pdc(ff,ii,j) = real(prec_ii*af_tmp(ii,j)/sqrt(a_j'*prec*a_j));
      end
  end
end

% Get Coherence and Phase
for k=1:d,
    for j=1:d,
        rkj=mar.P(:,k,j)./(sqrt(mar.P(:,k,k)).*sqrt(mar.P(:,j,j)));
        mar.C(:,k,j)=abs(rkj);
        l=atan(imag(rkj)./real(rkj));
        mar.L(:,k,j)=atan(imag(rkj)./real(rkj));
    end
end

% Get Geweke's formulation of Granger Causality in the Frequency domain
C=mar.noise_cov;
for j=1:d,
    for k=1:d,
        rjk=C(j,j)-(C(j,k)^2)/C(k,k);
        sk=abs(mar.P(:,k,k));
        hkj=abs(H(:,k,j)).^2;
        mar.pve(:,j,k)=rjk*hkj./sk;
        mar.gew(:,j,k)=-log(1-mar.pve(:,j,k));
    end
end

if show
    % Plot spectral estimates 
    h=figure;
    set(h,'name','Log Power Spectral Density');
    for k=1:d,
        for j=1:d,
            index=(k-1)*d+j;
            subplot(d,d,index);
            psd=real(mar.P(:,k,j)).^2;
            plot(mar.f,log(psd));
        end
    end
    
    h=figure;
    set(h,'name','Coherence');
    for k=1:d,
        for j=1:d,
            if ~(k==j)
                index=(k-1)*d+j;
                subplot(d,d,index);
                coh=real(mar.C(:,k,j)).^2;
                plot(mar.f,coh);
            end
        end
    end
    
    h=figure;
    set(h,'name','DTF');
    for k=1:d,
        for j=1:d,
            if ~(k==j)
                index=(k-1)*d+j;
                subplot(d,d,index);
                dtf=real(mar.dtf(:,k,j)).^2;
                plot(mar.f,dtf);
            end
        end
    end
    
    h=figure;
    set(h,'name','PDC');
    for k=1:d,
        for j=1:d,
            if ~(k==j)
                index=(k-1)*d+j;
                subplot(d,d,index);
                pdc=mar.pdc(:,k,j);
                plot(mar.f,pdc);
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