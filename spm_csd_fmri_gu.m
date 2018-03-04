function [Gu,Gn,w, dt] = spm_csd_fmri_gu(P,dt)
% spectra of neuronal fluctuations and noise
% FORMAT [Gu,Gn,Hz,dt] = spm_csd_fmri_gu(P,dt)
%
% P  - model parameters
% dt - sampling interval
%
% This routine returns the spectra of neuronal fluctuations and noise for a
% standard frequency range specified by the sampling interval
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_fmri_gu.m 7270 2018-03-04 13:08:10Z karl $


% compute log-spectral density
%==========================================================================

% frequencies of interest
%--------------------------------------------------------------------------
Hz1  = 1/128;
Hz2  = 1/(2*dt);
nw   = 32;
w    = linspace(Hz1,Hz2,nw)';

% number of nodes and endogenous (neuronal) fluctuations
%--------------------------------------------------------------------------
nn   = size(P.A,1);
nu   = nn;
form = '1/f';


% spectrum of neuronal fluctuations (Gu) and observation noise (Gn)
%==========================================================================

% experimental inputs
%--------------------------------------------------------------------------
Gu    = zeros(nw,nu,nu);
Gn    = zeros(nw,nn,nn);

% neuronal fluctuations (Gu) (1/f or AR(1) form)
%--------------------------------------------------------------------------
for i = 1:nu
    if strcmp(form,'1/f')
        G     = w.^(-exp(P.a(2,1)));
    else
        G     = spm_mar2csd(exp(P.a(2,1)),w);
    end
    Gu(:,i,i) = Gu(:,i,i) + exp(P.a(1,1))*G/sum(G);
end

% region specific observation noise (1/f or AR(1) form)
%--------------------------------------------------------------------------
for i = 1:nn
    if strcmp(form,'1/f')
        G     = w.^(-exp(P.b(2,1))/2);
    else
        G     = spm_mar2csd(exp(P.b(2,1))/2,w);
    end
    Gn(:,i,i) = Gn(:,i,i) + exp(P.c(1,i))*G/sum(G);
end


% global components
%--------------------------------------------------------------------------
if strcmp(form,'1/f')
    G = w.^(-exp(P.b(2,1))/2);
else
    G = spm_mar2csd(exp(P.b(2,1))/2,w);
end
for i = 1:nn
    for j = i:nn
        Gn(:,i,j) = Gn(:,i,j) + exp(P.b(1,1))*G/sum(G);
        Gn(:,j,i) = Gn(:,i,j);
    end
end


