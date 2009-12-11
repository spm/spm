function [pE,gE,pC,gC] = spm_phase_priors(DCM,fb,dipfit,freq_prior)
% Prior moments of DCM for phase coupling
% FORMAT [pE,gE,pC,gC] = spm_phase_priors(DCM,fb,dipfit,freq_prior)
%
% freq_prior   Priors on frequency: 'hard_freq' (default),'soft_freq'
% 
% Fields of DCM:
%
% As,Bs{m},Ac,Bc{m} - binary constraints (first two mandatory)
% dipfit       - prior forward model structure
%
% pE - prior expectation - f(x,u,P,M)
% gE - prior expectation - g(x,u,G,M)
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.As    - trial-invariant
%    pE.Bs{m} - trial-dependent
%    pE.Ac    - trial-invariant
%    pE.Bc{m} - trial-dependent
%
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging
 
% Will Penny
% $Id: spm_phase_priors.m 3637 2009-12-11 16:40:15Z will $
 
if nargin < 4 | isempty(freq_prior)
    freq_prior='hard_freq';
end

As=DCM.As;
Bs=DCM.Bs;

if isfield(DCM,'Ac')
    Ac=DCM.Ac;
end
if isfield(DCM,'Bc')
    Bc=DCM.Bc;
end

% orders
%--------------------------------------------------------------------------
n    = size(As,1);                                 % number of sources 

% parameters for electromagnetic forward model
%==========================================================================
G.L  = sparse(1,n);  
U.L  = ones(1,n);

p_exceed=0.001;
z_exceed=abs(spm_invNcdf(0.5*p_exceed,0,1));

% Prior variance of connectivity parameters
prior_sd=fb/z_exceed;
prior_var=prior_sd^2;

% intrinsic connectivity 
%--------------------------------------------------------------------------
E.As  = zeros(size(As));
V.As  = prior_var*As;
 
if isfield(DCM,'Ac')
    E.Ac  = zeros(size(Ac));
    V.Ac  = prior_var*Ac;
end

% Frequency priors
E.df=zeros(n,1);
switch freq_prior
    case 'hard_freq',
        V.df=1e-6*ones(n,1);
    case 'soft_freq',
        new_df=(DCM.options.Fdcm(2)-DCM.options.Fdcm(1))/2;
        new_sig=new_df/3;
        V.df=new_sig*ones(n,1);
    otherwise
        disp('Unknown prior on frequencies');
        return
end

% input-dependent
%--------------------------------------------------------------------------
for i = 1:length(Bs)
    E.Bs{i} = zeros(size(Bs{i}));
    V.Bs{i} = prior_var*Bs{i} & V.As;

end
if isfield(DCM,'Bc')
    for i = 1:length(Bc)
        E.Bc{i} = zeros(size(Bc{i}));
        V.Bc{i} = prior_var*Bc{i} & V.Ac;
    end
end

% prior moments
%--------------------------------------------------------------------------
pE   = E;
gE   = G;
pC   = diag(sparse(spm_vec(V)));
gC   = diag(sparse(spm_vec(U)));
