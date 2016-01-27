function [noise,M] = spm_mci_obsnoise (w,v,assign,noise,M,U,Y)
% Update observation noise
% FORMAT [noise,M] = spm_mci_obsnoise (w,v,assign,noise,M,U,Y)
%
% w         random effects
% v         fixed effects
% assign    for dynamical models this structure specifies whether init
%           states, flow and o/p params are random, fixed or known
% noise     observation noise structure
% M         model structures
% U         input structures
% Y         data structures
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_obsnoise.m 6697 2016-01-27 14:57:28Z spm $

N=length(M);

err=[];
for n=1:N,
    if ~isempty(w), wsub=w(:,n); else, wsub=[]; end
    try, ind=Y{n}.ind; catch, ind=1:M{n}.N; end
    if isfield(M{n},'IS')
        % Other model type
        yhat = feval(M{n}.IS,wsub,M{n},U{n});
        err=[err; Y{n}.y-yhat];
    else
        % Differential equation model
        [Pinit,Pflow] = spm_mci_init_flow (assign,wsub,v,M{n});
        M{n}.x0=Pinit;
        yhat = spm_mci_fwd (Pflow,M{n},U{n});
        err=[err; Y{n}.y-yhat(ind,:)];
    end
end
NT=size(err,1);
noise.cN=noise.c0+0.5*NT;
noise.DN=noise.D0+0.5*NT*cov(err);
L=spm_wishrnd(noise.DN,noise.cN);
Ce=inv(L);
for n=1:N,
    M{n}.Ce=Ce;
end