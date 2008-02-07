function [G,G1,G2,G3] = spm_vb_get_Gn (Y,slice,n)
% Compute Gn for VB-GLM-AR modelling 
% FORMAT [G,G1,G2,G3] = spm_vb_get_Gn (Y,slice,n)
%
% Y             [T x N] time series 
% slice         data structure 
% n             voxel number
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_get_Gn.m 1143 2008-02-07 19:33:33Z spm $

p=slice.p;
k=slice.k;

a_mean_t=slice.ap_mean(:,n);
a_mean=a_mean_t';
block_n = [(n-1)*k+1:n*k];
w_mean=slice.w_mean(block_n,1);
w2=slice.w2{n};
a2=slice.a2{n};

% Equations 77 to 81 of paper VB1 but implemented 
% efficiently using cross-covariance method described in paper VB3
G11=slice.y2(n);
G12=ones(1,p)*(a2.*slice.I.Gy(:,:,n))*ones(p,1);
G13=-2*slice.I.gy(:,n)'*a_mean_t;
G1=G11+G12+G13;

G21=ones(1,k)*(w2.*slice.I.Gx)*ones(k,1);
G22=trace(slice.I.A2_tilde(:,:,n)*slice.w_cov{n});
G23=ones(1,k)*((w_mean*w_mean').*slice.I.A2_tilde(:,:,n))*ones(k,1);
G24=2*ones(1,k)*(w2.*slice.I.A3a_tilde(:,:,n))*ones(k,1);
G2=G21+G22+G23+G24;

G31=-2*w_mean'*slice.I.gxy(:,n);
G32=2*a_mean*slice.I.rxy(:,:,n)*w_mean;
G33=2*a_mean*slice.I.Gxy(:,:,n)*w_mean;
G34=-2*ones(1,p)*(a2.*slice.I.W_tilde(:,:,n))*ones(p,1);
G3=G31+G32+G33+G34;

G=G1+G2+G3;


