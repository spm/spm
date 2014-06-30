function [voxel] = spm_vb_get_Ab(Y,slice)
% Get A and b quantities - average prediction errors from AR model
% FORMAT [voxel] = spm_vb_get_Ab(Y,slice)
% 
% Y      - [T x N] time series
% slice  - data structure (see spm_vb_glmar)
% 
% voxel(n).A  
% voxel(n).b
%
% The above quantities are estimated using pre-computed
% cross-covariance matrices
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_get_Ab.m 6079 2014-06-30 18:25:37Z spm $

k = slice.k;
N = slice.N;
    
for n=1:N
    % Equation 63 of paper VB1 but implemented 
    % efficiently using cross-covariance method described in paper VB3
    if isfield(slice.I,'A2_tilde')
        A2_tilde  = slice.I.A2_tilde(:,:,n);
    else
        A2_tilde  = reshape(slice.I.S*slice.a2{n}(:),k,k);
    end
    if isfield(slice.I,'A3a_tilde')
        A3a_tilde = slice.I.A3a_tilde(:,:,n);
    else
        A3a_tilde = -reshape(slice.I.R1*slice.ap_mean(:,n),k,k);
    end
    voxel(n).A    = slice.I.xtx+A2_tilde+A3a_tilde+A3a_tilde';
    
    % Equation 64 of paper VB1 but implemented 
    % efficiently using cross-covariance method described in paper VB3
    b2_tilde   = -slice.I.rxy(:,:,n)'*slice.ap_mean(:,n);
    b3_tilde   = -slice.I.Gxy(:,:,n)'*slice.ap_mean(:,n);
    b4_tilde   = slice.I.D(:,:,n)*slice.a2{n}(:);
    voxel(n).b = slice.I.gxy(:,n)+b2_tilde+b3_tilde+b4_tilde;
end
