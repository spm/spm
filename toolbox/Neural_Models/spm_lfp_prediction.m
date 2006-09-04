function [G,w] = spm_lfp_prediction(P,M,varargin)
% prediction for log-spectral density of a NMM
% FORMAT [G,w] = spm_lfp_prediction(P,M)
%
% P - parameters
% M - neural mass model stucture
%
% G - ln(G(w))
%__________________________________________________________________________


% compute log-spectral density
%==========================================================================
G  = log(spm_lfp_mtf(P,M));



