function x = min_Pn(n, p,s,u,S)
% used to find the value of n which gives probability p
% used by spm_Pkn
% FORMAT
% n	: spatial extent of the cluster
% u	: threshold
% p	: probability to match
% s 	: smoothness
% S	: volume of the search
%_______________________________________________________________________
% @(#)spm_min_Pn.m	1.1 Jean-Baptiste Poline 96/08/22

x = (spm_P(1,s,u,n,S) - p).^2;


