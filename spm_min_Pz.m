function x = spm_min_Pz(z, p,s,S)
% used to find the value of z which gives probability p
% used by spm_Pkn
% FORMAT
% z	: hight
% p	: probability to match
% s 	: smoothness
% S	: volume of the search
%_______________________________________________________________________
% @(#)spm_min_Pz.m	1.1 Jean-Baptiste Poline 96/08/22

x = (spm_P(1,s,z,0,S) - p).^2;


