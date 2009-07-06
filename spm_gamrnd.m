function r = spm_gamrnd(a,b,varargin)
% Random arrays from gamma distribution - a compiled routine
% FORMAT r = spm_gamrnd(a,b,m,n,...)
%
% a        - shape parameter
% b        - scale parameter
% m,n,...  - dimensions of the output array [optional]
%
% r        - array of random numbers chosen from the gamma distribution
%__________________________________________________________________________
%
% Reference
% 
% George Marsaglia and Wai Wan Tsang, "A Simple Method for Generating Gamma
% Variables": ACM Transactions on Mathematical Software, Vol. 26, No. 3,
% September 2000, Pages 363-372
% http://portal.acm.org/citation.cfm?id=358414
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_gamrnd.m 3251 2009-07-06 17:29:44Z guillaume $

%-This is merely the help file for the compiled routine
error('spm_gamrnd.c not compiled - see Makefile');
