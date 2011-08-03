function [p] = spm_dirichlet(x,alpha);
% Dirichlet distribution - deprecated
% 
% FORMAT [p] = dirichlet(x,alpha)
% 
% x     - vector of outcome/event probabilities
% alpha - vector of observed events
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% 
% Will Penny
% $Id: spm_dirichlet.m 4418 2011-08-03 12:00:13Z guillaume $

persistent runonce
if isempty(runonce)
    warning('spm_dirichlet is deprecated. Use spm_Dpdf instead.');
    runonce = 1;
end

p = spm_Dpdf(x,alpha);
