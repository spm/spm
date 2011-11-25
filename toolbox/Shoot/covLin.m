function [K1,lambda] = covLin(lambda0,settings,args,lab)
% Covariance function for linear regression/classification
% FORMAT [K1,lambda] = covLin(lambda0,settings,args,lab)
% No usage documentation yet
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: covLin.m 4573 2011-11-25 23:01:01Z john $

if ischar(lambda0) && strcmpi(lambda0,'init'),
    K            = settings;
    nK           = norm(K);
    settings     = args;
    if nargin<4
        settings.lab = true(size(K,1),1);
    else
        settings.lab = lab;
    end
    settings.K   = K;
    lambda       = zeros(2,1);
    K1           = settings;
    return;
end

K      = settings.K;
lab    = settings.lab;

lambda = zeros(2,1);%lambda(2)=-12;
lambda(1:numel(lambda0)) = lambda0(:);
K1     = min(exp(lambda(1)),1e12)*K(lab,lab) + exp(lambda(2));

