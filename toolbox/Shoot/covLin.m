function [K1,lambda] = covLin(lambda0,settings,args,lab)
% Covariance function for linear regression/classification
% FORMAT [K1,lambda] = covLin(lambda0,settings,args,lab)
% No usage documentation yet
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: covLin.m 6881 2016-09-19 09:48:54Z john $

if ischar(lambda0) && strcmpi(lambda0,'init'),
    K            = settings;
    settings     = args;
    if ~iscell(K), K = {K}; end;
    if nargin<4
        settings.lab = true(size(K{1},1),1);
    else
        settings.lab = lab;
    end
    lambda       = zeros(numel(K)+1,1);
    for i=1:numel(K)
        lambda(i) = -log(norm(K{i}));
    end
    lambda(end) = -log(norm(ones(size(K{1}))));

    settings.K   = K;
    K1           = settings;
    return;
end

K      = settings.K;
lab    = settings.lab;

lambda = zeros(numel(K)+1,1);
lambda(1:numel(lambda0)) = lambda0(:);
lambda = max(min(lambda,27),-27);
K1     = exp(lambda(end));
for i=1:numel(K)
    K1 = K1 + exp(lambda(i))*K{i}(lab,lab);
end

%fprintf('\t%f',lambda); fprintf('\n');
