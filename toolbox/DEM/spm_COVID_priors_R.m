function [pE,pC,str,erc] = spm_COVID_priors_R(data)
% Prior expectation and covariance of between region parameters
% FORMAT [pE,pC,str,erc] = spm_COVID_priors_R(data)
% 
% data(N)     - Meta data, including distance between regions
% 
% pE          - prior expectation (structure)
% pC          - prior covariances (structure)
% str.names   - parameter names
% str.regions - regional names
%
% This routine assembles the (Gaussian) and priors over the parameters of a
% generative model for COVID-19. This generative model is based upon a mean
% field approximation to ensemble of population dynamics, in which four
% marginal distributions are coupled through probability transition
% matrices. The marginal distributions correspond to 4 factors;
% namely,location, infection, clinical and diagnostic or testing states.
% Please see spm_COVID_priors for details. This routine prepares the priors
% for the parameters that couple different regions (e.g., American states).
% These parameters include the (effective) connectivity that controls the
% flux of people from one region to another. The total population size in
% these models is included as a precise prior, while the number of initial
% cases is allowed to vary to accommodate differential onset times.
% Finally, there is a federal parameter that determines the balance between
% region specific and federal densities in mediating lockdown or social
% distancing strategies.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% set up
%--------------------------------------------------------------------------
if nargin < 1
    data.state = 'region';
    data.dis   = 0; 
end

% assemble priors
%==========================================================================
N      = numel(data);                   % number of regions
erc    = 1 - eye(N,N);                  % full connectivity (all)
ln     = @(x)log(x + exp(-32));

% supplement with (priors over) between-region parameters
%--------------------------------------------------------------------------
pE.n   = zeros(N,1);                    % initial (log) infection count
pE.erc = ln(erc) - 8;                   % effective (regional) connectivity
pE.fed = ln(0);                         % social distancing (federal)

% prior variances
%--------------------------------------------------------------------------
pC.n   = ones(N,1);                     % uninformative priors on onset
pC.erc = erc;                           % effective (regional) connectivity
pC.fed = 0;                             % social distancing (federal)

% augment strings with parameter field names
%--------------------------------------------------------------------------
pnam  = fieldnames(pE);
name  = {};
for i = 1:numel(pnam)
    for j = 1:numel(pE.(pnam{i}))
        name{end + 1} = pnam{i};
    end
end
str.names   = name;
str.regions = {data.state};


return

% notes for different models of sparse coupling between regions
%==========================================================================
% ons  = [data.ons]; ons = ons(1:N,1:N);
% dis  = [data.dis]; dis = dis(1:N,1:N);
% erc  = ons < 64;                      % timing  distance
% erc  = dis < 64;                      % spatial distance
% erc  = ones(N,N) - eye(N,N);          % full connectivity (all)

% [d,i] = sort([data.lat], 'descend');  % latitude
% [d,i] = sort([data.first],'ascend');  % date of first case
% [d,i] = sort([data.cum], 'descend');  % cumulative deaths
% [d,i] = sort([data.lon], 'descend');  % longitude
% 
% erc   = full(spm_speye(N,N,1) + spm_speye(N,N,-1));
% erc   = erc(i,i);

% effective regional connectivity assessed through Bayesian model
% comparison
%--------------------------------------------------------------------------
% all DCM{1}.F -3.3017e+03
% all DCM{1}.F -3.4440e+03
% cumulative deaths DCM{1}.F -2.8567e+03
% cumulative deaths DCM{2}.F -2.9955e+03
% date of first DCM{1}.F -2.9640e+03
% date of first DCM{2}.F -2.9540e+03
% latitude DCM{1}.F -2.8883e+03
% latitude DCM{2}.F -2.9777e+03
% longitude DCM{1}.F -2.6833e+03
% longitude DCM{2}.F -2.8286e+03
% timing distance DCM{1}.F -3.0856e+03
% timing distance DCM{2}.F  -3.2311e+03
% spatial distance DCM{1}.F -3.0549e+03
% spatial distance DCM{2}.F -3.4772e+03
%--------------------------------------------------------------------------



