function [F,Ep,Cp,pE,pC,Eh] = spm_immune(Y,U,pE,pC,hC)
% Variational inversion of immune model
% FORMAT [F,Ep,Cp,pE,pC,Eh] = spm_immune(Y,U,pE,pC,hC)
% Y   - timeseries data
% pE  - prior expectation of parameters
% pC  - prior covariances of parameters
% hC  - prior covariances of precisions
% 
% F   - log evidence (negative variational free energy)
% Ep  - posterior expectation of parameters
% Cp  - posterior covariances of parameters
% pE  - prior expectation of parameters
% pC  - prior covariances of parameters
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging
 
% Thomas Parr
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% Gaussian priors over model parameters
%--------------------------------------------------------------------------
if nargin < 4
    [pE,pC] = spm_immune_priors;
end
if nargin < 5
    hC = 1/512;
end

% complete model specification
%--------------------------------------------------------------------------
M.G   = @spm_immune_gen;          % generative function
M.FS  = @spm_immune_FS;           % feature selection (link function)
M.pE  = pE;                       % prior expectations (parameters)
M.pC  = pC;                       % prior covariances  (parameters)
M.hE  = 0;                        % prior expectation  (log-precision)
M.hC  = hC;                       % prior covariances  (log-precision)
M.T   = max(U);                   % Simulate up to last measurement
U     = U/24;                     % Convert from hours to days
                       

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);


function Y = spm_immune_FS(Y)
% feature selection for COVID model
% FORMAT Y = spm_COVID_FS(Y)
% P    - model parameters
% M    - model structurel
% U    - inputs or control variables
% Y    - outputs or response variables
%
% This auxiliary function takes appropriate gradients and performs a square
% root transform
%__________________________________________________________________________

% ensure all counts are greater than zero
%--------------------------------------------------------------------------
i      = logical(Y < 1); Y(i) = 1;

% square root transform
%--------------------------------------------------------------------------
Y      = sqrt(Y);
