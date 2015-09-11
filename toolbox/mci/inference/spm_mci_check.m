function [corr] = spm_mci_check (M)
% Check model structure M is correctly specified
% FORMAT [corr] = spm_mci_check (M)
%
% corr      1 for correctly specified model
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_check.m 6548 2015-09-11 12:39:47Z will $

corr=1;
if ~isfield(M,'N')
    disp('Number of time points M.N not specified');
    corr=0;
end
if ~isfield(M,'T')
    disp('Total time M.T not specified');
    corr=0;
end
if ~isfield(M,'t')
    disp('Time point vector M.t not specified');
    corr=0;
end
if ~isfield(M,'m')
    disp('Number of inputs M.m not specified');
    corr=0;
end
if ~isfield(M,'n')
    disp('Number of states M.n not specified');
    corr=0;
end
if ~isfield(M,'l')
    disp('Number of outputs M.l not specified');
    corr=0;
end
if ~isfield(M,'x0')
    disp('Initial state M.x0 not specified');
    corr=0;
end
if ~isfield(M,'x')
    disp('Reshaped state template M.x not specified');
    corr=0;
end
if ~isfield(M,'f')
    disp('Flow function M.f not specified');
    corr=0;
end
if ~isfield(M,'g')
    disp('Observation function M.g not specified');
    corr=0;
end
if ~isfield(M,'int')
    disp('Integration method M.int not specified');
    corr=0;
end
if ~isfield(M,'Np')
    disp('Number of parameters M.Np not specified');
    corr=0;
end
if ~isfield(M,'pE')
    disp('Prior mean M.pE not specified');
    corr=0;
end
if ~isfield(M,'pC')
    disp('Prior covariance M.pC not specified');
    corr=0;
end
if ~isfield(M,'ipC')
    disp('Inverse prior covariance M.ipC not specified');
    corr=0;
end
if ~isfield(M,'L')
    disp('Likelihood function M.L not specified');
    corr=0;
end
if ~isfield(M,'Ce')
    disp('Observation noise covariance M.Ce not specified');
    corr=0;
end

if corr
    disp('Model is correctly specified');
end