function mci_setpath(spm_path, sundials_path)
% Set MATLAB path for MONTE CARLO INFERENCE (MCI) toolbox
% FORMAT mci_setpath(spm_path, sundials_path)
%
% SUNDIALS (SUite of Nonlinear and DIfferential/ALgebraic equation Solvers)
% is available from:
% http://computation.llnl.gov/casc/sundials/main.html
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_setpath.m 6277 2014-12-04 12:16:52Z guillaume $


if isdeployed, return; end

% Set MATLAB path for MCI toolbox
disp('Adding MCI toolbox to MATLAB path');
disp(' ');

mci_path = fileparts(mfilename('fullpath'));
addpath(fullfile(mci_path,'demos'));
addpath(fullfile(mci_path,'gradients'));
addpath(fullfile(mci_path,'inference'));
addpath(fullfile(mci_path,'models'));
addpath(fullfile(mci_path,'plotting'));

% Set MATLAB path for SPM12
disp('Adding SPM12 to MATLAB path');
disp(' ');

if nargin < 1 || isempty(spm_path)
    spm_path = fullfile(mci_path,'..','..');
end
addpath(spm_path);
addpath(fullfile(spm_path,'toolbox','dcm_meeg'));
addpath(fullfile(spm_path,'toolbox','spectral'));

% Set MATLAB path for Sundials
disp('Adding Sundials to MATLAB path');
disp(' ');

if nargin < 2
    sundials_path = fullfile(spm_path,'external','sundials-2.5.0','sundialsTB');
end
addpath(fullfile(sundials_path,'cvodes'));
addpath(fullfile(sundials_path,'cvodes','cvm'));
