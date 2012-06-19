function [MVB] = spm_mvb_ui(xSPM,SPM,MVB)
% multivariate Bayes (Bayesian decoding of a contrast)
% FORMAT [MVB] = spm_mvb_ui(xSPM,SPM,MVB)
%
% Sets up, evaluates and saves an MVB structure:
%
% MVB.contrast            % contrast structure
% MVB.name                % name
% MVB.c                   % contrast weight vector
% MVB.M                   % MVB model (see below)
% MVB.X                   % subspace of design matrix
% MVB.Y                   % multivariate response
% MVB.X0                  % null space of design
% MVB.XYZ                 % location of voxels (mm)
% MVB.V                   % serial correlation in response
% MVB.K                   % whitening matrix
% MVB.VOX                 % voxel scaling
% MVB.xyzmm               % centre of VOI (mm)
% MVB.Space               % VOI definition
% MVB.Sp_info             % parameters of VOI
% MVB.Ni                  % number of greedy search steps
% MVB.sg                  % size of reedy search split
% MVB.priors              % model (spatial prior)
% MVB.fSPM                % SPM analysis (.mat file)
%
% where MVB.M contains the following fields:
%
%                F: log-evidence [F(0), F(1),...]
%                G: covariance partition indices
%                h: covariance hyperparameters
%                U: ordered patterns
%               qE: conditional expectation of voxel weights
%               qC: conditional variance of voxel weights
%               Cp: prior covariance (ordered  pattern space)
%               cp: prior covariance (original pattern space)
%
%--------------------------------------------------------------------------
% This routine uses a multivariate Bayesian (MVB) scheme to decode or
% recognise brain states from neuroimages. It resolves the ill-posed
% many-to-one mapping, from voxel values or data features to a target
% variable, using a parametric empirical or hierarchical Bayesian model.
% This model is inverted using standard variational techniques, in this
% case expectation maximisation, to furnish the model evidence and the
% conditional density of the model's parameters. This allows one to compare
% different models or hypotheses about the mapping from functional or
% structural anatomy to perceptual and behavioural consequences (or their
% deficits). The aim of MVB is not to predict (because the outcomes are
% known) but to enable inference on different models of structure-function
% mappings; such as distributed and sparse representations. This allows one
% to optimise the model itself and produce predictions that outperform
% standard pattern classification approaches, like support vector machines.
% Technically, the model inversion and inference uses the same empirical
% Bayesian procedures developed for ill-posed inverse problems (e.g.,
% source reconstruction in EEG).
%
% CAUTION: MVB should not be used to establish a significant mapping
% between brain states and some classification or contrast vector. Its use
% is limited to comparison of different models under the assumption
% (hyperprior) that this mapping exists. To ensure the mapping exists, use
% CVA or compute the randomisation p-value (see spm_mvb_p)
%
% See: spm_mvb and
%
% Bayesian decoding of brain images.
% Friston K, Chu C, Mourao-Miranda J, Hulme O, Rees G, Penny W, Ashburner J.
% Neuroimage. 2008 Jan 1;39(1):181-205
% 
% Multiple sparse priors for the M/EEG inverse problem.
% Friston K, Harrison L, Daunizeau J, Kiebel S, Phillips C, Trujillo-Barreto 
% N, Henson R, Flandin G, Mattout J.
% Neuroimage. 2008 Feb 1;39(3):1104-20.
% 
% Characterizing dynamic brain responses with fMRI: a multivariate approach.
% Friston KJ, Frith CD, Frackowiak RS, Turner R.
% Neuroimage. 1995 Jun;2(2):166-72.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb_ui.m 4770 2012-06-19 13:24:40Z guillaume $
 
 
%-Get figure handles and set title
%--------------------------------------------------------------------------
Finter     = spm_figure('FindWin','Interactive');
spm_results_ui('Clear');
spm_input('!DeleteInputObj');
header     = get(Finter,'Name');
Fmvb       = spm_figure('GetWin','MVB');
spm_clf(Fmvb);
 
%-Get contrast: only the first line of F-contrast 
%--------------------------------------------------------------------------
try
    contrast = SPM.xCon(xSPM.Ic).name;
    c        = SPM.xCon(xSPM.Ic).c(:,1);
catch
    contrast = MVB.contrast;
    c        = MVB.c;
end
 
%-Get VOI name
%--------------------------------------------------------------------------
try
    name = MVB.name;
catch
    name = spm_input('name','-8','s',contrast);
end
name     = strrep(name,' ','_');
name     = ['MVB_' name];
 
%-Get current location {mm}
%--------------------------------------------------------------------------
try
    xyzmm = MVB.xyzmm;
catch
    xyzmm = spm_results_ui('GetCoords');
end

%-Specify search volume
%--------------------------------------------------------------------------
try
    xY   = MVB.xY;
    MVB  = rmfield(MVB,'xY');
catch
    xY   = [];
end
xY.xyz   = xyzmm;
Q        = ones(1,size(SPM.xVol.XYZ,2));
XYZmm    = SPM.xVol.M(1:3,:)*[SPM.xVol.XYZ; Q];

[xY,XYZ,j] = spm_ROI(xY, XYZmm);
 
% Get explanatory variables (data)
%--------------------------------------------------------------------------
XYZ      = XYZmm(:,j);
Y        = spm_get_data(SPM.xY.VY,SPM.xVol.XYZ(:,j));
 
% Check there are intracranial voxels
%--------------------------------------------------------------------------
if isempty(Y)
    spm('alert*',{'No voxels in this VOI';'Please use a larger volume'},...
        'Multivariate Bayes');
end
 
%-Get model[s]
%--------------------------------------------------------------------------
try
    priors = lower(MVB.priors);
catch
    str    = {'compact','sparse','smooth','support'};
    Ip     = spm_input('model (spatial prior)','!+1','m',str);
    priors = str{Ip};
end

%-Number of iterations
%--------------------------------------------------------------------------
try
    sg     = MVB.sg;
catch
    str    = 'size of successive subdivisions';
    sg     = spm_input(str,'!+1','e',.5);
end
 
% MVB is now specified
%==========================================================================
spm('Pointer','Watch')
 
%-Get target and confounds
%--------------------------------------------------------------------------
X   = SPM.xX.X;
X0  = X*(speye(length(c)) - c*pinv(c));
try
    % accounting for multiple sessions
    %----------------------------------------------------------------------
    tmpX0  = [];
    for ii = 1:length(SPM.xX.K)
        tmp   = zeros(sum(SPM.nscan),size(SPM.xX.K(ii).X0,2));
        tmp(SPM.xX.K(ii).row,:) = SPM.xX.K(ii).X0;
        tmpX0 = [tmpX0 tmp];
    end
    X0 = [X0 tmpX0];
end
X   = X*c;
 
% serial correlations
%--------------------------------------------------------------------------
V   = SPM.xVi.V;
 
% invert
%==========================================================================
VOX      = diag(abs(SPM.xVol.M));
U        = spm_mvb_U(Y,priors,X0,XYZ,VOX);
try
    Ni   = MVB.Ni;
catch
    str  = 'Greedy search steps';
    Ni   = spm_input(str,'!+1','i',max(8,ceil(log(size(U,2))/log(1/sg))));
end
M        = spm_mvb(X,Y,X0,U,V,Ni,sg);
M.priors = priors;
 
% assemble results
%--------------------------------------------------------------------------
MVB.contrast = contrast;                    % contrast of interest
MVB.name     = name;                        % name
MVB.c        = c;                           % contrast weight vector
MVB.M        = M;                           % MVB model (see below)
MVB.X        = X;                           % subspace of design matrix
MVB.Y        = Y;                           % multivariate response
MVB.X0       = X0;                          % null space of design
MVB.XYZ      = XYZ;                         % location of voxels (mm)
MVB.V        = V;                           % serial correlation in repeosne
MVB.K        = full(V)^(-1/2);              % whitening matrix
MVB.VOX      = SPM.xVol.M;                  % voxel scaling
MVB.xyzmm    = xyzmm;                       % centre of VOI (mm)
MVB.Space    = xY.def;                      % VOI definition
MVB.Sp_info  = xY.spec;                     % parameters of VOI
MVB.Ni       = Ni;                          % number of greedy search steps
MVB.sg       = sg;                          % size of reedy search split
MVB.priors   = priors;                      % model (spatial prior)
MVB.fSPM     = fullfile(SPM.swd,'SPM.mat'); % SPM analysis (.mat file)
 
% display
%==========================================================================
spm_mvb_display(MVB)
 
% save
%--------------------------------------------------------------------------
save(fullfile(SPM.swd,[name '.mat']),'MVB', spm_get_defaults('mat.format'))
assignin('base','MVB',MVB)

%-Reset title
%--------------------------------------------------------------------------
set(Finter,'Name',header)
spm('Pointer','Arrow')
