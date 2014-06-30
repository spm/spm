function [F] = spm_vb_regionF (Y,xY,SPM)
% Get log model evidence over a region of data for a GLM
% FORMAT [F] = spm_vb_regionF (Y,xY,SPM)
%
% Y     Matrix of fMRI data (eg. from spm_summarise.m)
% xY    Coordinates etc from region (eg. from spm_voi.m)
% SPM   SPM data structure (this must be loaded in from an 
%       SPM.mat file). If this field is not specified this function
%       will prompt you for the name of an SPM.mat file
%
% F     Log model evidence (single number for whole region)
%
% Importantly, the design matrix is normalised so that when you compare
% models their regressors will be identically scaled.
%
% Valid model comparisons also require that the DCT basis set used in high 
% pass filtering, as specified in SPM.xX.K.X0, is the same for all models 
% that are to be compared.
% 
% W. Penny, G. Flandin, and N. Trujillo-Barreto. (2007). Bayesian Model 
% Comparison of Spatially Regularised General Linear Models. Human 
% Brain Mapping, 28(4):275-293.
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_regionF.m 6079 2014-06-30 18:25:37Z spm $

try 
    SPM;
catch
    [Pf, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, return; end
    swd = spm_file(Pf,'fpath');
    load(fullfile(swd,'SPM.mat'))
end

s=length(SPM.Sess);
if s>1
    disp('spm_vb_regionF: only works for single session data');
    return
end

% Normalise design matrices to have columns of unit
% variance - except for last
X=SPM.xX.X;
X=norm_design(X);
   
% Make residual forming matrix (to later remove drifts)
X0=SPM.xX.K.X0;
N=size(X0,1);
iX0=pinv(X0);
R=eye(N)-X0*iX0;

% Remove drifts
Y=R*Y;

% Normalise Y to unit dev 
sigma=std(Y);
N=size(Y,1);
Y=Y./(ones(N,1)*sigma);

xyz  = xY.XYZmm;
N    = size(xyz,2);
xyz  = SPM.xVol.M\[xyz;ones(1,size(xyz,2))];
xyz  = xyz(1:3,:);
vxyz = spm_vb_neighbors(xyz',1);

priors.W='Voxel - Shrinkage';
priors.A='Voxel - Shrinkage';

block=spm_vb_init_volume(X,1);
block=spm_vb_set_priors(block,priors,vxyz);
block.maxits = 64;
block.update_alpha  = 1;
block.update_w      = 1;
block.update_lambda = 1;
block.update_F      = 1;
block=spm_vb_glmar(Y,block);

F=block.F;


function [Xout] = norm_design (Xin)
% Normalise design matrix to have unit variance columns
% FORMAT [Xout] = norm_design (Xin)

X=Xin(:,1:end-1);
s=std(X,[],1);
iS=diag(1./s);
Xout=[X*iS,Xin(:,end)];