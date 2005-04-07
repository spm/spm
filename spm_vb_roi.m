function [F,pm,Fabs] = spm_vb_roi (VOI_fname,SPM,model)
% Compare sp_glm_ar models for a cluster of interest
% FORMAT [F,pm,Fabs] = spm_vb_roi (VOI_fname,SPM,model)
%
% VOI_fname     VOI filename
% SPM           SPM data structure
% model         data structure specifying which basis functions to compare
%
% F             model evidence
% pm            posterior model probability
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny 
% $Id$

expr=['load ',VOI_fname];
eval(expr);

if nargin == 1
    %-Get SPM.mat 
    swd     = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat'),'H');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd = swd;
end

Y=xY.y;
xyz=xY.XYZmm;
N=size(xyz,2);
M=diag(SPM.xVol.M);
m=1./abs(M(1:3));
xyz=(m*ones(1,N)).*xyz;
vxyz = spm_vb_neighbors(xyz',1);

% Set number of AR coefficients
try 
    SPM.PPM.AR_P;
catch
    SPM.PPM.AR_P = 3;
end

% Specify type of prior for regression coefficients
try
    SPM.PPM.priors.W;
catch
    if N==1
        SPM.PPM.priors.W = 'Voxel - Shrinkage';
    else
        SPM.PPM.priors.W = 'Spatial - GMRF';
    end
end

% Specify type of prior for AR coefficients
try
    SPM.PPM.priors.A;
catch
    if N==1
        SPM.PPM.priors.A = 'Voxel - Shrinkage';
    else
        SPM.PPM.priors.A = 'Spatial - GMRF';
    end
end
    
% Get matrices that will remove low-frequency drifts 
% if high pass filters have been specified
s=1;
sess_nScan=length(SPM.xX.K(s).row);
if size(SPM.xX.K(s).X0,2) > 0
    X0=SPM.xX.K(s).X0;
    hpf(s).R0=eye(sess_nScan)-X0*pinv(X0);
else
    hpf(s).R0=eye(sess_nScan);
end    
% Filter data to remove low frequencies
R0Y=hpf(s).R0*Y(SPM.Sess(s).row,:);

% Set optimisation parameters
try
    SPM.PPM.maxits;
catch
    SPM.PPM.maxits=16;
end
try
    SPM.PPM.tol;
catch
    SPM.PPM.tol=0.00001;
end

% Specify basis functions
try
    model
catch
    model(1).name='hrf';
    model(1).sname='Inf-1';
    model(1).order=1;
    model(1).length=32; % length in seconds
    
    model(2).name='hrf (with time derivative)';
    model(2).sname='Inf-2';
    model(2).order=2;
    model(2).length=32;
    
    model(3).name='hrf (with time and dispersion derivatives)';
    model(3).sname='Inf-3';
    model(3).order      = 3;  
    model(3).length     = 32;                
    
    model(4).name='Fourier set';
    model(4).sname='F    ';
    model(4).order=5;
    model(4).length=20;
    
    model(5).name='Fourier set (Hanning)';
    model(5).sname='FH   ';
    model(5).order=5;
    model(5).length=20;
    
    model(6).name='Gamma functions';
    model(6).sname='Gamm3';
    model(6).order      = 3;  
    model(6).length     = 32;                
    
    model(7).name='Finite Impulse Response';
    model(7).sname='FIR  ';
    model(7).order      = 10; 
    model(7).length     = 20;                
end

M=length(model);
F=[];

sname=[];
% Fit models
for m=1:M,
    
    % Switch basis functions
    SPM.xBF.name=model(m).name;
    SPM.xBF.order=model(m).order;
    SPM.xBF.length=model(m).length;
    SPM.xBF = spm_get_bf(SPM.xBF);
    
    SPM = spm_fmri_design(SPM,0);  % 0 to override writing of SPM.mat file
    slice = spm_vb_init_volume (SPM.xX.X,SPM.PPM.AR_P);
    
    slice.maxits=SPM.PPM.maxits;
    slice.tol=SPM.PPM.tol;
    slice.compute_det_D=1;
    slice.verbose=1;
    slice.update_w=1;
    slice.update_lambda=1;
    slice.update_F=1;
    
    slice = spm_vb_set_priors(slice,SPM.PPM.priors,vxyz);
    slice = spm_vb_glmar(R0Y,slice);
    F=[F,slice.F];
    sname=[sname; model(m).sname];
end
Fabs=F;

F=F-mean(F);
figure
pm=exp(F)./sum(exp(F));
if sum(isnan(pm))
    [pF,pi]=max(F);
    pm=zeros(M,1);
    pm(pi)=1;
end
bar(pm);
set(gca,'XTickLabel',sname);
ylabel('Posterior model probability');
set(gca,'FontSize',12);

figure
F=F-min(F);
bar(F);
set(gca,'XTickLabel',sname);
ylabel('Log Evidence');
set(gca,'FontSize',12);

