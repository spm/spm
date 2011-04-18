function [post,model] = spm_vb_roi_anova (VOI_fname,SPM,factor)
% Bayesian ANOVA for a region of interest
% FORMAT [post,model] = spm_vb_roi_anova (VOI_fname,SPM,factor)
%
% VOI_fname     VOI filename
% SPM           SPM data structure
% factor        data structure relating conditions to levels of factors
%
% model         data structure describing models
%               (m).F             model evidence
%               (m).X             design matrix
% post          Posterior probabilities of
%               .factor1        main effect of factor 1
%               .factor2        main effect of factor 2
%               .interaction    interaction
%               .average        average
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_roi_anova.m 4310 2011-04-18 16:07:35Z guillaume $     

if nargin == 1
    %-Get SPM.mat 
    swd     = spm_str_manip(spm_select(1,'SPM.mat','Select SPM.mat'),'H');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd = swd;
end

expr=['load ',VOI_fname];
eval(expr);

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
    SPM.PPM.AR_P = 0;
end

% Specify type of prior for regression coefficients
try
    SPM.PPM.priors.W;
catch
    if N==1
        SPM.PPM.priors.W = 'Voxel - Shrinkage';
    else
        SPM.PPM.priors.W = 'Spatial - LORETA';
    end
end

% Specify type of prior for AR coefficients
try
    SPM.PPM.priors.A;
catch
    if N==1
        SPM.PPM.priors.A = 'Voxel - Shrinkage';
    else
        SPM.PPM.priors.A = 'Spatial - LORETA';
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
% SPM.xBF.name='hrf';
% SPM.xBF.order=1;
SPM.xBF.name='hrf (with time derivative)';
SPM.xBF.order=2;
SPM.xBF.length=32;
SPM.xBF = spm_get_bf(SPM.xBF);




nf=length(factor);

model = spm_vb_models (SPM,factor);

% Fit models
for m=1:6,
    
    if nf==2 || (nf==1 && (m==1 || m==2 || m==6))
        % fit model
        SPM=original_SPM;
        
        if ~(m==1 || m==6)
            % Get design matrix for relevant input set
            SPM.Sess(1).U=model(m).U;
            SPM.Sess(1).U=spm_get_ons(SPM,1);
            SPM=spm_fmri_design(SPM,0);     % Don't write SPM.mat file
            model(m).X=SPM.xX.X;
        end
        slice = spm_vb_init_volume (model(m).X,SPM.PPM.AR_P);
        
        slice.maxits=SPM.PPM.maxits;
        slice.tol=SPM.PPM.tol;
        slice.compute_det_D=1;
        slice.verbose=1;
        slice.update_w=1;
        slice.update_lambda=1;
        slice.update_F=1;
        slice = spm_vb_set_priors(slice,SPM.PPM.priors,vxyz);
        slice = spm_vb_glmar(R0Y,slice);
        
        model(m).F=slice.F;
        
        model(m).slice=slice;
    end
end

if nf==2
    
    F=[model(3).F,model(2).F];
    F=F-mean(F);
    post.factor1=exp(F(1))/sum(exp(F));
    
    F=[model(4).F,model(2).F];
    F=F-mean(F);
    post.factor2=exp(F(1))/sum(exp(F));
    
    F=[model(6).F,model(5).F];
    F=F-mean(F);
    post.interaction=exp(F(1))/sum(exp(F));
    
    F=[model(2).F,model(1).F];
    F=F-mean(F);
    post.average=exp(F(1))/sum(exp(F));
    
elseif nf==1
    
    F=[model(2).F,model(1).F];
    F=F-mean(F);
    post.average=exp(F(1))/sum(exp(F));
    
    F=[model(6).F,model(2).F];
    F=F-mean(F);
    post.factor1=exp(F(1))/sum(exp(F));
    
end




