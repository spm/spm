function [post,model] = spm_vb_anova (VOI_fname,SPM,factor)
% Bayesian ANOVA for a region of interest
% FORMAT [post,model] = spm_vb_anova (VOI_fname,SPM,factor)
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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny 
% $Id$     

if nargin == 1
    %-Get SPM.mat 
    swd     = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat'),'H');
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

original_SPM=SPM;
NC=length(SPM.Sess(1).U); % Number of conditions

dt=SPM.Sess(1).U(1).dt;
P=SPM.Sess(1).U(1).P;

nf=length(factor);

% Create full model (same as original)
model(6).name='Full';
model(6).U=SPM.Sess(1).U;
model(6).X=SPM.xX.X;

% Create null model
model(1).name='NULL';
null_cols=[SPM.xX.iH,SPM.xX.iB,SPM.xX.iG];
model(1).X=SPM.xX.X(:,null_cols);

% Create average model
model(2).name='Average';
model(2).U(1).name{1}='Average';
model(2).U(1).dt=dt;
ons=[];
dur=[];
for i=1:NC,
    ons=[ons;SPM.Sess(1).U(i).ons];
    dur=[dur;SPM.Sess(1).U(i).dur];
end
model(2).U(1).ons=ons;
model(2).U(1).dur=dur;
model(2).U(1).P=P;
    
if nf==2
    % For two factor models
    
    % Create model for factor A
    model(3).name=factor(1).name;
    nA=length(factor(1).level);
    for i=1:nA,
        % Has as many inputs as levels of factor A
        ons=[];dur=[];
        nAi=length(factor(1).level(i).conditions);
        for j=1:nAi,
            % Concatenate inputs from relevant conditions
            c=factor(1).level(i).conditions(j);
            ons=[ons;original_SPM.Sess(1).U(c).ons];
            dur=[dur;original_SPM.Sess(1).U(c).dur];
        end
        model(3).U(i).name{1}=['Level ',int2str(i)];
        model(3).U(i).dt=dt;
        model(3).U(i).ons=ons;
        model(3).U(i).dur=dur;
        model(3).U(i).P=P;
    end
    
    % Create model for factor B
    model(4).name=factor(2).name;
    nB=length(factor(2).level);
    for i=1:nB,
        % Has as many inputs as levels of factor A
        ons=[];dur=[];
        nBi=length(factor(2).level(i).conditions);
        for j=1:nBi,
            % Concatenate inputs from relevant conditions
            c=factor(2).level(i).conditions(j);
            ons=[ons;original_SPM.Sess(1).U(c).ons];
            dur=[dur;original_SPM.Sess(1).U(c).dur];
        end
        model(4).U(i).name{1}=['Level ',int2str(i)];
        model(4).U(i).dt=dt;
        model(4).U(i).ons=ons;
        model(4).U(i).dur=dur;
        model(4).U(i).P=P;
    end
    
    % Create model for both factors
    model(5).name='Both factors';
    nboth=nA+nB;
    i=1;
    for k=1:nA,
        model(5).U(i)=model(3).U(k);
        i=i+1;
    end
    % Note: to avoid rank deficiency don't estimate 
    % average response to level nB of factor 2 as 
    % this is given by a linear combination of responses in
    % other levels eg. for 2by2: B2= A1+A2-B1
    for k=1:nB-1,
        model(5).U(i)=model(4).U(k);
        i=i+1;
    end
end

% Fit models
for m=1:6,
    
    if nf==2 | (nf==1 & (m==1 | m==2 | m==6))
        % fit model
        SPM=original_SPM;
        
        if ~(m==1 | m==6)
            % Get design matrix for relevant input set
            SPM.Sess(1).U=model(m).U;
            SPM.Sess(1).U=spm_get_ons(SPM,1);
            SPM=spm_fmri_design(SPM);
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

% Put back original SPM
SPM=original_SPM;
save SPM SPM;


