function [slice] = spm_vb_set_priors (slice,priors,vxyz)
% Set up precision parameters and optimisation options to correspond to priors
% FORMAT [slice] = spm_vb_set_priors (slice,priors,vxyz)
%
% slice     VB-GLM-AR data structure
% priors    For regression coefficients:
%
%           .W= 'Spatial - GMRF' : coeffs that are spatially regularised
%               'Spatial - LORETA' : as above but different spatial prior
%               'Voxel - Shrinkage' : that are shrunk voxel-wise
%               'Voxel - Uninformative' : coeffs without prior
%
%           For AR coefficients:
%
%           .A= 'Spatial - GMRF' : coeffs that are spatially regularised
%               'Spatial - LORETA' : as above but different spatial prior
%               'Voxel - Shrinkage' : that are shrunk voxel-wise
%               'Voxel - Uninformative' : coeffs without prior
%               'Voxel - Limiting' :  Voxel-specific coeffs as limiting case of
%                                   a LORETA prior
%               'Slice - Limiting' :   Slice-specific coeffs as limiting case of
%                                   a LORETA prior
%               'Discrete'  :   Different coeffs as function of
%                                   grey/white/CSF or other masks
%
% vxyz      locations of voxels 
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny 
% $Id: spm_vb_set_priors.m 209 2005-08-02 17:09:13Z will $

if ~isfield(slice,'verbose')
    slice.verbose=0;
end

N=size(vxyz,1);
k=size(slice.X,2);
p=slice.p;

slice.update_alpha=1; 
slice.update_beta=1;

% Set c coefficients of gamma priors and default mean
slice.c_alpha_prior=0.1*ones(k,1);
slice.c_beta_prior=0.1*ones(p,1);
slice.c_lambda_prior=0.1*ones(N,1);

slice.mean_alpha=ones(k,1);
slice.mean_beta=1000*ones(p,1);
slice.mean_lambda=ones(N,1);

switch priors.W,
    case 'Spatial - GMRF',
        slice.Dw=spm_vb_spatial_precision (vxyz,'Spatial - GMRF');
        
    case 'Spatial - LORETA',
        slice.Dw=spm_vb_spatial_precision (vxyz,'Spatial - LORETA');
        
    case 'Spatial - LORETAu',
        slice.Dw=spm_vb_spatial_precision (vxyz,'Spatial - LORETAu');
        
    case 'Voxel - Shrinkage',
        slice.Dw=speye(N);
        slice.update_alpha=1;
        
    case 'Voxel - Uninformative',
        slice.Dw=speye(N);
        slice.mean_alpha=0.000001*ones(k,1);
        slice.update_alpha=0;
end

if slice.verbose
    fprintf('\n');
    disp(sprintf('%s priors for W',priors.W));
end

if slice.p > 0
    slice.update_a=1;
end

switch priors.A,
    
    case 'Spatial - GMRF',
        slice.Da=spm_vb_spatial_precision (vxyz,'Spatial - GMRF');
        
    case 'Spatial - LORETA',
        slice.Da=spm_vb_spatial_precision (vxyz,'Spatial - LORETA');
        
    case 'Spatial - LORETAu',
        slice.Da=spm_vb_spatial_precision (vxyz,'Spatial - LORETAu');
        
    case 'Voxel - Shrinkage',
        slice.Da=speye(N);
        slice.update_beta=1;
        
    case 'Voxel - Uninformative',
        slice.Da=speye(N);
        slice.mean_beta=0.000001*ones(p,1);
        slice.update_beta=0;
        
    case 'Voxel - Limiting',
        % Very low (fixed) spatial precision ensures spatial prior is 
        % effectively ignored
        slice.Da=spm_vb_spatial_precision (vxyz,'Spatial - LORETA');
        slice.update_beta=0;
        slice.mean_beta=1e-4*ones(p,1);
        
    case 'Slice - Limiting',
        % Very high (fixed) spatial precision ensures AR coefficients are
        % (nearly) identical over the slice
        slice.Da=spm_vb_spatial_precision (vxyz,'Spatial - LORETA');
        slice.update_beta=0;
        slice.mean_beta=1e4*ones(p,1);
    
    case 'Discrete',
        disp('Different AR coeff priors for masked regions');
        if ~(isfield(priors,'gamma'))
            disp('Must enter class labels');
            return
        end
        priors.S=size(priors.gamma,2);
        priors.N=sum(priors.gamma(:,:));
        for s=1:priors.S
            priors.voxel(s).i=find(priors.gamma(:,s)==1);  
        end
        %slice.mean_beta=1000*ones(p,1);
        %slice.mean_beta=100*ones(p,1);
        
    otherwise
        % Estimate spatial precision from data
end

if slice.verbose
    fprintf('\n');
    disp(sprintf('%s priors for A',priors.A));
end

% Set b coefficients of Gamma priors
slice.b_alpha_prior=slice.mean_alpha./slice.c_alpha_prior;
slice.b_beta_prior=slice.mean_beta./slice.c_beta_prior;
slice.b_lambda_prior=slice.mean_lambda./slice.c_lambda_prior;

slice.priors=priors;