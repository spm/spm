function [block] = spm_vb_set_priors(block,priors,vxyz)
% Set up precision parameters and optimisation options to correspond to priors
% FORMAT [block] = spm_vb_set_priors(block,priors,vxyz)
%
% block   - VB-GLM-AR data structure (see spm_vb_glmar)
% priors  - For regression coefficients:
%
%           .W= 'Spatial - UGL' : coeffs that are spatially regularised
%               'Spatial - GMRF' : as above but different spatial prior
%               'Spatial - LORETA' : as above but different spatial prior
%               'Spatial - WGL' : as above but different spatial prior
%               'Voxel - Shrinkage' : that are shrunk voxel-wise
%               'Voxel - Uninformative' : coeffs without prior
%
%           For AR coefficients:
%
%           .A= 'Spatial - UGL' : coeffs that are spatially regularised
%               'Spatial - GMRF' : as above but different spatial prior
%               'Spatial - LORETA' : as above but different spatial prior
%               'Voxel - Shrinkage' : that are shrunk voxel-wise
%               'Voxel - Uninformative' : coeffs without prior
%               'Voxel - Limiting' : Voxel-specific coeffs as limiting case
%                                    of a LORETA prior
%               'block - Limiting' : block-specific coeffs as limiting case
%                                    of a LORETA prior
%               'Discrete' : Different coeffs as function of grey/white/CSF
%                            or other masks
%
% vxyz    - locations of voxels 
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Lee Harrison
% $Id: spm_vb_set_priors.m 6079 2014-06-30 18:25:37Z spm $

if ~isfield(block,'verbose')
    block.verbose = 0;
end

N = size(vxyz,1);
k = size(block.X,2);
p = block.p;

block.update_alpha = 1; 
block.update_beta  = 1;

% Set c coefficients of gamma priors and default mean
block.c_alpha_prior  = 0.1*ones(k,1);
block.c_beta_prior   = 0.1*ones(p,1);
block.c_lambda_prior = 0.1*ones(N,1);

block.mean_alpha     = ones(k,1);
block.mean_beta      = 1000*ones(p,1);
block.mean_lambda    = ones(N,1);

switch priors.W
    case {'Spatial - UGL','Spatial - GMRF','Spatial - LORETA'}
        block.Dw = spm_vb_spatial_precision(priors.W,vxyz);  
  
    case 'Spatial - WGL'
        % spm_vb_spatial_precision called in spm_vb_init_block 
        % as Dw requires wk_ols
        block.Dw.vxyz = vxyz;    
        
    case 'Voxel - Shrinkage'
        block.Dw = speye(N);
        block.update_alpha = 1;
        
    case 'Voxel - Uninformative'
        block.Dw=speye(N);
        block.mean_alpha=0.000001*ones(k,1);
        block.update_alpha=0;
end

if block.verbose
    fprintf('\n');
    fprintf('%s priors for W\n',priors.W);
end

if block.p > 0
    block.update_a = 1;
end

switch priors.A

    case {'Spatial - UGL','Spatial - GMRF','Spatial - LORETA'}
        block.Da = spm_vb_spatial_precision (priors.A,vxyz);

    case 'Voxel - Shrinkage'
        block.Da = speye(N);
        block.update_beta = 1;

    case 'Voxel - Uninformative'
        block.Da = speye(N);
        block.mean_beta   = 0.000001*ones(p,1);
        block.update_beta = 0;

    case 'Voxel - Limiting'
        % Very low (fixed) spatial precision ensures spatial prior is
        % effectively ignored
        block.Da = spm_vb_spatial_precision ('Spatial - LORETA',vxyz);
        block.update_beta = 0;
        block.mean_beta = 1e-4*ones(p,1);

    case 'block - Limiting'
        % Very high (fixed) spatial precision ensures AR coefficients are
        % (nearly) identical over the block
        block.Da = spm_vb_spatial_precision ('Spatial - LORETA',vxyz);
        block.update_beta = 0;
        block.mean_beta = 1e4*ones(p,1);

    case 'Discrete'
        fprintf('Different AR coeff priors for masked regions\n');
        if ~isfield(priors,'gamma')
            error('Must enter class labels');
        end
        priors.S = size(priors.gamma,2);
        priors.N = sum(priors.gamma(:,:));
        for s=1:priors.S
            priors.voxel(s).i = find(priors.gamma(:,s)==1);
        end
        %block.mean_beta = 1000*ones(p,1);
        %block.mean_beta = 100*ones(p,1);

    otherwise
        % Estimate spatial precision from data
end

if block.verbose
    fprintf('\n');
    fprintf('%s priors for A\n',priors.A);
end

% Set b coefficients of Gamma priors
block.b_alpha_prior  = block.mean_alpha ./ block.c_alpha_prior;
block.b_beta_prior   = block.mean_beta ./ block.c_beta_prior;
block.b_lambda_prior = block.mean_lambda ./ block.c_lambda_prior;

block.priors = priors;
