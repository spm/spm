function [slice] = spm_vb_set_priors (slice,priors,vxyz)
% Set up precision parameters and optimisation options to correspond to priors
% FORMAT [slice] = spm_vb_set_priors (slice,priors,vxyz)
%
% slice     VB-GLM-AR data structure
% priors    For Regression and AR coefficients:
%
%           .WA='Spatial - GMRF' : W and A coeffs that are spatially regularised
%               'Spatial - LORETA' : as above but different spatial prior
%               'Voxel - Shrinkage' : W and A coeffs that are shrunk voxel-wise
%               'Voxel - Uninformative' : W and A coeffs without prior
%
%           For AR coefficients, override options:
%
%           .overA='Voxel' : Override spatial prior to have voxel-specific coeffs
%                 ='Slice' : Override spatial prior to have slice-specific coeffs
%
% vxyz      locations of voxels 
%
% Because slice.D, the spatial precision matrix, is the same for W and A,
% the options for the priors are currently limited to the above
%
% %W% Will Penny %E%

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
slice.mean_beta=ones(p,1);
slice.mean_lambda=ones(N,1);

switch priors.WA,
    % Get relevant spatial precision or roughness matrix
    
    case 'Spatial - GMRF',
        % Gaussian Markov Random Field with geometric boundary conditions
        % See equation in Woolrich et al. (ref 26 in paper VB2) 
        % Also described in dicussion in paper VB2
        % Warning: D for this prior can be singular !
        disp('\n');
        disp('Using GMRF');
        voxels=size(vxyz,1);
        slice.D=speye(voxels);
        for j=1:voxels,
            num_neighbors(j)=length(nonzeros(vxyz(j,:)));
        end
        for j=1:voxels,
            neighbors=nonzeros(vxyz(j,:));
            geom_means=sqrt(num_neighbors(j)*num_neighbors(neighbors));
            slice.D(j,neighbors)=-1./geom_means;
        end
        
%     case 'Spatial - LORETA',
%         % This produces the 
%         % Laplacian used in LORETA, with biased boundary conditions
%         % - see discussion in section 2 of paper VB2
%         disp('');
%         disp('Using LORETA prior');
%         voxels=size(vxyz,1);
%         L=4*speye(voxels);
%         for j=1:voxels,
%             L(j,nonzeros(vxyz(j,:)))=-1;
%         end
%         slice.D=L'*L;
        
    case 'Spatial - LORETA',
        % Unbiased LORETA PRIOR
        % Ensures normalisation is correct for edges/corners 
        % - see discussion in section 2 of paper VB2
        disp('\n');
        disp('Using LORETA prior');
        voxels=size(vxyz,1);
        L=4*speye(voxels);
        for j=1:voxels,
            L(j,nonzeros(vxyz(j,:)))=-1;
            L(j,j)=length(nonzeros(vxyz(j,:)));
        end
        slice.D=L'*L;
        
    case 'Spatial - RTPS',
        % Regularised Thin Plate Splines - Buckley, 1994.
        
    case 'Voxel - Shrinkage',
        disp('\n');
        disp('Using Voxel Shrinkage priors for W and A');
        slice.D=speye(N);
        
    case 'Voxel - Uninformative',
        disp('\n');
        disp('Using Voxel Uninformative priors for W and A');
        slice.D=speye(N);
        slice.mean_alpha=0.000001*ones(k,1);
        slice.update_alpha=0;
        slice.mean_beta=0.000001*ones(k,1);
        slice.update_beta=0;
        
end

if slice.p > 0
    slice.update_a=1;
end

if (strcmp(priors.WA,'Spatial - GMRF')) | (strcmp(priors.WA,'Spatial - LORETA'))
    switch priors.overA
        case 'Voxel',
            % Very low (fixed) spatial precision ensures spatial prior is 
            % effectively ignored
            slice.update_beta=0;
            slice.mean_beta=1e-4*ones(p,1);
            disp('Voxel-wise AR coeffs');
        case 'Slice',
            % Very high (fixed) spatial precision ensures AR coefficients are
            % (nearly) identical over the slice
            slice.update_beta=0;
            slice.mean_beta=1e4*ones(p,1);
            disp('Slice-wise AR coeffs');
        otherwise
            % Estimate spatial precision from data
            
    end
end

% Set b coefficients of Gamma priors
slice.b_alpha_prior=slice.mean_alpha./slice.c_alpha_prior;
slice.b_beta_prior=slice.mean_beta./slice.c_beta_prior;
slice.b_lambda_prior=slice.mean_lambda./slice.c_lambda_prior;

slice.priors=priors;