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
slice.update_alpha=1; 
slice.update_beta=1;

switch priors.WA,
    
case 'Spatial - GMRF',
    disp('');
    disp('Using GMRF');
    slice.D=spm_vb_geom (vxyz);
case 'Spatial - LORETA',
    disp('');
    disp('Using LORETA prior');
    L     = spm_vb_laplacian(vxyz);
    slice.D=L'*L;
case 'Voxel - Shrinkage',
    disp('');
    disp('Using Voxel Shrinkage priors for W and A');
    slice.D=speye(N);
case 'Voxel - Uninformative',
    disp('');
    disp('Using Voxel Uninformative priors for W and A');
    slice.D=speye(N);
    slice.mean_alpha=0.000001*ones(size(slice.X,2),1);
    slice.update_alpha=0;
    slice.mean_beta=0.000001*ones(size(slice.X,2),1);
    slice.update_beta=0;

end
  
if slice.p==0
    return
else
    slice.update_a=1;
end

if (priors.WA=='Spatial - GMRF') | (priors.WA=='Spatial - LORETA')                
    switch priors.overA
    case 'Voxel',
        % Very low (fixed) spatial precision ensures spatial prior is 
        % effectively ignored
        slice.update_beta=0;
        slice.mean_beta=1e-6*ones(slice.p,1);
        disp('Voxel-wise AR coeffs');
    case 'Slice',
        % Very high (fixed) spatial precision ensures AR coefficients are
        % (nearly) identical over the slice
        slice.update_beta=0;
        slice.mean_beta=1e6*ones(slice.p,1);
        disp('Slice-wise AR coeffs');
    otherwise
        % Estimate spatial precision from data
        return
    end
end