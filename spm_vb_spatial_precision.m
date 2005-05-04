function [D] = spm_vb_spatial_precision (vxyz,prior_type)
% Compute spatial precision matrix appropriate to prior
% FORMAT [D] = spm_vb_spatial_precision (vxyz,prior_type)
%
% vyyz          List of voxels
% prior_type    Type of prior
%
% D             Spatial precision matrix
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_spatial_precision.m 112 2005-05-04 18:20:52Z john $

switch prior_type,
    
    case 'Spatial - GMRF',
        % Gaussian Markov Random Field with geometric boundary conditions
        % See equation in Woolrich et al. (ref 26 in paper VB2) 
        % Also described in dicussion in paper VB2
        % Warning: D for this prior can be singular !
        voxels=size(vxyz,1);
        D=speye(voxels);
        if voxels > 1
            for j=1:voxels,
                num_neighbors(j)=length(nonzeros(vxyz(j,:)));
            end
            for j=1:voxels,
                neighbors=nonzeros(vxyz(j,:));
                geom_means=sqrt(num_neighbors(j)*num_neighbors(neighbors));
                D(j,neighbors)=-1./geom_means;
            end
        end
        
    case 'Spatial - LORETA',
        % Unbiased LORETA PRIOR
        % Ensures normalisation is correct for edges/corners 
        % - see discussion in section 2 of paper VB2
        voxels=size(vxyz,1);
        max_neighbors=size(vxyz,2);
        L=max_neighbors*speye(voxels);
        for j=1:voxels,
            L(j,nonzeros(vxyz(j,:)))=-1;
            L(j,j)=length(nonzeros(vxyz(j,:)));
        end
        D=L'*L;
        
    case 'Spatial - LORETAu',
        % This produces the 
        % Laplacian used in LORETA, with biased boundary conditions
        % - see discussion in section 2 of paper VB2
        voxels=size(vxyz,1);
        max_neighbors=size(vxyz,2);
        L=max_neighbors*speye(voxels);
        for j=1:voxels,
            L(j,nonzeros(vxyz(j,:)))=-1;
        end
        D=L'*L;    
        
    otherwise
        disp('Error in spm_vb_spatial_precision: unknown precision type');
        return
end