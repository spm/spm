function [S] = spm_vb_spatial_precision (prior_type,vxyz,img)
% Compute spatial precision matrix appropriate to prior
% FORMAT [S] = spm_vb_spatial_precision (vxyz,prior_type)
%
% vxyz          List of voxels
% prior_type    Type of prior
% img           used to construct weights of WGL
%
% S             Spatial precision matrix
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny, Nelson Trujillo-Barreto and Lee Harrison
% $Id: spm_vb_spatial_precision.m 2451 2008-11-10 16:20:32Z lee $

switch prior_type,
    
    case 'Spatial - UGL',
        % (un-normalized) Unweighted Graph Laplacian (UGL). 
        % L = A'*A, where A is the "incidence" matrix that has dimensions
        % (no. edges,no.nodes). A and A' are the discrete analogues of 
        % grad and div operators.
        % This is equaivalent to L = D - W, where W is the
        % "adjacency" matrix and D = diag(sum(W,2)). Note that in general W 
        % contains weighted edges (see WGL), but for UGL these are unity  
        N   = size(vxyz,1);
        [edges,weights] = spm_vb_edgeweights(vxyz); % weights equal to unity
        A   = spm_vb_incidence(edges,N);
        L   = A'*A; % UGL
        S   = L;
        
    case 'Spatial - GMRF',
        % Gaussian Markov Random Field with geometric boundary conditions
        % See equation in Woolrich et al. (ref 26 in paper VB2) 
        % Also described in dicussion in paper VB2
        % Warning: S for this prior can be singular !
        % This is equivalent to a normalized UGL, i.e. Z*UGL*Z, where
        % Z = D.^(-1/2), D = diag(sum(W,2)), and W is the adjacency matrix
        N   = size(vxyz,1);
        [edges,weights] = spm_vb_edgeweights(vxyz); % weights equal to unity
        A   = spm_vb_incidence(edges,N);
        L   = A'*A; %UGL
        W   = spm_vb_adjacency(edges,weights,N);
        Z   = spdiags((sum(W,2)+eps).^(-1/2),0,N,N);
        S   = Z*L*Z; % normalized UGL
        
    case 'Spatial - LORETA',
        % Unbiased LORETA PRIOR
        % Ensures normalisation is correct for edges/corners 
        % - see discussion in section 2 of paper VB2
        % this is equaivalent to a bi-Laplacian (using a UGL), i.e. UGL^2
        N   = size(vxyz,1);
        [edges,weights] = spm_vb_edgeweights(vxyz); % weights equal to unity
        Ne  = size(edges,1);
        A   = spm_vb_incidence(edges,N);
        L   = A'*A; % UGL
        S   = L'*L; % bi-Laplacian

    case 'Spatial - WGL',
        % Weighted graph-Laplacian, L. see Chung 1997 "Spectral graph theory"
        % and Strang 2007 "Computational Science and Engineering"
        % L = A'*C*A
        % where A and C are the edge-node "incidence" and "constitutive"
        % matrices resp. C is diagonal and contains edge weights, which here are
        % a function of the OLS estimates of regressors. See Eqn.8 in
        % Harrison et al 2008 NeuroImage 41 pp 408-423, but note that here
        % we DO NOT use the diffusion kernel, i.e. expm(-L*tau), as the
        % prior spatial covariance matrix, but L as the spatial precision
        % matrix, as in Penny et al 2005. 
        N   = size(vxyz,1);
        [edges,weights] = spm_vb_edgeweights(vxyz,img);
        Ne  = size(edges,1);
        A   = spm_vb_incidence(edges,N);
        C   = sparse(1:Ne,1:Ne,weights,Ne,Ne); % constitutive matrix 
        L   = A'*C*A; % (un-normalized) weighted graph-Laplacian
        S   = L;
       
    otherwise
        disp('Error in spm_vb_spatial_precision: unknown precision type');
        return
end