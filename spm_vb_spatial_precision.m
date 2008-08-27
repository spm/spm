function [D] = spm_vb_spatial_precision (vxyz,prior_type)
% Compute spatial precision matrix appropriate to prior
% FORMAT [D] = spm_vb_spatial_precision (vxyz,prior_type)
%
% vyyz          List of voxels
% prior_type    Type of prior
%
% D             Spatial precision matrix
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_spatial_precision.m 2022 2008-08-27 11:26:29Z lee $

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

    case 'Spatial - WGL',
        % Weighted graph-Laplacian, L. see Chung 1997 "Spectral graph theory"
        % and Strang 2007 "Computational Science and Engineering"
        % L = A'*C*A
        % where A and C are the edge-node incidence and "constitutive"
        % matrices resp. A has dimensions (no. edges,no.nodes). A and A' 
        % are the discrete analogues of grad and div operators, i.e. L=A'A 
        % is the standard unnormalized Laplacian (same as L in 'Spatial -
        % LORETA'). C is diagonal and contains edge weights, which here are
        % a function of the OLS estimates of regressors. See Eqn.8 in
        % Harrison et al 2008 NeuroImage 41 pp 408-423, but note that here
        % we DO NOT use the diffusion kernel, i.e. expm(-L*tau), as the
        % prior spatial covariance matrix, but L as the spatial precision
        % matrix, as in Penny et al 2005. 
        slice = vxyz;
        clear vxyz
        N = slice.N; % no. of nodes in the graph, i.e. voxels
        % norm of spatial gradient of regressors is w.r.t Hf 
        M = (1/N)*slice.wk_ols*ones(N);
        C = (slice.wk_ols-M)*(slice.wk_ols-M)';
        Hf = inv(C);
        % ka scales squared distance used for edge weights
        ka = 1;
        % get list of edges - from c to v
        vxyz = slice.Dw.vxyz;
        [r,c,v]=find(vxyz');
        edges = [c,v];
        % undirected graph, so only need store upper [lower] triangle
        i = find(edges(:,2) > edges(:,1));
        edges = edges(i,:);
        Ne = size(edges,1); % no. of edges in the graph
        % edge-node incidence matrix
        A = sparse([1:Ne,1:Ne],[edges(:,1),edges(:,2)],[ones(Ne,1),-ones(Ne,1)],Ne,N);
        % spatial gradients of ols estimates
        dB = slice.wk_ols*A';
        % edge weights
        dg2 = sum((dB'*Hf).*dB',2); % squared norm of spatial gradient of regressors
        ds2 = 1 + dg2; % squared distance in space is 1 as use only nearest neighbors
        weights = exp(-ds2/ka);
        % constitutive matrix 
        C = sparse(1:Ne,1:Ne,weights,Ne,Ne);
        % (un-normalized) weighted graph-Laplacian
        % Note this is equivalent to L = degree - (weight matrix), i.e.
        % W = sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],[weights;weights],N,N);
        % deg = sum(W,2) + eps;
        % L = sparse(1:N,1:N,deg,N,N) - W;
        L = A'*C*A;
        D = L;
        
    otherwise
        disp('Error in spm_vb_spatial_precision: unknown precision type');
        return
end