function varargout = spm_shoot_greens(varargin)
% Build and apply FFT of Green's function (to map from momentum to velocity)
% FORMAT v = spm_shoot_greens(m,K,prm)
% m    - Momentum field n1*n2*n3*3 (single prec. float)
% K    - Fourier transform representation of Green's function
%        - either size n1*n2*n3 or n1*n2*n3*3*3
% prm  - Differential operator parameters (3 voxel sizes, 5 hyper-parameters)
%        - only needed when K is of size n1*n2*n3, in which case, voxel sizes
%          are necessary for dealing with each component individually
% v    - velocity field
%
% FORMAT K = spm_shoot_greens('kernel',dm,prm)
% dm  - dimensions n1*n2*n3
% prm - Differential operator parameters (3 voxel sizes, 5 hyper-parameters)
% K   - Fourier transform representation of Green's function
%        - either size n1*n2*n3 or n1*n2*n3*3*3
%
%________________________________________________________
% (c) Wellcome Trust Centre for NeuroImaging (2012)

% John Ashburner
% $Id: spm_shoot_greens.m 6008 2014-05-22 12:08:01Z john $

spm_diffeo('boundary',0);

if nargin==3 && isa(varargin{1},'char') && strcmp(varargin{1},'kernel'),
    d   = varargin{2};
    prm = varargin{3};

    F = spm_diffeo('kernel',d,prm);
    if size(F,4) == 1,
        % The differential operator is symmetric, so the Fourier transform should be real
        F = 1./real(fftn(F));
        if prm(4)==0, F(1,1,1) = 0; end;
    else
        for j=1:size(F,5),
            for i=1:size(F,4),
                % The differential operator is symmetric, so the Fourier transform should be real
                F(:,:,:,i,j) = real(fftn(F(:,:,:,i,j)));
            end
        end
        for k=1:size(F,3),
            % Compare the following with inverting a 3x3 matrix...
            A   = F(:,:,k,:,:);
            dt  = A(:,:,:,1,1).*(A(:,:,:,2,2).*A(:,:,:,3,3) - A(:,:,:,2,3).*A(:,:,:,3,2)) +...
                  A(:,:,:,1,2).*(A(:,:,:,2,3).*A(:,:,:,3,1) - A(:,:,:,2,1).*A(:,:,:,3,3)) +...
                  A(:,:,:,1,3).*(A(:,:,:,2,1).*A(:,:,:,3,2) - A(:,:,:,2,2).*A(:,:,:,3,1));
            msk     = find(dt==0);
            dt      = 1./dt;
            dt(msk) = 0;
            if prm(4)==0 && k==1, dt(1,1,1) = 0; end;
            F(:,:,k,1,1) = (A(:,:,:,2,2).*A(:,:,:,3,3) - A(:,:,:,2,3).*A(:,:,:,3,2)).*dt;
            F(:,:,k,2,1) = (A(:,:,:,2,3).*A(:,:,:,3,1) - A(:,:,:,2,1).*A(:,:,:,3,3)).*dt;
            F(:,:,k,3,1) = (A(:,:,:,2,1).*A(:,:,:,3,2) - A(:,:,:,2,2).*A(:,:,:,3,1)).*dt;

            F(:,:,k,1,2) = (A(:,:,:,1,3).*A(:,:,:,3,2) - A(:,:,:,1,2).*A(:,:,:,3,3)).*dt;
            F(:,:,k,2,2) = (A(:,:,:,1,1).*A(:,:,:,3,3) - A(:,:,:,1,3).*A(:,:,:,3,1)).*dt;
            F(:,:,k,3,2) = (A(:,:,:,1,2).*A(:,:,:,3,1) - A(:,:,:,1,1).*A(:,:,:,3,2)).*dt;

            F(:,:,k,1,3) = (A(:,:,:,1,2).*A(:,:,:,2,3) - A(:,:,:,1,3).*A(:,:,:,2,2)).*dt;
            F(:,:,k,2,3) = (A(:,:,:,1,3).*A(:,:,:,2,1) - A(:,:,:,1,1).*A(:,:,:,2,3)).*dt;
            F(:,:,k,3,3) = (A(:,:,:,1,1).*A(:,:,:,2,2) - A(:,:,:,1,2).*A(:,:,:,2,1)).*dt;
        end
    end
    varargout{1} = F;
else
    % Convolve with the Green's function via Fourier methods
    m = varargin{1};
    F = varargin{2};
    v = zeros(size(m),'single');
    if size(F,4) == 1,
        % Simple case where convolution is done one field at a time
        prm = varargin{3};
        for i=1:3,
            v(:,:,:,i) = ifftn(F.*fftn(m(:,:,:,i))*prm(i)^2,'symmetric');
        end
    else
        % More complicated case for dealing with linear elasticity, where
        % convolution is not done one field at a time
        for i=1:3,
            m(:,:,:,i) = fftn(m(:,:,:,i));
        end
        for k=1:size(m,3),
            a = m(:,:,k,:);
            m(:,:,k,:) = 0;
            for j=1:3,
                for i=1:3,
                    m(:,:,k,j) = m(:,:,k,j) + F(:,:,k,j,i).*a(:,:,:,i);
                end
            end
        end
        for i=1:3,
            v(:,:,:,i) = ifftn(m(:,:,:,i),'symmetric');
        end
    end
    varargout{1} = v;
end
%__________________________________________________________________________________

