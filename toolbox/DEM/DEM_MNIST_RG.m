function mdp = DEM_MNIST_RG
% Demo of structure learning (i.e.,disentanglement) of MNIST digits
%__________________________________________________________________________
%
% This routine uses the MNIST digit classification problem to explore
% hierarchical generative models in the spirit of the renormalisation
% group. Its focus is to demonstrate the generation of image like content
% using discrete state spaces. Specifically, a hierarchical decomposition
% in which content and location are factorise. In other words, things can
% move around at any scale — and their movement of spatial transformations
% are removed before passing to the next level for grouping and implicit
% composition. This can be regarded as the R (reduction) and G (grouping)
% operators of an RG flow — a flow that converts place coded content into
% object coded latent states (with accompanying course graining of place
% information).
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_surveillance.m 8313 2022-09-30 18:33:43Z karl $


%% set up and preliminaries
%==========================================================================

% train
%==========================================================================
% Specify first hierarchical level
%--------------------------------------------------------------------------
spm_figure('GetWin','train'); clf;

% Primitives (basis functions - lines)
%--------------------------------------------------------------------------
% Load (pre-processed) MNIST dataset (uncomment to prepare spm_MNIST)
%--------------------------------------------------------------------------
cd 'C:\Users\karl\Dropbox\matlab'
load mnist
Y    = test.images(:,:,1:1000);
Y    = reshape(Y,28,28,1,1000);

n    = [28,28];
c    = 2;
d    = 5;

% window
%--------------------------------------------------------------------------
w    = exp(-((1:d) - mean((1:d))).^2/(2*(d/6)^2));
W    = diag(kron(w(:),w(:)));

% Displacement operators
%--------------------------------------------------------------------------
T     = {};
r     = [0,1,-1];
for i = 1:3
    for j = 1:3
        Di     = spm_speye(d,d,r(i),c);
        Dj     = spm_speye(d,d,r(j),c);
        T{i,j} = kron(Dj,Di)';
    end
end

% Hierarchical decomposition
%--------------------------------------------------------------------------
for h = 1:4

    % Indices of groups (partition)
    %----------------------------------------------------------------------
    ij    = {};
    ji    = {};
    for i = 0:d:(n(1) - d)
        ij{end + 1} = (1:d) + i;
    end
    for j = 0:d:(n(2) - d)
        ji{end + 1} = (1:d) + j;
    end


    % Grouping (G)
    %----------------------------------------------------------------------
    g     = {};
    for i = 1:numel(ij)
        for j = 1:numel(ji)
            for f = 1:size(Y,3)
                for t = 1:size(Y,4)
                    g{i,j}(:,f,t) = spm_vec(Y(ij{i},ji{j},f,t));
                end
            end
        end
    end
    
    a     = cell(size(g));
    V     = cell(size(g));
    for j = 1:numel(g)

        ng    = d*d;
        nf    = size(g{j},2);
        
        y     = reshape(g{j},ng*nf,[]);
        y     = kron(eye(nf,nf),W)*y;
        C     = y*y';

        % eigenmodes and spatial derivatives
        %------------------------------------------------------------------
        V     = cell(1,8);
        A     = cell(1,8);
        for i = 1:8

            S(i)   = trace(C);
            u      = spm_svd(C);
            U      = u(:,1);

            if max(U) < 0, U = -U; end
            V{i}   = U;

            for ii = 1:size(T,1)
                for jj = 1:size(T,1)
                    A{i}(:,ii,jj) = full(kron(eye(nf,nf),T{ii,jj})*U);
                end
            end

            % residual forming matrix
            %--------------------------------------------------------------
            R = A{i}(:,:);
            R = eye(size(R,1)) - R*pinv(R);
            C = R*C*R;

        end

        % Retain principle modes
        %------------------------------------------------------------------
        i    = S > max(S)*exp(-4);
        a{j,h} = A(i);
        v{j,h} = V(i);

        % eigenmodes
        %------------------------------------------------------------------
        for j = 1:sum(i)
            subplot(4,4,j)
            imagesc(reshape(V{j},d,d)),
            axis image, drawnow
        end

    end

end


% Subroutines
%==========================================================================


function spm_show_a(a)
% display learned images in likelihood tensors
% FORMAT spm_show_a(a)
% a{g} - likelihood tensors
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% illustrate images
%--------------------------------------------------------------------------
for n = 1:size(a{1},3)
    for u = 1:size(a{1},2)
        for g = 1:numel(a)
            Y{g}   = spm_dir_norm(a{g}(:,u,n));      % normalised
            % Y{g} = a{g}(:,u,n);                    % unnormalised
        end
        YY{1,u}  = spm_o2MNIST(Y);
    end
    subplot(10,1,n)
    if numel(YY) > 31
        imagesc(spm_cat(reshape(YY(1:32),2,16)))
    else
        imagesc(spm_cat(YY))
    end
    axis image, axis off
    drawnow
end

return

function [O,D] = spm_MNIST2o(mnist,id,d)
% convert image into a probabilistic output
% FORMAT [O,D] = spm_MNIST2o(mnist,id,[d])
% mnist - struct
% id    - content
% d     - location (1 to 25) [optional for mnist.mat images]
%
% O     - cell array of outcome probabilities
% D     - label (prior)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% create array of outcomes (one modality for each pixel)
%--------------------------------------------------------------------------
if isfield(mnist,'I')
    I = mnist.I(:,id);
    D = mnist.D(:,id);
    o = [I,1 - I]';
    O = num2cell(o,1);
    O = O(:);
    return
end


% If we are dealing with raw MNIST data, pre-process on the fly
%==========================================================================
h    = mnist.height;
w    = mnist.width;
n    = h*w;
I    = mnist.images(:,:,id);

% move image around, if requested
%--------------------------------------------------------------------------
if nargin > 2
    c    = zeros(5,5);
    c(d) = 1;
    I    = conv2(I,c,'same');
end

% smooth and histogram equalisation
%--------------------------------------------------------------------------
I     = spm_conv(I,2,2);                       % convolve
[m,i] = sort(I(:));                            % histogram equalisation
I(i)  = exp(-((1:n) - n).^2/(32*n));
for i = 1:h
    for j = 1:w
        o           = I(i,j);
        O{i,j}(1,1) = o;
        O{i,j}(2,1) = 1 - o;
    end
end

% return a cell array for this image
%--------------------------------------------------------------------------
O = O(:);
D = sparse(mnist.labels(id) + 1,1,1,10,1);

return


function [I] = spm_o2MNIST(Y)
% Converts a probabilistic output into an image
% FORMAT [I] = spm_o2MNIST(Y)
% Y      - posterior prediction
% I      - image
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% create array from (output) structures
%--------------------------------------------------------------------------
for i = 1:numel(Y)
    D(i) = Y{i}(1);
end

% assign pixel values
%--------------------------------------------------------------------------
[n,i] = sort(D(:));
n     = numel(D);
D(i)  = exp(-((1:n) - n).^2/(24*n));

% re-arrange into image
%--------------------------------------------------------------------------
if numel(D) ~= (28^2)
    load spm_MNIST pixels;
    I         = zeros(28,28);
    I(pixels) = D;
else
    I         = reshape(D,28,28);
end

return

