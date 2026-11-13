function [gg] = spm_MDP_conv(n,d)
% Multidimensional dot (inner) product
% FORMAT [gg] = spm_MDP_conv(n,d)
%
% n(1,2) - Size of image matrix
% d      - Number of pixel displacements
%
% gg{g}  - Cell array of output indices for each modality g
%
% This illustrative routine returns a cell array of output indices that map
% from modalities (i.e., the likelihood tensors) to outputs. This
% effectively allows one to implement a simple form of weight sharing by
% changing the mapping from modalities to outputs in a systematic way.
% Here, we consider translations as exemplar affine mappings from
% modalities to outputs.
% 
% The resulting cell array is then used during inversion to estimate the
% latent states encoding a fine transformations. The dimensions of the cell
% array (gg) correspond to orthogonal transformations and, implicitly, the
% shape of domain states (ff) on which they depend. By using an explicit
% geometry the transformations are strictly divergence; in the sense that
% one modality can generate multiple outputs; however, each transformation
% or mapping generates the same number of outputs (and is thereby informed
% by the same number).
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging


% Preliminaries
%==========================================================================
if nargin < 1, n = [32,32]; end
if nargin < 2, d = 4;       end

% Create Euclidean coordinates
%--------------------------------------------------------------------------
[x,y]  = ndgrid(1:n(1),1:n(2));
x      = x - mean(x,'all');
y      = y - mean(y,'all');
X      = [x(:), y(:)];
X(:,3) = 1;            % Coordinates (augmented with a constant)


% Shape of implicit domain factors
%--------------------------------------------------------------------------
f{1}  = -d:d;          % translation (x)
f{2}  = -d:d;          % translation (y)
for i = 1:numel(f)
    s{i} = numel(f{i});
end

% One cell tensor for each modality
%--------------------------------------------------------------------------
Ng    = size(X,1);
gg    = cell(Ng,1);
for g = 1:Ng
    gg{g} = cell(s{:});
end

% Assemble domain cell tensor of affine transformations
%==========================================================================
for f1 = 1:s{1}
    for f2 = 1:s{2}

        % translation (x)
        %----------------------------------------------------------
        M{2}      = eye(3,3);
        M{2}(3,1) = f{1}(f1);

        % translation (y)
        %----------------------------------------------------------
        M{1}      = eye(3,3);
        M{1}(3,2) = f{2}(f2);

        % affine mapping
        %----------------------------------------------------------
        Y  = X;
        for i = 1:numel(M)
            Y = Y*M{i};
        end

        % discrete mapping (closest neighbour)
        %----------------------------------------------------------
        for g = 1:Ng

            % Find nearest neighbours
            %------------------------------------------------------
            D     = minus(X,Y(g,:));
            [~,o] = min(sum(D.^2,2));

            % save in domain tensor
            %------------------------------------------------------
            gg{o}{f1,f2}(end + 1) = g;

        end

    end
end

% Primitives (basis functions - lines)
%--------------------------------------------------------------------------
if nargout
    return
end

% Illustrate affine transformations, with a circle
%==========================================================================
spm_figure('GetWin','Affine transformations'); clf;

% Create an image
%--------------------------------------------------------------------------
I     = zeros(n(1),n(2));
for i = 1:n(2)
    for j = 1:n(2)
        if (i - n(1)/2)^2 + (j - n(1)/2)^2 > (n(1)/3)^2
            I(i,j) = 1;
        end
    end
end

% And show set of transformations
%--------------------------------------------------------------------------
for f = 1:numel(gg{1})
    J     = NaN(size(I));
    for g = 1:Ng
        if numel(gg{g}{f})
            J(gg{g}{f}) = I(g);
        end
    end

    subplot(2,2,1)
    imagesc(reshape(J(:),n(1),n(2))), axis image
    drawnow
end

return

