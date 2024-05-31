function [O,L,RGB] = spm_rgb2O(I,RGB)
% RGB image format to discrete outcomes (O)
% FORMAT [O,L,RGB] = spm_rgb2O(I);
% 
% I    - RGB image (pixels x pixels x 3 colours x time): uint8
%
% RGB.nd  - Diameter of tiles in pixels [default: 32]
% RGB.nb  - Number of discrete singular variates  [9: must be odd]  
% RGB.mm  - Maximum number of singular modes      [32]
% RGB.su  - Singular value (normalised) threshold [16]
% RGB.R   - temporal resampling                   [2]
%
% O    - Cell array of probabilistic outputs for pixel partition
% L    - Location of each probabilistic output modality
%
% RBG  - structure containing image size, pixel locations, singular vectors
%        and bins for singular variates
%
% RGB.M   - locations of groups
% RGB.N   - image size
% RGB.G   - cell array of outcome indices
% RGB.V   - singular vectors
% RGB.A   - singular variates
% RGB.O   - group index
%
% This routine effectively learns a basis set or receptive fields for a
% sequence of images supplied in RGB like format. It uses singular value
% decomposition, retaining the principal singular vectors for subsequent
% compression and reconstruction. The singular variates discretised over
% the range of each variate. The ensuing amortisation is saved in RGB
% structure for image reconstruction. This version uses overlapping
% receptive fields during the grouping and applies a Gaussian window around
% the centroid of each group location.
%
% Crucially, these spatial grouping operators can be generalised to
% spatiotemporal grouping by specifying a temporal resampling rate. In
% other words, one can discretise short image sequences.
%
% If the basis set has already been evaluated, this routine will return the
% discretised space-time voxels of any image timeseries, based upon the
% pre-evaluated spatial grouping and singular vectors.
%__________________________________________________________________________


% defaults
%--------------------------------------------------------------------------
try nd = RGB.nd; catch, nd = 32; end            % Diameter of tiles 
try nb = RGB.nb; catch, nb = 9;  end            % Number of variates  
try mm = RGB.mm; catch, mm = 32; end            % Maximum number of modes
try su = RGB.su; catch, su = 32; end            % Singular value threshold
try R  = RGB.R;  catch, R  = 2;  end            % temporal resampling

% ensure number of bins is odd (for symmetry)
%--------------------------------------------------------------------------
nb     = 2*fix(nb/2) + 1;
RGB.nb = nb;

% Permute for grouping
%--------------------------------------------------------------------------
T      = R*fix(size(I,4)/R);                    % length of time parition
I      = I(:,:,:,1:T);                          % truncate time series
I      = permute(I,[3,4,1,2]);                  % reshape time
N      = size(I,[1,2,3,4]);                     % size
I      = reshape(I,3*R,[],N(3),N(4));           % time 2nd
I      = permute(I,[2,1,3,4]);                  % time leading dimension
N      = size(I,[2,3,4]);                       % truecolours
T      = size(I,1);                             % time

% if V is provided
%--------------------------------------------------------------------------
if isfield(RGB,'V')

    % for each group
    %----------------------------------------------------------------------
    O     = {1,T};
    L     = [1,1];
    o     = 1;
    for g = 1:numel(RGB.G)

        % singular variates for this group
        %------------------------------------------------------------------
        Nm     = size(RGB.V{g},2);
        if Nm
            Y  = times(double(I(:,RGB.G{g})),RGB.W{g});
            u  = Y*RGB.V{g};
        end

        % generate (probability over discrete) outcomes
        %------------------------------------------------------------------
        for m = 1:Nm

            % dicretise singular variates
            %--------------------------------------------------------------
            for t = 1:T
                [~,U]  = min(abs(u(t,m) - RGB.A{o}));
                O{o,t} = sparse(U,1,1,RGB.nb,1);
            end

            % record locations and group for this outcome
            %--------------------------------------------------------------
            L(o,:)     = RGB.M(g,:);
            o          = o + 1;
        end
    end

    % discretization complete
    %----------------------------------------------------------------------
    return
end

% Spatial grouping
%--------------------------------------------------------------------------
L       = spm_combinations(N);
[G,M,W] = spm_tile(L(:,[2 3]),nd);

RGB.M = M;
RGB.N = N;
RGB.G = G;
RGB.W = W;
RGB.R = R;

O     = {1,T};
L     = [1,1];
o     = 1;
for g = 1:numel(G)

    % singular value decomposition for this group
    %----------------------------------------------------------------------
    Y        = times(double(I(:,G{g})),W{g});
    [u,s,v]  = spm_svd(Y,1/su);
    Nm       = min(length(s),mm);
    if Nm
        RGB.V{g} = v(:,1:Nm);
        u        = u(:,1:Nm)*s(1:Nm,1:Nm);
    end

    % generate (probability over discrete) outcomes
    %----------------------------------------------------------------------
    for m = 1:Nm
        
        % dicretise singular variates
        %------------------------------------------------------------------
        d     = max(abs(u(:,m)));
        a     = linspace(-d,d,nb);
        for t = 1:T
            [~,U]  = min(abs(u(t,m) - a));
            O{o,t} = sparse(U,1,1,nb,1);
        end

        % record locations and group for this outcome
        %------------------------------------------------------------------
        L(o,:)     = RGB.M(g,:);
        RGB.O(o)   = g;
        RGB.A{o}   = a;
        o          = o + 1;
    end
end

return
