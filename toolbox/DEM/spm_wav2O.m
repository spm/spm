function [O,L,WAV] = spm_wav2O(I,WAV)
% WAV image format to discrete outcomes (O)
% FORMAT [O,L,WAV] = spm_wav2O(I,WAV)
% 
% I    - WAV image (frequency bins x time)
%
% WAV.nd  - Size of tiles in frequency bins       [default: 4]
% WAV.nb  - Number of discrete singular variates  [9]  
% WAV.mm  - Maximum number of singular modes      [32]
% WAV.su  - Singular value (normalised) threshold [16]
% WAV.R   - temporal resampling                   [512]
%
% O    - Cell array of probabilistic outputs for partition
% L    - Location of each probabilistic output modality
%
% WAV  - structure containing image size, pixel locations, singular vectors
%        and bins for singular variates
%
% WAV.M   - locations of groups
% WAV.N   - image size
% WAV.G   - cell array of outcome indices
% WAV.V   - singular vectors
% WAV.A   - singular variates
% WAV.O   - group index
%
% This routine effectively learns a basis set or receptive fields for a
% sequence of images supplied in WAV like format. It uses singular value
% decomposition, retaining the principal singular vectors for subsequent
% compression and reconstruction. The singular variates discretised over
% the range of each variate. The ensuing amortisation is saved in WAV
% structure for sonogram reconstruction. This version uses overlapping
% receptive fields during the grouping and applies a Gaussian window around
% the centroid of each group location.
%
% Crucially, these spatial grouping operators can be generalised to
% spatiotemporal grouping by specifying a temporal resampling rate. In
% other words, one can discretise short sequences. This version is
% specialised for one-dimensional images; namely, frequency bins times
% time bins.
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
try nd = WAV.nd; catch, nd = 4;   end            % Diameter of tiles 
try nb = WAV.nb; catch, nb = 9;   end            % Number of variates  
try mm = WAV.mm; catch, mm = 32;  end            % Maximum number of modes
try su = WAV.su; catch, su = 16;  end            % Singular value threshold
try R  = WAV.R;  catch, R  = 512; end            % temporal resampling

% ensure number of bins is odd (for symmetry)
%--------------------------------------------------------------------------
nb     = 2*fix(nb/2) + 1;
WAV.nb = nb;

% permute for spatial grouping
%--------------------------------------------------------------------------
T      = R*fix(size(I,2)/R);                     % length of time parition
I      = I(:,1:T);                               % truncate time series
I      = permute(I,[2,1]);                       % reshape time
N      = size(I,[1,2]);                          % size
I      = reshape(I,R,[],N(2));                   % time 2nd dimension
I      = permute(I,[2,1,3]);                     % time leading dimension
N      = size(I,[2,3]);                          % voxels
T      = size(I,1);                              % time

% if V is provided
%--------------------------------------------------------------------------
if isfield(WAV,'V')

    % for each group
    %----------------------------------------------------------------------
    o     = 1;
    L     = 1;
    O     = {1,T};
    for g = 1:numel(WAV.G)

        % singular variates for this group
        %------------------------------------------------------------------
        Nm    = size(WAV.V{g},2);
        if Nm
            Y = times(double(I(:,WAV.G{g})),WAV.W{g});
            u = Y*WAV.V{g};
        end

        % generate (probability over discrete) outcomes
        %------------------------------------------------------------------
        for m = 1:Nm

            % dicretise singular variates
            %--------------------------------------------------------------
            for t = 1:T
                [~,U]  = min(abs(u(t,m) - WAV.A{o}));
                O{o,t} = sparse(U,1,1,WAV.nb,1);
            end

            % record locations and group for this outcome
            %--------------------------------------------------------------
            L(o,:)     = WAV.M(g,:);
            o          = o + 1;
        end
    end

    % discretization complete
    %----------------------------------------------------------------------
    return
end

% Frequency Grouping
%--------------------------------------------------------------------------
L       = spm_combinations(N);
[G,M,W] = spm_tile(L(:,2),nd);
WAV.M = M;
WAV.N = N;
WAV.G = G;
WAV.W = W;

o     = 1;
L     = 1;
O     = {1,T};
for g = 1:numel(G)

    % singular value decomposition for this group
    %----------------------------------------------------------------------
    Y        = times(double(I(:,G{g})),W{g});
    [u,s,v]  = spm_svd(Y,1/su);
    Nm       = min(length(s),mm);
    if Nm
        WAV.V{g} = v(:,1:Nm);
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
        L(o,:)     = WAV.M(g,:);
        WAV.O(o)   = g;
        WAV.A{o}   = a;
        o          = o + 1;
    end
end

return
