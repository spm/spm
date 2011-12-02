% Slow pursuit under active inference:
%__________________________________________________________________________
% This demo illustrates
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: ADEM_salience.m 4580 2011-12-02 20:22:19Z karl $


% hidden causes and states
%==========================================================================
% x    - hidden states:
%   o(1) - oculomotor angle
%   o(2) - oculomotor angle
%   x(1) - relative amplitude of visual hypothesis 1
%   x(2) - relative amplitude of visual hypothesis 2
%   x(3) - ...
%
% v    - hidden causes
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception - x)
%   g(2) - oculomotor angle (proprioception - y)
%   g(3) - retinal input - channel 1
%   g(4) - retinal input - channel 2
%   g(5) - ...
%--------------------------------------------------------------------------


% mapp images and get retinal precision
%--------------------------------------------------------------------------
global STIM
dim    = 16;
ns     = dim*dim;
nh     = 3;

H{1}   = spm_vol('face.nii');
H{2}   = spm_vol('face_null.nii');
H{3}   = spm_vol('face_inverted.nii');

ns     = dim*dim;                 % number of sensory (visual) channels
nh     = 3;                       % number of hypothesis

% rescale to a maximum of one
%--------------------------------------------------------------------------
for i = 1:nh
    H{i}.pinfo(1) = 1/max(max(spm_read_vols(H{i})));
end
R      = hanning(dim)*hanning(dim)';

STIM.V = H{1};                    % Stimulus
STIM.H = H;                       % Hypotheses
STIM.R = R;                       % Retinal precision


% parameters mapping from (unstable) point attractors to visual space
%--------------------------------------------------------------------------
x.o    = [0;0];                              % oculomotor angle
x.x    = -log(nh)*ones(nh,1);                % hypotheses


% Recognition model
%==========================================================================
M(1).E.s = 1/2;                               % smoothness
M(1).E.n = 4;                                 % order of
M(1).E.d = 2;                                 % generalised motion


% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f  = 'spm_fx_dem_salience';             % plant dynamics
M(1).g  = 'spm_gx_dem_salience';             % prediction

M(1).x  = x;                                 % hidden states
M(1).V  = diag(exp([6 6  4*ones(1,ns)]));    % error precision (g)
M(1).W  = diag(exp([8 8  4*ones(1,nh)]));    % error precision (f)

M(1).Ra = [1;2];                             % proprioceptive action


% level 2:
%--------------------------------------------------------------------------
M(2).v  = [0;0];                             % priors
M(2).V  = exp(16);


% generative model
%==========================================================================

% first level
%--------------------------------------------------------------------------
G(1).f  = 'spm_fx_adem_salience';
G(1).g  = 'spm_gx_adem_salience';
G(1).x  = [0;0];                               % hidden states
G(1).V  = exp(12);                             % error precision
G(1).W  = exp(12);                             % error precision

% second level
%--------------------------------------------------------------------------
G(2).v  = 0;                                   % exogenous forces
G(2).a  = [0; 0];                              % action forces
G(2).V  = exp(16);

% Salience map
%==========================================================================
DIM   = STIM.V.dim;
ng    = 32;
X     = linspace(-DIM(1)/2,DIM(1)/2,ng);
Y     = linspace(-DIM(2)/2,DIM(2)/2,ng);
[Y,X] = meshgrid(X,Y);
L     = [X(:) Y(:)]'/16;



% generate and invert
%==========================================================================
N     = 16;                                  % length of data sequence
DEM.G = G;
DEM.M = M;
DEM.C = sparse(1,N);
DEM.U = sparse(2,N);

for k = 1:4
    
    % solve and save saccade
    %----------------------------------------------------------------------
    DEM     = spm_ADEM(DEM);
    DEM     = spm_ADEM_update(DEM);
    
    % overlay true values
    %----------------------------------------------------------------------
    spm_DEM_qU(DEM.qU,DEM.pU)
    
    % ccompute salience
    %----------------------------------------------------------------------
    M     = DEM.M;
    for i = 1:length(L)
        M(2).v   = L(:,i);               % hidden cause (control)
        M(1).x.o = L(:,i);               % fictive hidden state
        qP       = spm_DEM_qC(M);
        S(i,1)   = spm_logdet(qP)/2;
    end
    
    % save
    %----------------------------------------------------------------------
    DEM.S   = reshape(S,ng,ng);
    ADEM{k} = DEM;
    
    % optimise and specify prior belief
    %----------------------------------------------------------------------
    [i j]   = max(S);
    DEM.U   = L(:,j)*ones(1,N);
    
end

% create movie in extrinsic and intrinsic coordinates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_dem_search_plot(ADEM)

% create movie in extrinsic and intrinsic coordinates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_dem_search_movie(ADEM)

return






% create images for spm (memory mapped) smapling
%==========================================================================
load DEM_IMG

F   = sum(F((1:128) + 32,(1:128) + 8,:),3);
DIM = [size(F) 1];
M   = eye(4,4);

%-Initialise new mask name: current mask & conditions on voxels
%----------------------------------------------------------------------
V   = struct(...
    'fname',  'face.nii',...
    'dim',    DIM,...
    'mat',    M,...
    'pinfo',  [1 0 0]',...
    'descrip','image of face');
V   = spm_create_vol(V);
V   = spm_write_plane(V,F,1);





