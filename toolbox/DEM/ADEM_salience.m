function ADEM_salience
% Saccadic eye movements under active inference
%__________________________________________________________________________
% This demo illustrates exploration or visual search in terms of optimality
% principles based on straightforward ergodic or allostatic principles.
% In other words, to maintain the constancy of our external milieu, it is
% sufficient to expose ourselves to predicted and predictable stimuli.
% Being able to predict what is currently seen also enables us to predict
% fictive sensations that we will experience from another viewpoint. This
% provides a principled way in which to explore and sample the world for
% example with visual searches using saccadic eye movements. These
% theoretical considerations are remarkably consistent with a number
% of compelling heuristics; most notably the Infomax principle or the
% principle of minimum redundancy, signal detection theory and recent
% formulations of salience in terms of Bayesian surprise. The example
% here uses saliency (the posterior precision associated with fictive
% sampling of sensory data) to simulate saccadic eye movements under
% active inference.
%__________________________________________________________________________
% Copyright (C) 2011-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: ADEM_salience.m 6901 2016-10-08 13:21:41Z karl $


pth = fileparts(mfilename('fullpath'));


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



% mapp images and get hypotheses
%--------------------------------------------------------------------------
global STIM
DEMO = 0;

try
    
    STIM.H{1}   = spm_vol(fullfile(pth,'face_R.nii'));
    STIM.H{2}   = spm_vol(fullfile(pth,'face_rot_R.nii'));
    STIM.H{3}   = spm_vol(fullfile(pth,'face_inv_R.nii'));
    
    STIM.S{1}   = spm_vol(fullfile(pth,'face.nii'));
    STIM.S{2}   = spm_vol(fullfile(pth,'face_rot.nii'));
    STIM.S{3}   = spm_vol(fullfile(pth,'face_inv.nii'));
    
catch
    
    error('Images not found.');
    
end


% set-up:
%--------------------------------------------------------------------------
dim    = 16;                                  % dimension of visual sample
ns     = dim*dim;                             % number of sensory channels
nh     = length(STIM.H);                      % number of hypotheses
STIM.R = spm_hanning(dim)*spm_hanning(dim)';  % Retinal precision

STIM.V = spm_vol(fullfile(pth,'Nefertiti_R.nii')); % Stimulus (filtered)
STIM.U = spm_vol(fullfile(pth,'Nefertiti.nii'));   % Stimulus (unfiltered)

STIM.V = spm_vol(fullfile(pth,'face_R.nii'));      % Stimulus (filtered)
STIM.U = spm_vol(fullfile(pth,'face.nii'));        % Stimulus (unfiltered)


% hidden states
%--------------------------------------------------------------------------
x.o    = [0;0];                               % oculomotor angle
x.x    = -log(nh)*ones(nh,1);                 % hypotheses


% Recognition model
%==========================================================================
M(1).E.s = 1/2;                               % smoothness
M(1).E.n = 4;                                 % order of
M(1).E.d = 2;                                 % generalised motion


% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f  = 'spm_fx_dem_salience';              % plant dynamics
M(1).g  = 'spm_gx_dem_salience';              % prediction

M(1).x  = x;                                  % hidden states
M(1).V  = exp([8 8 4*ones(1,ns)]);            % error precision (g)
M(1).W  = exp(8);                             % error precision (f)


% level 2:
%--------------------------------------------------------------------------
M(2).v  = [0;0];                              % priors
M(2).V  = exp(16);


% generative model
%==========================================================================

% first level
%--------------------------------------------------------------------------
G(1).f  = 'spm_fx_adem_salience';
G(1).g  = 'spm_gx_adem_salience';
G(1).x  = [0;0];                              % hidden states
G(1).V  = exp(16);                            % error precision
G(1).W  = exp(8);                             % error precision
G(1).U  = [exp(8) exp(8) zeros(1,ns)];        % gain

% second level
%--------------------------------------------------------------------------
G(2).v  = 0;                                  % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = exp(16);


% generate and invert
%==========================================================================
N     = 16;                                   % length of data sequence
nr    = 32;                                   % size of salience map
a     = 1/2;                                  % autoregression for salience
s     = STIM.V.dim(1)/16/16;                  % width of IOR
R     = sparse(nr*nr,1);                      % salience map with IOR

DEM.G = G;
DEM.M = M;
DEM.C = sparse(1,N);
DEM.U = sparse(2,N);

if DEMO
    
    % (k) saccades
    %----------------------------------------------------------------------
    for k = 1:8
        
        % solve and save saccade
        %------------------------------------------------------------------
        DEM     = spm_ADEM(DEM);
        DEM     = spm_ADEM_update(DEM);
        
        % overlay true values
        %------------------------------------------------------------------
        spm_DEM_qU(DEM.qU,DEM.pU)
        
        % compute salience
        %------------------------------------------------------------------
        [S L]   = spm_salience_map(DEM.M,nr);
        
        % optimise prior belief
        %------------------------------------------------------------------
        S       = (S - min(S)).*(1 - R);
        [i j]   = max(S);
        DEM.U   = L(:,j)*ones(1,N);
        DEM.S   = reshape(S,nr,nr);
        
        % inhibition of return (IOR)
        %------------------------------------------------------------------
        D       = exp(-sum((L - L(:,j)*ones(1,nr*nr)).^2)/(2*s*s))';
        R       = a*R + D;
        
        % store
        %------------------------------------------------------------------
        ADEM{k} = DEM;
        
    end
    
    % save
    %----------------------------------------------------------------------
    save ADEM_saccades ADEM
    
end

% load
%--------------------------------------------------------------------------
load ADEM_saccades
    
    
% create movie in extrinsic and intrinsic coordinates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_dem_search_plot(ADEM(1:end))

% create movie in extrinsic and intrinsic coordinates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_dem_search_trajectory(ADEM)

% create movie in extrinsic and intrinsic coordinates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');
STIM.R = spm_hanning(32)*spm_hanning(32)';
spm_dem_search_movie(ADEM)


return


% illustrate salience for
%==========================================================================
nr          = 64;
STIM.H{1}   = spm_vol(fullfile(pth,'Nefertiti_R.nii'));
M(1).x.x    = 0;
S           = spm_salience_map(M,nr);

subplot(2,1,1)
imagesc(reshape(exp(S/6),nr,nr))
axis image
subplot(2,1,2)
imagesc(spm_read_vols(STIM.H{1}))
axis image



% create images for (memory mapped) sampling - for a given image F
%==========================================================================
fname = 'Nefertiti'
DIM = [size(F) 1];
M   = eye(4,4);
F   = F/max(max(F));

% write volume structure
%--------------------------------------------------------------------------
V   = struct(...
    'fname',  fullfile(pth,[fname '.nii']),...
    'dim',    DIM,...
    'mat',    M,...
    'pinfo',  [1 0 0]',...
    'descrip',fname);
V   = spm_create_vol(V);
V   = spm_write_plane(V,F,1);


% create images for (memory mapped) sampling - with Gabor filtering
%==========================================================================
for i = 1:length(H)
    
    fname = spm_file(H{i}.fname,'suffix','_R');
    s     = spm_read_vols(H{i});
    s     = s/max(max(s));
    s     = (spm_conv(s,1) - spm_conv(s,4));
    DIM   = [size(s) 1];
    M     = eye(4,4);
    
    % write volume structure
    %----------------------------------------------------------------------
    V   = struct(...
        'fname',  fname,...
        'dim',    DIM,...
        'mat',    M,...
        'pinfo',  [1 0 0]',...
        'descrip',spm_file(fname,'basename'));
    V   = spm_create_vol(V);
    V   = spm_write_plane(V,s,1);
    
end

