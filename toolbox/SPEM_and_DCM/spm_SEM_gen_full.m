function [y,DEM] = spm_SEM_gen_full(P,MM,U)
% Slow/saccadic eye movement prediction scheme - for model
% FORMAT [y,DEM] = spm_SEM_gen_full(P,M,U)
%
%   P - parameters
%   M - (meta) model structure
%   U - trial-specific parameter deviates
%
%   y   - {[ns,nx];...} - predictions for nx states {trials}
%                       - for ns samples (normalised lag)
%
% This generative routine is the same as spm_SEM_gen but includes an extra
% hierarchical level to infer the phase of underlying target motion. this
% sort of generative model is required when characterising violation or
% omission responses due to departures from the expected trajectory.
%
% see also: spm_SEM_gen
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_SEM_gen_full.m 6014 2014-05-23 15:00:35Z guillaume $
 
% trial-specific initial states and parameters
%==========================================================================
if nargin > 2
    
    % cycle over different conditions
    %----------------------------------------------------------------------
    M       = MM;
    x       = M.x;
    [nu,nc] = size(U);
    y       = cell(nu);
    DEM     = cell(nu);
    for i = 1:nc
        
        % target location
        %------------------------------------------------------------------
        M.C   = M.u{i};
        
        % initial states
        %------------------------------------------------------------------
        M.x   = x(i);

        % condition-specific parameters (encoded in U)
        %------------------------------------------------------------------
        Q = spm_vec(P.A);
        for j = 1:nu
            Q = Q + U(j,i)*spm_vec(P.B{j});
        end
        Q      = spm_unvec(Q,P.A);
        [Y,D]  = spm_SEM_gen(Q,M);
        y{i}   = Y;
        DEM{i} = D;
        
    end
    return
        
end
 
 
% INFERENCE MODEL: DEM.M and DEM.G
%==========================================================================
 
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor velocity
%   x.x(1) - target angle    - extrinsic coordinates
%   x.x(2) - target velocity - extrinsic coordinates
%
% v    - causal states: force on target
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception)
%   g(2) - oculomotor velocity
%   g(:) - visual input - intrinsic coordinates
%--------------------------------------------------------------------------
 
 
% Parameters and precisions
%==========================================================================
 
% Experimental input (C) and (correct) priors (U)
%==========================================================================
 
% angular frequency of target motion (per time bin)
%--------------------------------------------------------------------------
N  = MM.ns;                                   % length of data sequence
w  = 2*pi/N;                                  % frequency (cycles/bin)
 
% Sine wave prior and cause
%--------------------------------------------------------------------------
U  = sparse(1,N) + w;                         % frequency of target motion
try
    C   = MM.C;
catch
    C   = cos((1:N)*w);                       % sinusoidal target motion
end
 
% Phase
%--------------------------------------------------------------------------
try
    phi = phase(MM.x.x(1) - 1j*MM.x.x(2)/w);    % phase (radians) of target
catch
    phi = 0;
end
 
 
% lag (in radians) and smoothness of random fluctuations
%--------------------------------------------------------------------------
try, lag = exp(P.l)*2*pi/32; catch, lag = 2*pi/32; end
try, s   = exp(P.s)/2;       catch, s   = 1/2;     end
 
% other parameters and precisions
%--------------------------------------------------------------------------
pE.m = 4;
pE.a = 4;
pE.b = 4;
pE.c = 32;
pE.d = 1;
hE.u = 4;
hE.v = 4;
hE.w = 3;
 
try
    pE.m = pE.m.*exp(P.m);
    pE.a = pE.a.*exp(P.a);
    pE.b = pE.b.*exp(P.b);
    pE.c = pE.c.*exp(P.c);
    pE.c = pE.c.*exp(P.c);
    hE.u = hE.u.*exp(P.u);
    hE.v = hE.v.*exp(P.v);
    hE.w = hE.w.*exp(P.w);
end
 
% Model set-up
%==========================================================================
M(1).E.s = s;                                    % smoothness
M(1).E.n = 4;                                    % order of
M(1).E.d = 1;                                    % generalised motion
 
% occlusion and visual mapping functions
%--------------------------------------------------------------------------
occ = MM.occ;
vis = @(x) exp(-((-8:8)' - x.x(1) + x.o(1)).^2)*occ(x.x(1));
 
% forces (accelerations)
%----------------------------------------------------------------------
a_o = @(x,v,P)(v - x.o(1))/P.m - x.o(2)/P.b - x.o(1)/P.c + occ(v)*(x.x(1) - x.o(1))/P.a; % occulomotor
a_x = @(x,v,P)(v - x.x(1))/P.m - x.x(2)/P.d;     % target
 
 
% Generative model (M) (sinusoidal movement)
%==========================================================================
% Endow the model with internal dynamics (a simple oscillator) so that is
% recognises and remembers the trajectory to anticipate jumps in rectified
% sinusoidal motion.
 
% slow pursuit following with (second order) generative model
%--------------------------------------------------------------------------
try
    x   = MM.x;
catch
    x.o = [cos(phi); -w*sin(phi)];               % motor  angle & velocity
    x.x = [cos(phi); -w*sin(phi)];               % target angle & velocity
end
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f  = @(x,v,P) [x.o(2); a_o(x,v,P); x.x(2); a_x(x,v,P)];
M(1).g  = @(x,v,P) [x.o; vis(x)];
M(1).x  = x;                                     % hidden states
M(1).V  = exp(hE.u);                             % error precision
M(1).W  = exp(hE.w);                             % error precision
M(1).pE = pE;
 
 
% level 2: With hidden (memory) states
%--------------------------------------------------------------------------
M(2).f  = @(x,v,P)[x(2); -x(1)]*v;
M(2).g  = @(x,v,P) x(1);
M(2).x  = [cos((phi + lag)); -sin((phi + lag))]; % hidden states
M(2).V  = exp(hE.v);                             % error precision
M(2).W  = exp(16);                               % error precision
 
% level 3: Encoding frequency of memory states (U)
%--------------------------------------------------------------------------
M(3).v  = U(1);
M(3).V  = exp(16);
 
 
% generative model (G) (with no noise)
%==========================================================================
 
% slow pursuit following with (second order) generative model
%--------------------------------------------------------------------------
x   = x.o;                                       % motor  angle & velocity
v   = C(1);                                      % target angle
 
% occlusion and visual mapping functions
%--------------------------------------------------------------------------
vis = @(x,v) exp(-((-8:8)' - v + x(1)).^2)*occ(x(1));
 
% first level
%--------------------------------------------------------------------------
G(1).f  = @(x,v,a,P) [x(2); a/4 - x(2)];
G(1).g  = @(x,v,a,P) [x; vis(x,v)];
G(1).x  = x;                                     % hidden states
G(1).U  = exp([8 8 -ones(1,19)*16]);             % motor gain
 
 
% second level
%--------------------------------------------------------------------------
G(2).v  = v;                                     % exogenous force
G(2).a  = 0;                                     % action force
 
 
% Check generative model
%--------------------------------------------------------------------------
M       = spm_DEM_M_set(M);
 
 
% Generate prediction
%==========================================================================
 
% Bayesian inference (spm_DEM, spm_LAP, spm_ADEM or spm_ALAP)
%--------------------------------------------------------------------------
DEM.db = 0;
DEM.M  = M;
DEM.G  = G;
DEM.U  = U;
DEM.C  = C;
DEM    = spm_ADEM(DEM);
 
% (Hidden) behavioural or electrophysiological response (EXAMPLE)
%--------------------------------------------------------------------------
y      = DEM.pU.x{1}(1,:)';                      % ocular displacement
try
    u  = DEM.pU.x{1}(3,:)';                      % target displacement
catch
    u  = DEM.pU.v{2}(1,:)';                      % target displacement
end
y      = (y - u)/max(u);                         % normalised lag
