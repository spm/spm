function [y,DEM] = spm_SEM_gen(P,MM,U)
% Slow/saccadic eye movement prediction scheme
% FORMAT [y,DEM] = spm_SEM_gen(P,M,U)
%
%   P - parameters
%   M - (meta) model structure
%   U - trial-specific parameter deviates
%
%   y - {[ns,nx];...} - predictions for nx states {trials}
%                     - for ns samples (normalised lag)
%
% This smooth pursuit eye movement routine generates one cycle of motion
% under prior beliefs about a sinusoidal trajectory with variable phase.
%
% see also: spm_SEM_gen_full
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_SEM_gen.m 6014 2014-05-23 15:00:35Z guillaume $
 
% trial-specific initial states and parameters
%==========================================================================
if nargin > 2
    
    % cycle over different conditions
    %----------------------------------------------------------------------
    M       = MM;
    x       = M.x;
    [nu,nc] = size(U);
    y       = cell(1,nc);
    DEM     = cell(1,nc);
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
        [Y D]  = spm_SEM_gen(Q,M);
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
 
% angular frequency and phase (radians) of target motion (per time bin)
%--------------------------------------------------------------------------
N         = MM.ns;                               % length of data sequence
w         = 2*pi/N;                              % frequency (cycles/bin)
phi       = 0;
try, phi  = phase(MM.x.x(1) - 1j*MM.x.x(2)/w); end

% lag (in radians) of prior beliefs (U) about trajectory of attractor
%--------------------------------------------------------------------------
amp       = 1;
lag       = pi/16;
try, amp  = amp*exp(P.u(1)); end
try, lag  = lag*exp(P.u(2)); end
U         = amp*cos((1:N)*w + phi + lag);

% smoothness fluctuations
%--------------------------------------------------------------------------
s      = 1/4;
try, s = s*exp(P.s); end

% Experimental input (C)
%==========================================================================
try, C = MM.C; catch, C = cos((1:N)*w + phi); end


% Model set-up
%==========================================================================
M(1).E.s = s;                                    % smoothness
M(1).E.n = 4;                                    % order of
M(1).E.d = 1;                                    % generalised motion

% other parameters and precisions
%--------------------------------------------------------------------------
pE = [1/4 1/2 1/2 1/32 1/32 1/4];
hE = [4 4 4];

try, pE = pE + P.k; end
try, hE = hE + P.h; end

 
% occlusion and visual mapping functions
%--------------------------------------------------------------------------
occ = MM.occ;
vis = @(x) exp(-((-8:8)' - x.x(1) + x.o(1)).^2)*occ(x.x(1));
 
% forces (accelerations)
%--------------------------------------------------------------------------
a_x = @(x,v,P)(v - x.x(1))/4 - x.x(2)*P(6);
a_o = @(x,v,P)   - x.o(2)*P(2) +  ... 
              (v - x.o(1))*(P(1) - P(4)*(occ(v) | occ(x.x(1)))) + ...
         (x.x(1) - x.o(1))*(P(3) + P(5)*(occ(v) | occ(x.x(1))));

   

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
M(1).V  = exp(hE(1));                            % error precision
M(1).W  = exp(hE(2));                            % error precision
M(1).pE = pE;
 
 
% level 2: With hidden (memory) states
%--------------------------------------------------------------------------
M(2).v  = U(1);
M(2).V  = exp(hE(3));
 
% generative model (G) (with no noise)
%==========================================================================
 
% slow pursuit following with (second order) generative model
%--------------------------------------------------------------------------
x       = x.o;                                   % motor  angle & velocity
v       = C(1);                                  % target angle
 
% occlusion and visual mapping functions
%--------------------------------------------------------------------------
vis     = @(x,v) exp(-((-8:8)' - v + x(1)).^2)*occ(v);
 
% first level
%--------------------------------------------------------------------------
G(1).f  = @(x,v,a,P) [x(2); a - x(2)];
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
u      = DEM.pU.v{2}(1,:)';                      % target displacement
y      = y - u;                                  % lag
