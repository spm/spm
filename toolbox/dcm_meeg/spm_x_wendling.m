function [x] = spm_x_wendling(P)
% Initial state of a Wendling neural mass model
% FORMAT [x] = spm_x_wendling(P)
% P - parameters
%
% x        - x(0)
%
% 10 states per source:
%   x(:,1)  - excitatory PSP on pyramidal (y0)
%   x(:,2)  - excitatory PSP on excitatory interneurons (y1)
%   x(:,3)  - slow inhibitory PSP on pyramidal (y2)
%   x(:,4)  - excitatory PSP on slow inhibitory interneurons (y3)
%   x(:,5)  - fast inhibitory PSP on pyramidal (y4)
%   x(:,6:10) - derivatives of y0-y4
%__________________________________________________________________________
%
% Wendling F, Bartolomei F, Bellanger JJ, Chauvel P (2002)
%__________________________________________________________________________

% array of states
%--------------------------------------------------------------------------
n  = length(P.A{1});                          % number of sources
m  = 10;                                      % number of states
x  = sparse(n,m);
