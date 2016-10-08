function [g] = spm_gx_dem_salience(x,v,P)
% returns the prediction for visual search
% FORMAT [g] = spm_gx_dem_salience(x,v,P)
%
% x    - hidden states:
%   o(1) - oculomotor angle
%   o(2) - oculomotor angle
%   x(1) - relative amplitude of visual hypothesis 1
%   x(2) - relative amplitude of visual hypothesis 2
%   x(3) - ...
%
% v    - hidden causes
% P    - parameters
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception - x)
%   g(2) - oculomotor angle (proprioception - y)
%   g(3) - retinal input - channel 1
%   g(4) - retinal input - channel 2
%   g(5) - ...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_dem_salience.m 6901 2016-10-08 13:21:41Z karl $
 
 
% sensory input sampled from image
%--------------------------------------------------------------------------
global STIM

% retinotopic predictions
%--------------------------------------------------------------------------
s     = 0;
p     = spm_softmax(x.x(:));
for i = 1:min(length(STIM.H),length(x.x))
    s = s + p(i)*ADEM_sample_image(STIM.H{i},x.o,STIM.R);
end

 
% add proprioceptive angles in intrinsic coordinates
%--------------------------------------------------------------------------
g   = spm_vec(x.o,s);

