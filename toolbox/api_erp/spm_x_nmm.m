function [x] = spm_x_nmm(P,GE,GI)
% initialises a state structure for a mean field model
% FORMAT [x] = spm_x_nmm(P,GE,GI)
%
% P - parameter structure
% GE - extrinsic connections (excitatory)
% GI - extrinsic connections (inhibitory)
%
% x        - array of states
% x(i,j,k) - k-th state of j-th population on i-th source
%
%   population: 1 - excitatory spiny stellate cells (input cells)
%               2 - inhibitory interneurons
%               3 - excitatory pyramidal cells      (output cells)
%
%        state: 1 V  - voltage
%               2 gE - conductance (excitatory)
%               3 gI - conductance (inhibitory)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_x_nmm.m 1228 2008-03-18 21:28:04Z karl $
 
% intrinsic connections (specifying the number of populations per source)
%==========================================================================
try, GE; catch
 
    % intrinsic connections (np x np) - excitatory
    %----------------------------------------------------------------------
    GE   = [0   0   1/2;
            0   0   1;
            1   0   0];
end
try, GI; catch
 
    % intrinsic connections (np x np) - inhibitory
    %----------------------------------------------------------------------
    GI   = [0   1/2 0;
            0   0   0;
            0   2   0];
end
 
% get initialisation from full mean-field model
%--------------------------------------------------------------------------
[x M] = spm_x_mfm(P,GE,GI);
 
% remove dispersion and fix the covariance of the states (Cx)
%--------------------------------------------------------------------------
M.x   = x{1};
M.Cx  = x{2}(:,:,1,1);
x     = x{1};
