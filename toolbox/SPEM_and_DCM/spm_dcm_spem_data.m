function xY = spm_dcm_spem_data(xY)
% Prepare (decimate and normalise) data DCM for SPEM
% FORMAT xY = spm_dcm_spem_data(xY)
%
%   xY.Y{i}  - original data
%   xY.C{i}  - original target
%   xY.DT    - original timing
%
% creates:
%
%   xY.y{i} - normalised (decimated) lag (data - target)
%   xY.u{i} - normalised (decimated) target
%   xY.R(i) - decimation
%   xY.x(i) - intial states
%   xY.dt   - mean normalised (decimated) timing
%
%  This auxiliary routine  decimates and normalises eye movement data to a
%  single period of a (negative) cosine wave - of unit amplitude.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_spem_data.m 6014 2014-05-23 15:00:35Z guillaume $
 
% Get frequency (in bins) and target trajectory
%==========================================================================
N  = 64;                        % number of time bins after decimation
 
% cycle over conditions or trials
%--------------------------------------------------------------------------
for i = 1:length(xY.Y)
    
    % decimate
    %----------------------------------------------------------------------
    R         = ceil(length(xY.Y{i})/N);
    Cmax      = max(abs(xY.C{i}));
    xY.y{i}   = decimate(spm_vec(xY.Y{i})/Cmax,R);
    xY.u{i}   = decimate(spm_vec(xY.C{i})/Cmax,R);
    
    % take differences (lag in relation to target)
    %----------------------------------------------------------------------
    xY.y{i}   = xY.y{i} - xY.u{i};
    
    % initial states and velocities (o - occulomotor; x - target)
    %----------------------------------------------------------------------
    xY.x(i).o = [xY.y{i}(1) + xY.u{i}(1); xY.y{i}(2) - xY.y{i}(1)];
    xY.x(i).x = [xY.u{i}(1)             ; xY.u{i}(2) - xY.u{i}(1)];
    
    % time bins
    %----------------------------------------------------------------------
    xY.dt(i)  = xY.DT*R;
    xY.R(i)   = R;
 
end

% average time bin
%--------------------------------------------------------------------------
xY.dt = mean(xY.dt);
