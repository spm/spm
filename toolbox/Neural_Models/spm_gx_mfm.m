function [m] = spm_gx_mfm(x,u,P,M)
% observer for a mean-field model (spiking)
% FORMAT [m] = spm_gx_mfm(x,u,P,M)
% x      - state vector
% m      - spiking activity (ns x np)
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_mfm.m 2941 2009-03-24 17:45:56Z maria $
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
if iscell(x)
    mfm = 1;                                    % mean-field model
else
    mfm = 0;
    x = {x};                                    % neural-mass model
end
ns   = size(x{1},1);                            % number of sources
np   = size(x{1},2);                            % number of populations
 
% Voltages
%--------------------------------------------------------------------------
VR   =  -40;                                   % threshold potential
 
% mean-field effects
%==========================================================================
if mfm
    
    Vx  = squeeze(x{2}(1,1,:,:));            % population variance (mV^2)
    Vx  = reshape(Vx,ns,np);                 % of voltage
   
else
    
    % neural-mass approximation to covariance of states
    %----------------------------------------------------------------------
    Cx = [   75.3843    0.1746    0.7487;
             0.1746     0.0040         0;
             0.7487          0    0.0160];
    Cx = exp(P.S)*Cx;
    Vx = Cx(1,1);  
    
end
 
% mean population firing and afferent extrinsic input
%--------------------------------------------------------------------------
m     = spm_Ncdf_jdw(x{1}(:,:,1),VR,Vx);  
