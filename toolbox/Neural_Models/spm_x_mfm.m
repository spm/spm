function [x] = spm_x_mfm(P)
% intialises a state structure for a mean field model
% FORMAT [x] = spm_x_mfm(P)
%
% P - parameter structure
% x - state structure
%__________________________________________________________________________

% dimensions
%--------------------------------------------------------------------------
ns   = size(P.A{1},1);                           % number of sources
np   = 3;                                        % number of populations


% create 
%--------------------------------------------------------------------------
for i = 1:ns
    for j = 1:np
        x{i,j}.V  = -8;
        x{i,j}.gE = 0;
        x{i,j}.gI = 0;
    end
end


% steady-state solution 
%==========================================================================

% create LFP model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_mfm';
M.x   = x;
M.pE  = P;
M.n   = length(spm_vec(x));
M.m   = size(P.C,2);
M.l   = size(P.C,1);

% solve for fixed point 
%--------------------------------------------------------------------------
U.u   = sparse(64,1);
U.dt  = 8/1000;
x     = spm_int_L(P,M,U);
x     = spm_unvec(x(end,:),M.x);
