function [x] = spm_x_mfm(P)
% intialises a state structure for a mean field model
% FORMAT [x] = spm_x_mfm(P)
%
% P - parameter structure
% x - state structure
%__________________________________________________________________________

% dimnesions
%--------------------------------------------------------------------------
ns   = size(P.A{1},1);                            % number of sources
np   = 3;                                       % number of populations


% create 
%--------------------------------------------------------------------------
for i = 1:ns
    for j = 1:np
        x{i,j}.V  = -8;
        x{i,j}.gE = 0;
        x{i,j}.gI = 0;
    end
end
