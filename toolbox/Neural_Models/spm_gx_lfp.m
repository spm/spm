function [y] = spm_gx_lfp(x,u,P)
% observer for a neural mass model of erps
% FORMAT [y] = spm_gx_lfp(x,u,P)
% x      - state vector
%   x(:,1) - voltage (spiny stellate cells)
%   x(:,2) - voltage (pyramidal cells) +ve
%   x(:,3) - voltage (pyramidal cells) -ve
%   x(:,4) - current (spiny stellate cells)    depolarizing
%   x(:,5) - current (pyramidal cells)         depolarizing
%   x(:,6) - current (pyramidal cells)         hyerpolarizing
%   x(:,7) - voltage (inhibitory interneurons)
%   x(:,8) - current (inhibitory interneurons) depolarizing
%   x(:,9) - voltage (pyramidal cells)
%
% y        - measured voltage
%__________________________________________________________________________
%
% This is a simplified version of spm_gx_erp
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%

global M

% get dimensions and configure state variables
%--------------------------------------------------------------------------
m    = length(x);
n    = length(P.A{1});
x    = reshape(x,n,m);

% slect mixture of neuronal states providing signal for each source
% default - pyramidal cell depolarization
%--------------------------------------------------------------------------
try
    P.M;
catch
    P.M = sparse(9,1,1,m,1);                   
end
x  = x*P.M;

% output (mixtures of sources specified by lead field)
%==========================================================================
try
    M.Spatial_type;
catch
    M.Spatial_type = 3;                   
end

% fixed lead field (note; L has already been rotated by S: P.L = S*L)
%--------------------------------------------------------------------------
if M.Spatial_type == 3

    y = P.L*(x.*P.K);

% parameterised lead field ECD/EEG
%--------------------------------------------------------------------------
elseif M.Spatial_type == 1

    % re-compute lead field only if any dipoles changed
    %----------------------------------------------------------------------
    Id = find(any([M.Lpos; M.Lmom] ~= [P.Lpos; P.Lmom]));

    if length(Id)

        % record new spatial parameters
        %------------------------------------------------------------------
        M.Lpos = P.Lpos;
        M.Lmom = P.Lmom;

        % note that in all DCM code, EEG coordinates are in MNI
        % space-orientation. Transform here to fieldtrip.
        % transformation matrix from MNI-oriented coordinate system
        % to fieldtrip
        %------------------------------------------------------------------
        iMt = [0 1 0 20;
              -1 0 0 0;
               0 0 1 10;
               0 0 0 1];
        iSt = [0 1 0 0;
              -1 0 0 0;
               0 0 1 0;
               0 0 0 1];
        P.Lpos = iMt*[P.Lpos; ones(1, size(P.Lpos, 2))];
        P.Lmom = iSt*[P.Lmom; ones(1, size(P.Lmom, 2))];
        P.Lpos = P.Lpos(1:3, :); P.Lmom = P.Lmom(1:3, :);

        for i = Id
            Lf = fieldtrip_eeg_leadfield4(P.Lpos(:,i), M.dipfit.elc, M.dipfit.vol);
            M.L(:,i) = M.E*Lf*P.Lmom(:,i);
        end
    end

    % output
    %----------------------------------------------------------------------
    y = M.S'*L*(x(:,9)*(2*10^4).*P.K);

% parameterised lead field ECD/MEG
%--------------------------------------------------------------------------
elseif M.Spatial_type == 2

    % re-compute lead field only if any dipoles changed
    %----------------------------------------------------------------------
    Id = find(any([M.Lpos; M.Lmom] ~= [P.Lpos; P.Lmom]));

    if length(Id)
        for i = Id
            Lf = fieldtrip_meg_leadfield(P.Lpos(:,i)', M.grad, M.dipfit.vol);
            M.L(:,i) = M.E*Lf(M.Ichannels, :)*P.Lmom(:,i);
        end

        % record new spatial parameters
        %------------------------------------------------------------------
        M.Lpos = P.Lpos;
        M.Lmom = P.Lmom;
    end

    % output
    %----------------------------------------------------------------------
    y = M.S'*L*(x(:,9)*(10^20).*P.K);

end
