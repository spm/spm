function [y] = spm_gx_erp(x,u,P)
% observer for a neural mass model of erps
% FORMAT [y] = spm_gx_erp(x,u,P)
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
%___________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%___________________________________________________________________________
% %W% Karl Friston %E%

global M

% get dimensions and configure state variables
%---------------------------------------------------------------------------
x     = x(2:end);
n     = M.Nareas;
x     = reshape(x,n,9);

% configure parameters: Af,...
%---------------------------------------------------------------------------
P    = spm_erp_pack(P,M.m,n);

% output (sources if Lead field is empty)
%===========================================================================
if isempty(M.L)
    y = x(:,9);
else
    if M.Spatial_type == 3
        
        % fixed lead field
        y = M.S'*P.L*(x(:,9).*P.K);
        
    elseif M.Spatial_type == 1
        % parameterised lead field ECD/EEG
        
        % leadfield
        L = M.L;

        % re-compute lead field only if any dipoles changed
        Id = find(any([M.Lpos; M.Lmom] ~= [P.Lpos; P.Lmom]));

        if ~isempty(Id)
            % store new parameters
            M.Lpos = P.Lpos;
            M.Lmom = P.Lmom;
            
            % note that in all DCM code, EEG coordinates are in MNI
            % space-orientation. Transform here to fieldtrip.

            % transformation matrix from MNI-oriented coordinate system
            % to fieldtrip
            iMt = [[0 1 0 20]; [-1 0 0 0]; [0 0 1 10]; [0 0 0 1]];
            iSt = [[0 1 0 0]; [-1 0 0 0]; [0 0 1 0]; [0 0 0 1]];
            P.Lpos = iMt*[P.Lpos; ones(1, size(P.Lpos, 2))];
            P.Lmom = iSt*[P.Lmom; ones(1, size(P.Lmom, 2))];
            P.Lpos = P.Lpos(1:3, :); P.Lmom = P.Lmom(1:3, :);

            for i = Id
                
                Lf = fieldtrip_eeg_leadfield4(P.Lpos(:,i), M.dipfit.elc, M.dipfit.vol);
                L(:,i) = M.E*Lf*P.Lmom(:,i);
            end

            % store new leadfield
            M.L = L;
        end


        % output
        y = M.S'*L*(x(:, 9)*(2*10^4).*P.K);
    elseif M.Spatial_type == 2
        % parameterised lead field ECD/MEG
        L = M.L;

        % re-compute lead field only if any dipoles changed
        Id = find(any([M.Lpos; M.Lmom] ~= [P.Lpos; P.Lmom]));

        if ~isempty(Id)
            for i = Id
                Lf = fieldtrip_meg_leadfield(P.Lpos(:,i)', M.grad, M.dipfit.vol);
                L(:,i) = M.E*Lf(M.Ichannels, :)*P.Lmom(:,i);
            end

            % store new leadfield and parameters
            M.L = L;
            M.Lpos = P.Lpos;
            M.Lmom = P.Lmom;
        end

        % output
        y = M.S'*L*(x(:, 9)*(10^20).*P.K);

        
    end
end
