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
    if ~isfield(M, 'dipfit')
        
        % static lead field
        y = P.L*(x(:,9).*P.K);
        
    else
        % parameterised lead field
        
        % leadfield
        L = M.L;

        % re-compute lead field only if any dipoles changed
        Id = find(any([M.Lpos; M.Lmom] ~= [P.Lpos; P.Lmom]));

        if ~isempty(Id)
            for i = Id
                Lf = eeg_leadfield4(P.Lpos(:,i), M.dipfit.elc, M.dipfit.vol);
                L(:,i) = M.E*Lf*P.Lmom(:,i);
            end

            % store new leadfield and parameters
            M.L = L;
            M.Lpos = P.Lpos;
            M.Lmom = P.Lmom;
        end


        % output
        y = M.S'*L*(x(:, 9)*(2*10^4).*P.K);

    end
end
