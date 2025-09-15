function [G, H] = spm_bsd_csd(P, M)
% Compute cross-spectral density (CSD) matrix G and optional modes/noise structure H
% for a Bayesian spectral density model.
%
% Inputs:
%   P - Parameter structure (contains spectral parameters)
%   M - Model structure (contains frequency info, dipfit info, etc.)
%
% Outputs:
%   G - Cross-spectral density matrix (nf x ns x ns)
%   H - (optional) Structure with spectral modes and noise components
%_______________________________________________________________________
% Copyright (C) 2024-2025 Wellcome Trust Centre for Neuroimaging

% Johan Medrano


    % Get dimensions
    %----------------------------------------------------------------------
    nf = length(M.Hz);           % Number of frequency bins
    ns = M.dipfit.Ns;            % Number of sources
    nk = size(P.f, 1);           % Number of spectral modes
    G = zeros(nf, ns, ns);       % Initialize CSD matrix

    % Convert parameters to spectral specification
    spec = spm_bsd_param2spec(P, M); 
    
    f = spec.freq;               % Mode frequencies
    a = spec.ampl;               % Mode amplitudes
    S = spec.fwhm.^2./(8 * log(2)); % Mode variances (from FWHM)
    b = spec.noise;              % Noise parameters

    % Compute mode basis functions (Hc)
    if M.sharefreqs 
        % Shared frequencies across sources
        Hc = (exp(-(M.Hz - f).^2./(2*S)));
    else
        % Separate frequencies for each source pair
        Hc = zeros(nk, nf, ns, ns); 
        for i = 1:ns
            for j = i:ns
                Hc(:, :, i, j) = (exp(-(M.Hz - f(:,i,j)).^2./(2*S(:,i,j))));
            end
        end
    end

    % Compute CSD matrix G
    for i = 1:ns
        for j = i:ns
            if M.sharefreqs
                % Shared frequencies
                Gc = Hc' * a(:,i,j) ... 
                + log(b(1,i,j)) - log(b(3,i,j)+ M.Hz'.^(b(2,i,j)));  
            else
                % Separate frequencies
                Gc =  Hc(:, :, i, j)'... 
                * a(:,i,j)  ...
                + log(b(1,i,j)) - log(b(3,i,j) + M.Hz'.^(b(2,i,j))); 
            end
            Gc = exp(Gc); % Exponentiate to get power

            G(:, i, j) = Gc;
            if i ~= j
                G(:, j, i) = conj(Gc); % Ensure Hermitian symmetry
            end
        end
    end
    

    % Optionally compute spectral modes and noise structure H
    if nargout > 1
        H = []; 
        H.modes = cell(1,size(f, 1)); 
        [H.modes{:}] = deal(zeros(nf, ns, ns)); 
        H.noise = zeros(nf, ns, ns); 
        for i = 1:ns
            for j = 1:ns 
                for k = 1:size(f, 1)
                    if M.sharefreqs
                        H.modes{k}(:,i,j) = Hc(k, :)' .* a(k,i,j); % Mode contribution
                    else
                        H.modes{k}(:,i,j) =  Hc(k, :, i, j)'.* a(k,i,j); % Mode contribution
                    end
                end
                H.noise(:, i, j) =  log(b(1,i,j)) - log(b(3,i,j) + M.Hz'.^(b(2,i,j))); % Noise contribution
                H = spm_unvec(exp(spm_vec(H)), H); % Reshape and exponentiate
            end
        end
    end

end
