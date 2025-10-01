function [pE, pC, pV] = spm_bsd_priors(fqs,ns,X,M,pE,pC)
% spm_bsd_priors.m
% Generate prior means, default values, and covariances for BSD model parameters.
% 
% Usage:
%   [pE, pC, pV] = spm_bsd_priors(fqs, ns, X, M, pE, pC)
%
% Inputs:
%   fqs - Frequency bands (cell array or numeric)
%   ns  - Number of sources
%   X   - External input matrix
%   M   - Model structure (must contain fields 'Hz' and 'sharefreqs')
%   pE  - (Optional) Prior means structure
%   pC  - (Optional) Prior covariances structure
%
% Outputs:
%   pE  - Prior means structure
%   pC  - Prior covariances structure
%   pV  - Default values structure
%
%_______________________________________________________________________
% Copyright (C) 2024-2025 Wellcome Trust Centre for Neuroimaging

% Johan Medrano


    % Set default values for pE and pC if not provided
    try pE; catch, pE = []; end
    try pC; catch, pC = []; end

    nf  = length(fqs); % Number of frequency bands

    % Get frequency means, default values, and band limits
    [f, S, W] = get_bsd_fS(fqs, M); 

    % Prior means, default values, and covariances for frequencies
    pE.f = 0.*f;
    pV.f = W;
    pC.f = zeros(nf, 1) + 1;
 
    % Prior means, default values, and covariances for width
    pE.S = f.*0; 
    pV.S = S; 
    pC.S = zeros(nf, 1) + 0.5; 
    
    % If frequencies are not shared, expand dimensions
    if ~M.sharefreqs
        pE.f = repmat(pE.f, 1, ns, ns); 
        pC.f = repmat(pC.f, 1, ns, ns); 
    
        pE.S = repmat(pE.S, 1, ns, ns); 
        pC.S = repmat(pC.S, 1, ns, ns); 
    end

    % Prior means and covariances for amplitude
    pE.a = zeros(nf, ns, ns) + 1; 
    pC.a = zeros(nf, ns, ns) + 1/16; 

    % Prior means and covariances for other parameters (b)
    pV.b = zeros(3, ns, ns) + [1; 2; 10];
    pE.b = 0.*pV.b;  
    pC.b = zeros(3, ns, ns) + [3; 3; 3] ; 

    % Zero out off-diagonal elements for non-shared frequencies and parameters
    for i = 1:ns
        for j = i+1:ns
            if ~M.sharefreqs
                pC.f(:, j, i) = 0;
                pC.S(:, j, i) = 0;
            end
            pC.a(:, j, i) = 0;
            pC.b(:, j, i) = 0;
        end
    end

    % Priors for B parameters (external inputs)
    nX = size(X, 2); 
    eB = cell(1,nX); 
    cB = cell(1,nX);

    [eB{:}] = deal(spm_unvec(0.*spm_vec(pE), pE));
    [cB{:}] = deal(spm_unvec(spm_vec(pC)*exp(4), pC)); 
    pE.B = eB; 
    pC.B = cB; 
end

% Helper function to extract frequency means, default values, and band limits
function [f, S, W] = get_bsd_fS(fqs, M)
    nf = length(fqs); 
    
    if iscell(fqs)
        f = []; 
        S = []; 
        W = zeros(nf, 2); 
        for i = 1:nf
            if isnumeric(fqs{i})
                bands = reshape(fqs{i}, 1, 2);
            else
                tks   = strsplit(fqs{i}, '+'); 
                bands = []; 
                
                for j = 1:numel(tks)
                    tk = strsplit(tks{j}, {'-', '_', '.', ':',',',';','/',' '});
                    if numel(tk) == 2 && ~isempty(intersect(tk{1}, {'high', 'low'}))
                        modifier = tk{1}; 
                        bandname = tk{2}; 
                    elseif numel(tk) == 1
                        modifier = ''; 
                        bandname = tk{1};
                    else 
                        error('Wrong frequency specifier: %s', fqs{i}); 
                    end
            
                    % Map band names to frequency ranges
                    switch(lower(bandname))
                        case 'delta'
                            band = [0 4]; 
                        case 'theta'
                            band = [4 8]; 
                        case {'alpha','mu'}
                            band = [8 12]; 
                        case 'beta'
                            band = [12 30]; 
                        case 'gamma'
                            band = [30 64];
                        otherwise
                            error('Frequency band not understood: %s', bandname); 
                    end
            
                    bwidth = band(2) - band(1);
                    
                    % Apply modifiers to band limits
                    switch modifier
                        case 'low'
                            band = [band(1) band(1)+bwidth/2]; 
                            bwidth = bwidth/2; 
                        case 'high'
                            band = [band(1)+bwidth/2 band(2)]; 
                            bwidth = bwidth/2;
                    end
                    if ~isempty(bands) 
                        if bands(1) ~= band(end) && bands(end) ~= band(1)
                            error('Frequency bands should be contiguous.')
                        end
    
                        bands = [min(band(1), bands(1)) max(band(end), bands(end))];
                    else 
                        bands = band;
                    end
    
                end
            end
        band = bands;

        % Compute mean and variance for the band
        bmean = (band(1) + band(2)) / 2; 
        bvar  = ((band(2) - band(1)))/16; 

        f = [f; bmean]; 
        S = [S; bvar]; 
        W(i,:) = band; 
    
        end
    end

    % Reshape outputs
    f = reshape(f, nf, 1); 
    S = reshape(S, nf, 1); 

    % Clip frequencies and bands to allowed range
    f = clip(f, min(M.Hz), max(M.Hz));
    W = clip(W, min(M.Hz), max(M.Hz));
end