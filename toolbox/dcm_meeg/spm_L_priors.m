function [pE,pC] = spm_L_priors(dipfit,pE,pC)
% prior moments for the lead-field parameters of ERP models
% FORMAT [pE,pC] = spm_L_priors(dipfit)
%
% dipfit    - forward model structure:
%
%    dipfit.type     - 'ECD', 'LFP' or 'IMG'
%    dipfit.symmetry - distance (mm) for symmetry constraints (ECD)
%    dipfit.location - allow changes in source location       (ECD)
%    dipfit.Lpos     - x,y,z source positions (mm)            (ECD)
%    dipfit.Nm       - number of modes                        (IMG)
%    dipfit.Ns       - number of sources
%    dipfit.Nc       - number of channels
%
% pE - prior expectation
% pC - prior covariance
%
% adds spatial parameters
%--------------------------------------------------------------------------
%    pE.Lpos - position                    - ECD
%    pE.L    - orientation                 - ECD
%              coefficients of local modes - Imaging
%              gain of electrodes          - LFP
%    pE.J    - contributing states (length(J) = number of states per source
%
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_L_priors.m 3738 2010-02-25 13:19:07Z vladimir $

% defaults
%--------------------------------------------------------------------------
try, model    = dipfit.model;    catch, model    = 'LFP'; end
try, type     = dipfit.type;     catch, type     = 'LFP'; end
try, location = dipfit.location; catch, location = 0;     end
try, symmetry = dipfit.symmetry; catch, symmetry = 0;     end
try, pC;                         catch, pC       = [];    end


% number of sources
%--------------------------------------------------------------------------
try
    n = dipfit.Ns;
    m = dipfit.Nc;
catch
    n = dipfit;
    m = n;
end

% location priors (4 mm)
%--------------------------------------------------------------------------
if location, V = 2^2; else, V = 0; end

% parameters for electromagnetic forward model
%==========================================================================
switch type
    
    case{'ECD'} % mean           and variance
        %------------------------------------------------------------------
        pE.Lpos = dipfit.Lpos;   Lpos = V*eye(3*n);     % positions
        pE.L    = sparse(3,n);   L    =   eye(3*n);     % orientations
        
    case{'IMG'}
        %------------------------------------------------------------------
        m       = dipfit.Nm;                            % number of modes
        pE.Lpos = sparse(3,0);   Lpos =   eye(0,0);     % positions
        pE.L    = sparse(m,n);   L    =   eye(m*n);     % modes
        
    case{'LFP'}
        %------------------------------------------------------------------
        pE.Lpos = sparse(3,0);   Lpos =   eye(0,0);     % positions
        pE.L    = ones(1,m);     L    =   eye(m,m);     % gains
        
    otherwise
        warndlg('Unknown spatial model')
        
end

% contributing states (encoded in J)
%==========================================================================
switch model
    
    case{'ERP','SEP'}
        %------------------------------------------------------------------
        pE.J = sparse(1,[1 7 9],[0.2 0.2 0.6],1,9);       % 9 states
        J    = diag(pE.J/64);
        
    case{'LFP'}
        %------------------------------------------------------------------
        pE.J = sparse(1,[1 7 9],[0.2 0.2 0.6],1,13);      % 13 states
        J    = diag(pE.J/64);
        
    case{'NMM'}
        %------------------------------------------------------------------
        pE.J = sparse(1,[1,2,3],[0.1 0.1 1],1,9);          % 9 states
        J    = sparse([1,2],[1,2],[1/128 1/128],9,9);
        
    case{'MFM'}
        %------------------------------------------------------------------
        pE.J = sparse(1,[1,2,3],[0.1 0.1 1],1,36);         % 9 states
        J    = sparse([1,2],[1,2],[1/128 1/128],36,36);    % 27 covariances
        
    case{'DEM'}
        %------------------------------------------------------------------
        pE.J = [];                                         % null
        J    = [];
        
    otherwise
        warndlg('Unknown neural model')
        
end


% Distance between homolgous sources (16mm)
%--------------------------------------------------------------------------
if symmetry, V = 16; else V = 0; end

% symmetry constraints (based on Euclidean distance from mirror image)
%==========================================================================
switch type
    
    case{'ECD'}
        if symmetry
            
            % correlation
            %------------------------------------------------------------------
            Mpos  = [-pE.Lpos(1,:); pE.Lpos(2,:); pE.Lpos(3,:)];
            D = Inf*ones(n);
            for i = 1:n
                for j = 1:n
                    if sign(pE.Lpos(1,i)) == sign(Mpos(1, j))
                        D(i,j) = sqrt(sum(pE.Lpos(:,i) - Mpos(:,j)).^2);
                    end
                end
            end
            D     = (D + D')/2;
            DD = zeros(n);
            for i = 1:n
                [M, I] = min(D(i, :));
                if M<V                   
                    DD(i, I) = 1;
                end
            end
            
            L     = L - kron(DD,diag([1 0 0])) + kron(DD,diag([0 1 1]));
        end
end

% prior covariance
%--------------------------------------------------------------------------
pC  = spm_cat(spm_diag({pC, Lpos, exp(8)*L, J}));
