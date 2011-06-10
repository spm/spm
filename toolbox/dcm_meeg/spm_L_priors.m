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
% $Id: spm_L_priors.m 4348 2011-06-10 20:50:23Z karl $

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
        pE.Lpos = dipfit.Lpos;   pC.Lpos = ones(3,n)*V;      % positions
        pE.L    = sparse(3,n);   pC.L    = ones(3,n)*exp(8); % orientations
        
    case{'IMG'}
        %------------------------------------------------------------------
        m       = dipfit.Nm;                                 % number modes
        pE.Lpos = sparse(3,0);   pC.Lpos = sparse(3,0);      % positions
        pE.L    = sparse(m,n);   pC.L    = ones(m,n)*exp(8); % modes
        
    case{'LFP'}
        %------------------------------------------------------------------
        pE.Lpos = sparse(3,0);   pC.Lpos = sparse(3,0);      % positions
        pE.L    = ones(1,m);     pC.L    = ones(1,m)*64;     % gains
        
    otherwise
        warndlg('Unknown spatial model')
        
end

% contributing states (encoded in J)
%==========================================================================
switch model
    
    case{'ERP','SEP'}
        %------------------------------------------------------------------
        pE.J = sparse(1,[1 7 9],[0.2 0.2 0.6],1,9);       % 9 states
        pC.J = pE.J/16;
        
    case{'CMC'}
        %------------------------------------------------------------------
        pE.J = sparse(1,[1 3 7],[0.2 0.8 0.2],1,8);       % 8 states
        pC.J = sparse(1,[3 7],1,1,8)/16;
        
    case{'LFP'}
        %------------------------------------------------------------------
        pE.J = sparse(1,[1 7 9],[0.2 0.2 0.6],1,13);      % 13 states
        pC.J = pE.J/16;
        
    case{'NMM'}
        %------------------------------------------------------------------
        pE.J = sparse(1,[1,2,3],[0.1 0.1 1],1,9);         % 9 states
        pC.J = sparse(1,[1,2],[1/32 1/128],1,9);
        
    case{'MFM'}
        %------------------------------------------------------------------
        pE.J = sparse(1,[1,2,3],[0.1 0.1 1],1,36);        % 36 states =
        pC.J = sparse(1,[1,2],[1/64 1/128],1,36);         % 9 1st + 27 2nd
        
    case{'DEM','NFM'}
        %------------------------------------------------------------------
        pE.J = [];                                        % null
        pC.J = [];
        
        
    otherwise
        warndlg('Unknown neural model')
        
end
