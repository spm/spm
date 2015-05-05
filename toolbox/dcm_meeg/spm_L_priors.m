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
% $Id: spm_L_priors.m 6427 2015-05-05 15:42:35Z karl $



% defaults
%--------------------------------------------------------------------------
try, model    = dipfit.model;    catch, model    = 'LFP'; end
try, type     = dipfit.type;     catch, type     = 'LFP'; end
try, location = dipfit.location; catch, location = 0;     end
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
        pE.Lpos = dipfit.Lpos;   pC.Lpos = ones(3,n)*V;    % positions
        pE.L    = zeros(3,n);    pC.L    = ones(3,n)*64;   % orientations
        
       Sc = find(~cellfun(@isempty,dipfit.silent_source)); % silence sources for CSD 
       if(Sc); pC.L(:,Sc) = pC.L(:,Sc)*0; end
        
    case{'IMG'}
        %------------------------------------------------------------------
        m       = dipfit.Nm;                               % number modes
        pE.Lpos = sparse(3,0);   pC.Lpos = sparse(3,0);    % positions
        pE.L    = zeros(m,n);    pC.L    = ones(m,n)*64;   % modes
        
       Sc = find(~cellfun(@isempty,dipfit.silent_source)); % silence sources for CSD 
       if(Sc); pC.L(:,Sc) = pC.L(:,Sc)*0; end
        
    case{'LFP'}
        %------------------------------------------------------------------
        pE.Lpos = sparse(3,0);   pC.Lpos = sparse(3,0);    % positions
        pE.L    = ones(1,m);     pC.L    = ones(1,m)*64;   % gains
        
    otherwise
        warndlg('Unknown spatial model')
        
end

% contributing states (encoded in J)
%==========================================================================
switch upper(model)
    
    case{'ERP','SEP'}
        %------------------------------------------------------------------
        pE.J = sparse(1,9,1,1,9);               % 9 states
        pC.J = sparse(1,[1 7],1/32,1,9);
        
    case{'CMC','TFM'}
        %------------------------------------------------------------------
        pE.J = sparse(1,3,1,1,8);               % 8 states
        pC.J = sparse(1,[1 7],1/32,1,8);
        
    case{'LFP'}
        %------------------------------------------------------------------
        pE.J = sparse(1,9,1,1,13);              % 13 states
        pC.J = sparse(1,[1 7],1/32,1,13);
        
    case{'NMM'}
        %------------------------------------------------------------------
        pE.J = sparse(1,3,1,1,9);               % 9 states
        pC.J = sparse(1,[1,2],1/32,1,9);
        
   case{'NMDA'}
        %------------------------------------------------------------------
        pE.J = sparse(1,3,1,1,12);               % 12 states
        pC.J = sparse(1,[1,2],1/32,1,12);
        
    case{'CMM'}
        %------------------------------------------------------------------
        pE.J = sparse(1,2,1,1,12);              % 12 states
        pC.J = sparse(1,[3,4],1/32,1,12);
        
    case{'CMM_NMDA'}
        %------------------------------------------------------------------
        pE.J = sparse(1,2,1,1,16);              % 12 states
        pC.J = sparse(1,[3,4],1/32,1,16);
        
    case{'MFM'}
        %------------------------------------------------------------------
        pE.J = sparse(1,3,1,1,36);              % 36 (9 + 27) states
        pC.J = sparse(1,[1,2],1/32,1,36);
        
    case{'DEM','NFM'}
        %------------------------------------------------------------------
        pE.J = [];                              % null
        pC.J = [];
        
    case{'NULL'}
        %------------------------------------------------------------------
        nx   = size(pE.A,1)/m;
        pE.J = ones(1,nx);                      % nx states
        pC.J = ones(1,nx);
        
    otherwise
        warndlg('Unknown neural model')
        
end
