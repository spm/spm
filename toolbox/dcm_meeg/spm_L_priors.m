function [pE,pC] = spm_L_priors(dipfit)
% prior moments for the lead-field parameters of ERP models
% FORMAT [pE,pC] = spm_L_priors(dipfit)
%
% dipfit    - forward model structure:
%
%    dipfit.type - 'ECD', 'LFP' or 'IMG'
%    dipfit.symm - distance (mm) for symmetry constraints (ECD)
%    dipfit.Lpos - x,y,z source positions (mm)            (ECD)
%    dipfit.Nm   - number of modes                        (IMG)
%    dipfit.Ns   - number of sources
%    dipfit.Nc   - number of channels           
%
% pE - prior expectation
% pC - prior covariance
%
% spatial parameters
%--------------------------------------------------------------------------
%    pE.Lpos - position                    - ECD
%    pE.L    - orientation                 - ECD
%              coefficients of local modes - Imaging
%              gain of electrodes          - LFP
%
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_L_priors.m 2208 2008-09-26 18:57:39Z karl $
 
% defaults
%--------------------------------------------------------------------------
try, type = dipfit.type; catch, type = 'LFP'; end
try, symm = dipfit.symm; catch, symm = 0;     end

% number of sources
%--------------------------------------------------------------------------
try
    n = dipfit.Ns;
    m = dipfit.Nc;
catch
    n = dipfit;
    m = n;
end
 
% parameters for electromagnetic forward model
%==========================================================================
switch type
 
    case{'ECD'}
        %------------------------------------------------------------------
        pE.Lpos = dipfit.Lpos;   Lpos = 0*eye(3*n);     % positions
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
 
% symmetry constraints (based on Euclidean distance from mirror image)
%==========================================================================
switch type

    case{'ECD'}

        % correlation
        %------------------------------------------------------------------
        Mpos  = [-pE.Lpos(1,:); pE.Lpos(2,:); pE.Lpos(3,:)];
        for i = 1:n
            for j = 1:n
                D(i,j) = sqrt(sum(pE.Lpos(:,i) - Mpos(:,j)).^2);
            end
        end
        D     = (D + D')/2;
        D     = D < symm;
        D     = D - diag(diag(D));
        L     = L - kron(D,diag([1 0 0])) + kron(D,diag([0 1 1]));
end

% prior covariance
%--------------------------------------------------------------------------
pC  = spm_cat(diag({Lpos, exp(8)*L}));
