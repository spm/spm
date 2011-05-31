function [L] = spm_erp_L(P,M)
% returns [projected] lead field L as a function of position and moments
% FORMAT [L] = spm_erp_L(P,M)
% P  - model parameters
% M  - model specification
% L  - lead field
%__________________________________________________________________________
%
% The lead field (L) is constructed using the specific parameters in P and,
% where necessary information in the dipole structure M.dipfit. For ECD
% models P.Lpos and P.L encode the position and moments of the ECD. The
% field M.dipfit.type:
%
%    'ECD', 'LFP' or 'IMG'
%
% determines whether the model is ECD or not. For imaging reconstructions
% the paramters P.L are a (m x n) matrix of coefficients that scale the
% contrition of n sources to m = M.dipfit.Nm modes encoded in M.dipfit.G.
%
% For LFP models (the default) P.L simply encodes the electrode gain for 
% each source contributing a LFP.
%
% see; Kiebel et al. (2006) NeuroImage
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_erp_L.m 4336 2011-05-31 16:49:34Z rosalyn $

% Create a persient variable that rembers the last locations
%--------------------------------------------------------------------------
persistent LastLpos LastL


% type of spatial model and modality
%==========================================================================
try,   type = M.dipfit.type;  catch, type = 'LFP'; end

switch type

    % parameterised lead field - ECD
    %----------------------------------------------------------------------
    case{'ECD'}
        
        % number of sources - n
        %------------------------------------------------------------------
        n  = size(P.L,2);

        % re-compute lead field only if any dipoles changed
        %----------------------------------------------------------
        try
            Id = find(any(LastLpos ~= P.Lpos));
        catch
            Id = 1:n;
        end

        % record new spatial parameters
        %----------------------------------------------------------
        LastLpos = P.Lpos;
        for i  = Id
            if any(P.Lpos(:,i)>=200)
                Lf = zeros(M.dipfit.Nc, 3);
            else
                Lf = ft_compute_leadfield(transform_points(M.dipfit.datareg.fromMNI, P.Lpos(:,i)'), M.dipfit.sens, M.dipfit.vol);
            end
            LastL(:,:,i) = Lf;
        end
        G     = spm_cond_units(LastL);
        for i = 1:n
            L(:,i) = G(:,:,i)*P.L(:,i);
        end
              
    % Imaging solution {specified in M.dipfit.G}
    %----------------------------------------------------------------------
    case{'IMG'}
        
        % number of sources - n
        %------------------------------------------------------------------
        n  = size(P.L,2);

        % re-compute lead field only if any coeficients have changed
        %------------------------------------------------------------------
        try
            Id = find(any(LastLpos ~= P.L));
        catch
            Id = 1:n;
        end
        for i = Id
            LastL(:,i) = M.dipfit.G{i}*P.L(:,i);
        end

        % record new spatial parameters
        %------------------------------------------------------------------
        LastLpos = P.L;
        L        = LastL;

    % LFP electrode gain
    %----------------------------------------------------------------------
    case{'LFP'}
        m     = length(P.L);
        try
            n = M.dipfit.Ns;
        catch
            n = m;
        end
       L     = sparse(1:m,1:m,P.L,m,n);

    otherwise
        warndlg('unknown spatial model')
end

% -------------------------------------------------------------------------
function new = transform_points(M, old)
old(:,4) = 1;
new = old * M';
new = new(:,1:3);
