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
%    'ECD'
%    'Imaging'
%    'LFP'
%
% determines whether the model is ECD or not. For imaging reconstructions
% the paramters P.L are a (m x n) matrix of coefficients that scale the
% contrition of n sources to m = M.dipfit.Nm modes encoded in M.dipfit.G.
%
% For LFP models (the default) P.L (1 x n) vector simply encodes the
% electrode gain for each of n sources.
%
% see; Kiebel et al. (2006) NeuroImage
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_erp_L.m 1535 2008-05-01 17:08:22Z vladimir $

% Create a persient variable that rembers the last locations
%--------------------------------------------------------------------------
persistent LastLpos LastL

% output
%==========================================================================

% number of sources - n
%--------------------------------------------------------------------------
n  = size(P.L,2);

% type of spatial model and modality
%--------------------------------------------------------------------------
try, type     = M.dipfit.type;     catch, type     = 'LFP'; end
try, modality = M.dipfit.modality; catch, modality = 'LFP'; end

switch type

    % parameterised lead field - ECD
    %----------------------------------------------------------------------
    case{'ECD'}

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

        for i = Id
            Lf = forwinv_compute_leadfield(P.Lpos(:,i), M.dipfit.sens, M.dipfit.vol);
            LastL(:,:,i) = Lf;
        end
        G    = spm_cond_units(LastL);
        for i = 1:n
            L(:,i) = G(:,:,i)*P.L(:,i);
        end
              
        % Imaging solution {specified in M.dipfit.G}
        %----------------------------------------------------------------------
    case{'Imaging'}

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

        % LFP or sources - no lead field required
        %----------------------------------------------------------------------
    case{'LFP'}

        L = sparse(diag(P.L));

    otherwise
        warndlg('unknown type of model')
end

