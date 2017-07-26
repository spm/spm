function [L] = spm_erp_L(P,dipfit)
% returns [projected] lead field L as a function of position and moments
% FORMAT [L] = spm_erp_L(P,dipfit)
% P       - model parameters
% dipfit  - spatial model specification
% L       - lead field
%__________________________________________________________________________
%
% The lead field (L) is constructed using the specific parameters in P and,
% where necessary information in the dipole structure dipfit. For ECD
% models P.Lpos and P.L encode the position and moments of the ECD. The
% field dipfit.type:
%
%    'ECD', 'LFP' or 'IMG'
%
% determines whether the model is ECD or not. For imaging reconstructions
% the paramters P.L are a (m x n) matrix of coefficients that scale the
% contrition of n sources to m = dipfit.Nm modes encoded in dipfit.G.
%
% For LFP models (the default) P.L simply encodes the electrode gain for 
% each source contributing a LFP.
%
% see; Kiebel et al. (2006) NeuroImage
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_erp_L.m 7142 2017-07-26 20:38:45Z karl $

% Create a persient variable that rembers the last locations
%--------------------------------------------------------------------------
persistent LastLpos LastL 


% type of spatial model and modality
%==========================================================================
if isfield(dipfit,'type'), type = dipfit.type; else, type = 'LFP'; end

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

        M = dipfit.datareg.fromMNI;
        
        if isfield(dipfit, 'siunits') && dipfit.siunits
            M = diag([1e-3 1e-3 1e-3 1])*M;
        end
        
        % record new spatial parameters
        %----------------------------------------------------------
        LastLpos = P.Lpos;
        for i  = Id
            if any(P.Lpos(:,i) >= 200)
                Lf = zeros(dipfit.Nc, 3);
            else
                Lf = ft_compute_leadfield(transform_points(M, P.Lpos(:,i)'), dipfit.sens, dipfit.vol);
            end
            LastL(:,:,i) = Lf;
        end
        G     = spm_cond_units(LastL);
        for i = 1:n
            L(:,i) = G(:,:,i)*P.L(:,i);
        end
              
    % Imaging solution {specified in dipfit.G}
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
            LastL(:,i) = dipfit.G{i}*P.L(:,i);
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
            n = dipfit.Ns;
        catch
            n = m;
        end
        L     = sparse(1:m,1:m,P.L,m,n);
        
        % assume common sources contribute to the last channel
        %------------------------------------------------------------------
        if isfield(dipfit,'common_source')
            L(m,m:n) = L(m,m);
        end

    otherwise
        warndlg('unknown spatial model')
end

% -------------------------------------------------------------------------
function new = transform_points(M,old)
old(:,4) = 1;
new      = old*M';
new      = new(:,1:3);
