function [P] = spm_dcm_fmri_graph_gen(x,v,P)
% Generates adjacency graph for spectral DCM for fMRI
% FORMAT [g] = spm_dcm_fmri_graph_gen(x,v,P)
%
% This routine computes the adjacency matrix (A) for spm_fx_fmri
%
% see also: spm_fx_fmri
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_graph_gen.m 5821 2013-12-31 14:26:41Z karl $


% compute bias for log connectivity using functional space
%==========================================================================

% spectral deomposition
%--------------------------------------------------------------------------
P.A   = full(P.A);

if isfield(P,'modes')
    
    % outer product
    %======================================================================
    n     = length(P.A);                          % number of nodes
    m     = length(v);                            % number of modes
    s     = [-exp(-v); (zeros(n - m,1) - 1)];
    P.A   = P.modes*diag(s)*P.modes';
    P.A   = full(P.A + diag(log(-2*diag(P.A)) - diag(P.A)));
    P     = rmfield(P,'modes');

    return
    
end


if isnumeric(v)
    
    % static modes (explicit negativity constraints on self excitation)
    %======================================================================
    % P.A   = v'*v;
    
    % dynamical modes (implicit negative definite constraints)
    %======================================================================
    P.A   = logm(v'*v + eye(size(v,2),size(v,2))*exp(-16))/16;
    P.A   = P.A + diag(log(-2*diag(P.A)) - diag(P.A));

    return
    
end

% Distance-based bias on (empirical) prior mean of log connectivity
%--------------------------------------------------------------------------
[n m]          = size(v.x);    
if size(P.A,3) == 1 && numel(v.a) == 1
    
    % one-state model of (MoG) connectivity
    %======================================================================
    for i = 1:m
        for j = (i + 1):m
            
            % Euclidean distance
            %--------------------------------------------------------------
            P.A(i,j) =  ...
                exp(v.a - sum((v.x(:,i) - v.x(:,j)).^2)/2)/4 - ... % excitatory
                exp(v.a - sum((v.x(:,i) - v.x(:,j)).^2)/8)/8;      % inhibitory
            P.A(j,i) = P.A(i,j);
            
        end
    end
    
elseif size(P.A,3) == 1 && numel(v.a) == 2
    
    % one-state model of centre-surround connectivity
    %======================================================================
    for i = 1:m
        for j = (i + 1):m
            
            % Euclidean distance
            %--------------------------------------------------------------
            D        = exp(-sum((v.x(:,i) - v.x(:,j)).^2)/2);
            P.A(i,j) = exp(v.a(1))*D/16 + v.a(2);
            P.A(j,i) = P.A(i,j);
            
        end
    end
    
elseif size(P.A,3) == 2
    
    % assume two-state model of log connectivity
    %======================================================================
    for i = 1:m
        for j = (i + 1):m
            
            % Euclidean distance
            %--------------------------------------------------------------
            P.A(i,j,1) = v.a - sum((v.x(:,i) - v.x(:,j)).^2)/2;
            P.A(j,i,1) = P.A(i,j,1);
            
            
            % hierarchical distance
            %--------------------------------------------------------------
            P.A(i,j,2) = (sqrt(sum(v.x(:,i).^2)) - sqrt(sum(v.x(:,j)).^2))/2;
            P.A(j,i,2) = -P.A(i,j,2);
            
        end
    end
    
end







