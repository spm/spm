function [P] = spm_dcm_fmri_graph_gen(x,v,P)
% Generates adjacency graph for DCM for CSD (fMRI)
% FORMAT [g] = spm_dcm_fmri_graph_gen(x,v,P)
%
% This routine computes the connecivity graph for DCM
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_graph_gen.m 5746 2013-11-14 20:28:50Z karl $


% compute bias for log connectivity using functional space
%==========================================================================

% Distance-based bias on (empirical) prior mean of log connectivity
%--------------------------------------------------------------------------
[n m] = size(v.x);
P.A   = full(P.A);

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
    
    % one-state model of (MoG) connectivity
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







