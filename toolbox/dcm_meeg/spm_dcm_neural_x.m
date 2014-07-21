function [x] = spm_dcm_neural_x(P,M)
% Returns the fixed point or steady-state of a neural mass DCM
% FORMAT [x] = spm_dcm_neural_x(P,M)
%
% P   - parameter structure
% M   - model structure
%
% x   - steady state solution
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_neural_x.m 6112 2014-07-21 09:39:53Z karl $

switch lower(char(M.f))
    
    % conductance based models
    %----------------------------------------------------------------------
    case lower({'spm_fx_cmm'})
        
        % solve for fixed point
        %------------------------------------------------------------------
        ns    = size(P.A{1},1);     % number of sources (endogenous inputs)
        a     = 2;                  % regulariser
        M.u   = sparse(ns,1);
        dnx   = 0;
        for i = 1:128
            
            % solve under locally linear assumptions
            %--------------------------------------------------------------
            [f,dfdx] = feval(M.f,M.x,M.u,P,M);
            dx       = - dfdx\f;
            
            % regularise
            %--------------------------------------------------------------
            ndx   = norm(dx,Inf);
            if ndx < dnx
                a = a/2;
            end
            dnx    = ndx;
            
            % update and convergence
            %--------------------------------------------------------------
            M.x    = spm_unvec(spm_vec(M.x) + exp(-a)*dx,M.x);
            if dnx < 1e-12, break, end
            
        end
        
    % conductance based models (rank deficient)
    %----------------------------------------------------------------------
    case lower({'spm_fx_mfm'})
        
        % solve for fixed point
        %------------------------------------------------------------------
        ns    = size(P.A{1},1);
        M.u   = sparse(ns,1);
        
        % solve for fixed point using ode113
        %------------------------------------------------------------------
        M.g   = {};
        U.u   = sparse(4,ns);
        U.dt  = 64/1000;
        x     = spm_int_ode(P,M,U);
        x     = spm_unvec(x(end,:),M.x);
        M.x   = x;

        
    % convolution based models
    %----------------------------------------------------------------------
    otherwise
        
end

% solution
%--------------------------------------------------------------------------
x    = M.x;
