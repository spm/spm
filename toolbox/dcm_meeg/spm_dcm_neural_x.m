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
% $Id: spm_dcm_neural_x.m 5019 2012-10-26 19:32:57Z karl $
 
 
% solve for fixed point
%--------------------------------------------------------------------------
model = M.f;                        % neural mass model
ns    = size(P.A{1},1);             % number of sources (endogenous inputs)
a     = 2;                          % regulariser
switch lower(model)
    
    % conductance based models
    %----------------------------------------------------------------------
    case lower({'spm_fx_mfm','spm_fx_cmm'})
        
        M.u   = sparse(ns,1);
        nx    = 0;
        for i = 1:128
            
            % solve under locally linear assumptions
            %--------------------------------------------------------------
            [f,dfdx] = feval(M.f,M.x,M.u,P,M);
            dx       = - spm_pinv(dfdx)*f;

            % regularise
            %--------------------------------------------------------------
            ndx   = norm(dx,Inf);
            if ndx < nx
                a = a/2;
            end
            nx    = ndx;
            
            % update and convergence
            %--------------------------------------------------------------
            M.x    = spm_unvec(spm_vec(M.x) + exp(-a)*dx,M.x);
            if nx < 1e-12, break, end
        end
        
        
    % convolution based models
    %----------------------------------------------------------------------
    otherwise

end

% solution
%--------------------------------------------------------------------------
x    = M.x;
