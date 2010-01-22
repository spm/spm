function [p dp] = spm_LAP_eval(M,qu,qh)
% evaluates precisions for a LAP model
% FORMAT [p dp] = spm_LAP_eval(M,qu,qh)
%
% p.h     - vector of precisions for causal states (v)
% p.g     - vector of precisions for hidden states (v)
%
% dp.h.dx - dp.h/dx
% dp.h.dv - dp.h/dv
% dp.h.dh - dp.h/dh
%
% dp.g.dx - dp.g/dx
% dp.g.dv - dp.g/dv
% dp.g.dg - dp.g/dg
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_LAP_eval.m 3694 2010-01-22 14:16:51Z karl $


% Get states {qu.v{1},qu.x{1}} in hierarchical form (v{i},x{i})
%--------------------------------------------------------------------------
v          = spm_unvec(qu.v{1},{M(1 + 1:end).v});
x          = spm_unvec(qu.x{1},{M(1:end - 1).x});
v{end + 1} = [];
x{end + 1} = [];


% precisions
%==========================================================================
for i = 1:length(M)

    % precision of causal and hidden states
    %----------------------------------------------------------------------
    h{i,1} = feval(M(i).ph,x{i},v{i},qh.h{i},M(i));
    g{i,1} = feval(M(i).pg,x{i},v{i},qh.g{i},M(i));

end

% Concatenate over hierarchical levels
%--------------------------------------------------------------------------
p.h  = spm_cat(h);
p.g  = spm_cat(g);

if nargout < 2, return, end

% gradients
%==========================================================================

% Determine which gradients to evaluate (default: all x,v and h)
%--------------------------------------------------------------------------
try
    method = M(1).E.precision;
catch
    method = 3;
end

switch method

    % gradients w.r.t. h only (no state-dependent noise)
    %----------------------------------------------------------------------
    case(1)

        for i = 1:length(M)

            % precision of causal and hidden states
            %--------------------------------------------------------------
            dhdh{i,i} = spm_diff(M(i).ph,x{i},v{i},qh.h{i},M(i),3);
            dgdg{i,i} = spm_diff(M(i).pg,x{i},v{i},qh.g{i},M(i),3);
           
        end

        % Concatenate over hierarchical levels
        %------------------------------------------------------------------
        dp.h.dh = spm_cat(dhdh);
        dp.g.dg = spm_cat(dgdg);
        
        % number of variables
        %--------------------------------------------------------------------------
        nx      = numel(spm_vec(x));
        nv      = numel(spm_vec(v));
        nh      = size(dp.h.dh,1);
        ng      = size(dp.g.dg,1);

        dp.h.dx = sparse(nh,nx);
        dp.h.dv = sparse(nh,nv);
        dp.g.dx = sparse(ng,nx);
        dp.g.dv = sparse(ng,nv);

        
    % gradients w.r.t. causal states v and h
    %----------------------------------------------------------------------
    case(2)

        for i = 1:length(M)

            % precision of causal states
            %--------------------------------------------------------------
            dhdv{i,i} = spm_diff(M(i).ph,x{i},v{i},qh.h{i},M(i),2);
            dhdh{i,i} = spm_diff(M(i).ph,x{i},v{i},qh.h{i},M(i),3);

            % precision of hidden states
            %--------------------------------------------------------------
            dgdv{i,i} = spm_diff(M(i).pg,x{i},v{i},qh.g{i},M(i),2);
            dgdg{i,i} = spm_diff(M(i).pg,x{i},v{i},qh.g{i},M(i),3);

        end

        % Concatenate over hierarchical levels
        %------------------------------------------------------------------
        dp.h.dv = spm_cat(dhdv);
        dp.h.dh = spm_cat(dhdh);
        dp.g.dv = spm_cat(dgdv);
        dp.g.dg = spm_cat(dgdg);
        
        % number of variables
        %--------------------------------------------------------------------------
        nx      = numel(spm_vec(x));
        nh      = size(dp.h.dh,1);
        ng      = size(dp.g.dg,1);

        dp.h.dx = sparse(nh,nx);
        dp.g.dx = sparse(ng,nx);

        
        
    % gradients w.r.t. hidden and causal states v, and h
    %----------------------------------------------------------------------
    case(3)

        for i = 1:length(M)

            % precision of causal states
            %--------------------------------------------------------------
            dhdx{i,i} = spm_diff(M(i).ph,x{i},v{i},qh.h{i},M(i),1);
            dhdv{i,i} = spm_diff(M(i).ph,x{i},v{i},qh.h{i},M(i),2);
            dhdh{i,i} = spm_diff(M(i).ph,x{i},v{i},qh.h{i},M(i),3);

            % precision of hidden states
            %--------------------------------------------------------------
            dgdx{i,i} = spm_diff(M(i).pg,x{i},v{i},qh.g{i},M(i),1);
            dgdv{i,i} = spm_diff(M(i).pg,x{i},v{i},qh.g{i},M(i),2);
            dgdg{i,i} = spm_diff(M(i).pg,x{i},v{i},qh.g{i},M(i),3);

        end

        % Concatenate over hierarchical levels
        %------------------------------------------------------------------
        dp.h.dx = spm_cat(dhdx);
        dp.h.dv = spm_cat(dhdv);
        dp.h.dh = spm_cat(dhdh);
        dp.g.dx = spm_cat(dgdx);
        dp.g.dv = spm_cat(dgdv);
        dp.g.dg = spm_cat(dgdg);     
        
    otherwise
        
        disp('Unknown method - spm_LAP_eval')
end


