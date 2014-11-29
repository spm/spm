function [E,dE,f,g] = spm_DEM_eval(M,qu,qp)
% evaluates state equations and derivatives for DEM schemes
% FORMAT [E dE f g] = spm_DEM_eval(M,qu,qp)
%
% M  - model structure
% qu - conditional mode of states
%  qu.v{i} - casual states
%  qu.x(i) - hidden states
%  qu.y(i) - response
%  qu.u(i) - input
% qp - conditional density of parameters
%  qp.p{i} - parameter deviates for i-th level
%  qp.u(i) - basis set
%  qp.x(i) - expansion point ( = prior expectation)
%
% E  - generalised errors  (i.e.., y - g(x,v,P); x[1] - f(x,v,P))
%
% dE:
%  dE.du   - de[1:n]/du
%  dE.dy   - de[1:n]/dy[1:n]
%  dE.dc   - de[1:n]/dc[1:d]
%  dE.dp   - de[1:n]/dp
%  dE.dup  - d/dp[de[1:n]/du
%  dE.dpu  - d/du[de[1:n]/dp
%
% where u = x{1:d]; v[1:d]
%
% To accelerate computations one can specify the nature of the model using
% the field:
%
% M(1).E.linear = 0: full        - evaluates 1st and 2nd derivatives
% M(1).E.linear = 1: linear      - equations are linear in x and v
% M(1).E.linear = 2: bilinear    - equations are linear in x, v & x*v
% M(1).E.linear = 3: nonlinear   - equations are linear in x, v, x*v, & x*x
% M(1).E.linear = 4: full linear - evaluates 1st derivatives (for generalised 
%                                  filtering, where parameters change)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_eval.m 6270 2014-11-29 12:04:48Z karl $
 
 
% get dimensions
%==========================================================================
nl    = size(M,2);                       % number of levels
ne    = sum(spm_vec(M.l));               % number of e (errors)
nv    = sum(spm_vec(M.m));               % number of x (causal states)
nx    = sum(spm_vec(M.n));               % number of x (hidden states)
np    = sum(spm_vec(M.p));               % number of p (parameters)


% evaluate functions at each hierarchical level
%==========================================================================
 
% Get states {qu.v{1},qu.x{1}} in hierarchical form (v{i},x{i})
%--------------------------------------------------------------------------
v     = spm_unvec(qu.v{1},{M(1 + 1:end).v});
x     = spm_unvec(qu.x{1},{M(1:end - 1).x});
for i = 1:(nl - 1)
    p      = spm_unvec(spm_vec(M(i).pE) + qp.u{i}*qp.p{i},M(i).pE);
    f{i,1} = feval(M(i).f,x{i},v{i},p);
    g{i,1} = feval(M(i).g,x{i},v{i},p);
end
 
 
% Get Derivatives
%==========================================================================
persistent D
try
    method = M(1).E.linear;
catch
    method = 0;
end

switch method
    
    % get derivatives at each iteration of D-step - full evaluation
    %----------------------------------------------------------------------
    case{0} 
        
        D     = spm_DEM_eval_diff(x,v,qp,M);
 
        % gradients w.r.t. states
        %------------------------------------------------------------------
        dedy  = D.dedy;
        dedc  = D.dedc;
        dfdy  = D.dfdy;
        dfdc  = D.dfdc;
        dgdx  = D.dgdx;
        dgdv  = D.dgdv;
        dfdv  = D.dfdv;
        dfdx  = D.dfdx;
        dgdxp = D.dgdxp;
        dfdxp = D.dfdxp;
        dgdvp = D.dgdvp;
        dfdvp = D.dfdvp;
    
        % gradients w.r.t. parameters
        %------------------------------------------------------------------
        dgdp  = D.dgdp;
        dfdp  = D.dfdp;
 
    % linear: assume equations are linear in x and v
    %----------------------------------------------------------------------      
    case{1}
        
        % get derivatives and store expansion point (states)
        %------------------------------------------------------------------
        if isempty(D)
            
            D     = spm_DEM_eval_diff(x,v,qp,M);
            D.x   = x;
            D.v   = v;
            
            % gradients w.r.t. states
            %--------------------------------------------------------------
            dedy  = D.dedy;
            dedc  = D.dedc;
            dfdy  = D.dfdy;
            dfdc  = D.dfdc;
            dgdx  = D.dgdx;
            dgdv  = D.dgdv;
            dfdv  = D.dfdv;
            dfdx  = D.dfdx;
            
            % gradients w.r.t. parameters (state-dependent)
            %--------------------------------------------------------------
            dgdxp = D.dgdxp;
            dfdxp = D.dfdxp;
            dgdvp = D.dgdvp;
            dfdvp = D.dfdvp;
        
            % gradients w.r.t. parameters
            %--------------------------------------------------------------
            dgdp  = D.dgdp;
            dfdp  = D.dfdp;
            
        % linear expansion for derivatives w.r.t. parameters
        %------------------------------------------------------------------
        else
 
            % gradients w.r.t. states
            %--------------------------------------------------------------
            dedy  = D.dedy;
            dedc  = D.dedc;
            dfdy  = D.dfdy;
            dfdc  = D.dfdc;
            dgdx  = D.dgdx;
            dgdv  = D.dgdv;
            dfdv  = D.dfdv;
            dfdx  = D.dfdx;
            dgdxp = D.dgdxp;
            dfdxp = D.dfdxp;
            dgdvp = D.dgdvp;
            dfdvp = D.dfdvp;
            
            % gradients w.r.t. parameters
            %--------------------------------------------------------------
            dx    = spm_vec(qu.x{1}) - spm_vec(D.x);
            dv    = spm_vec(qu.v{1}) - spm_vec(D.v);
            dgdp  = D.dgdp;
            dfdp  = D.dfdp;
            for p = 1:np
                dgdp(:,p) = D.dgdp(:,p) + D.dgdxp{p}*dx + D.dgdvp{p}*dv;
                if nx
                    dfdp(:,p) = D.dfdp(:,p) + D.dfdxp{p}*dx + D.dfdvp{p}*dv;
                end
            end
 
        end
        
    % bilinear: assume equations are linear in x and v and x*v
    %----------------------------------------------------------------------   
    case{2}
        
        % get derivatives and store expansion point (states)
        %------------------------------------------------------------------
        if isempty(D)
 
            % get high-order derivatives
            %--------------------------------------------------------------
            [Dv D] = spm_diff('spm_DEM_eval_diff',x,v,qp,M,2);
            
            for i = 1:nv, Dv{i} = spm_unvec(Dv{i},D); end
            D.x   = x;
            D.v   = v;
            D.Dv  = Dv;
           
            % gradients w.r.t. states
            %--------------------------------------------------------------
            dedy  = D.dedy;
            dedc  = D.dedc;
            dfdy  = D.dfdy;
            dfdc  = D.dfdc;
            dgdx  = D.dgdx;
            dgdv  = D.dgdv;
            dfdv  = D.dfdv;
            dfdx  = D.dfdx;
            dgdxp = D.dgdxp;
            dfdxp = D.dfdxp;
            dgdvp = D.dgdvp;
            dfdvp = D.dfdvp;
            
            % gradients w.r.t. parameters
            %--------------------------------------------------------------
            dgdp  = D.dgdp;
            dfdp  = D.dfdp;
            
        % linear expansion for derivatives w.r.t. parameters
        %------------------------------------------------------------------
        else
            
            % gradients w.r.t. causes and data
            %--------------------------------------------------------------
            dedy  = D.dedy;
            dedc  = D.dedc;
            dfdy  = D.dfdy;
            dfdc  = D.dfdc;
                        
            % states (relative to expansion point)
            %--------------------------------------------------------------
            dv    = spm_vec(qu.v{1}) - spm_vec(D.v);
            
            % gradients w.r.t. states
            %--------------------------------------------------------------
            dgdx  = D.dgdx;
            dgdv  = D.dgdv;
            dfdx  = D.dfdx;
            dfdv  = D.dfdv;
            for i = 1:nv; dgdx = dgdx + D.Dv{i}.dgdx*dv(i); end
            for i = 1:nv; dgdv = dgdv + D.Dv{i}.dgdv*dv(i); end
            for i = 1:nv; dfdx = dfdx + D.Dv{i}.dfdx*dv(i); end
            for i = 1:nv; dfdv = dfdv + D.Dv{i}.dfdv*dv(i); end
            
            
            % second-order derivatives
            %--------------------------------------------------------------
            dgdxp = D.dgdxp;
            dgdvp = D.dgdvp;
            dfdxp = D.dfdxp;
            dfdvp = D.dfdvp;
            for p = 1:np
                for i = 1:nv; dgdxp{p} = dgdxp{p} + D.Dv{i}.dgdxp{p}*dv(i); end
                for i = 1:nv; dgdvp{p} = dgdvp{p} + D.Dv{i}.dgdvp{p}*dv(i); end
                for i = 1:nv; dfdxp{p} = dfdxp{p} + D.Dv{i}.dfdxp{p}*dv(i); end
                for i = 1:nv; dfdvp{p} = dfdvp{p} + D.Dv{i}.dfdvp{p}*dv(i); end
            end
 
            
            % gradients w.r.t. parameters
            %--------------------------------------------------------------
            dgdp  = D.dgdp;
            dfdp  = D.dfdp;
            for p = 1:np
                Dgdxp     = (D.dgdxp{p} + dgdxp{p})/2;
                Dgdvp     = (D.dgdvp{p} + dgdvp{p})/2;
                Dfdxp     = (D.dfdxp{p} + dfdxp{p})/2;
                Dfdvp     = (D.dfdvp{p} + dfdvp{p})/2;
                dgdp(:,p) = dgdp(:,p) + Dgdvp*dv;
                dfdp(:,p) = dfdp(:,p) + Dfdvp*dv;
            end
 
        end
        
    % nonlinear: assume equations are bilinear in x and v
    %----------------------------------------------------------------------
    case{3}
        
        % get derivatives and store expansion point (states)
        %------------------------------------------------------------------
        if isempty(D)
 
            % get high-order derivatives
            %--------------------------------------------------------------
            [Dx D] = spm_diff('spm_DEM_eval_diff',x,v,qp,M,1,'q');
            [Dv D] = spm_diff('spm_DEM_eval_diff',x,v,qp,M,2,'q');
            
            for i = 1:nx, Dx{i} = spm_unvec(Dx{i},D); end
            for i = 1:nv, Dv{i} = spm_unvec(Dv{i},D); end
            D.x   = x;
            D.v   = v;
            D.Dx  = Dx;
            D.Dv  = Dv;
           
            % gradients w.r.t. states
            %--------------------------------------------------------------
            dedy  = D.dedy;
            dedc  = D.dedc;
            dfdy  = D.dfdy;
            dfdc  = D.dfdc;
            dgdx  = D.dgdx;
            dgdv  = D.dgdv;
            dfdv  = D.dfdv;
            dfdx  = D.dfdx;
            dgdxp = D.dgdxp;
            dfdxp = D.dfdxp;
            dgdvp = D.dgdvp;
            dfdvp = D.dfdvp;
            
            % gradients w.r.t. parameters
            %--------------------------------------------------------------
            dgdp  = D.dgdp;
            dfdp  = D.dfdp;
            
        % linear expansion for derivatives w.r.t. parameters
        %------------------------------------------------------------------
        else
            
            % gradients w.r.t. causes and data
            %--------------------------------------------------------------
            dedy  = D.dedy;
            dedc  = D.dedc;
            dfdy  = D.dfdy;
            dfdc  = D.dfdc;
                        
            % states (relative to expansion point)
            %--------------------------------------------------------------
            dx    = spm_vec(qu.x{1}) - spm_vec(D.x);
            dv    = spm_vec(qu.v{1}) - spm_vec(D.v);
            
            % gradients w.r.t. states
            %--------------------------------------------------------------
            dgdx  = D.dgdx;
            dgdv  = D.dgdv;
            dfdx  = D.dfdx;
            dfdv  = D.dfdv;
            for i = 1:nx; dgdx = dgdx + D.Dx{i}.dgdx*dx(i); end
            for i = 1:nv; dgdx = dgdx + D.Dv{i}.dgdx*dv(i); end
            for i = 1:nx; dgdv = dgdv + D.Dx{i}.dgdv*dx(i); end
            for i = 1:nv; dgdv = dgdv + D.Dv{i}.dgdv*dv(i); end
            for i = 1:nx; dfdx = dfdx + D.Dx{i}.dfdx*dx(i); end
            for i = 1:nv; dfdx = dfdx + D.Dv{i}.dfdx*dv(i); end
            for i = 1:nx; dfdv = dfdv + D.Dx{i}.dfdv*dx(i); end
            for i = 1:nv; dfdv = dfdv + D.Dv{i}.dfdv*dv(i); end
            
            
            % second-order derivatives
            %--------------------------------------------------------------
            dgdxp = D.dgdxp;
            dgdvp = D.dgdvp;
            dfdxp = D.dfdxp;
            dfdvp = D.dfdvp;
            for p = 1:np
                for i = 1:nx; dgdxp{p} = dgdxp{p} + D.Dx{i}.dgdxp{p}*dx(i); end
                for i = 1:nv; dgdxp{p} = dgdxp{p} + D.Dv{i}.dgdxp{p}*dv(i); end
                for i = 1:nx; dgdvp{p} = dgdvp{p} + D.Dx{i}.dgdvp{p}*dx(i); end
                for i = 1:nv; dgdvp{p} = dgdvp{p} + D.Dv{i}.dgdvp{p}*dv(i); end
                for i = 1:nx; dfdxp{p} = dfdxp{p} + D.Dx{i}.dfdxp{p}*dx(i); end
                for i = 1:nv; dfdxp{p} = dfdxp{p} + D.Dv{i}.dfdxp{p}*dv(i); end
                for i = 1:nx; dfdvp{p} = dfdvp{p} + D.Dx{i}.dfdvp{p}*dx(i); end
                for i = 1:nv; dfdvp{p} = dfdvp{p} + D.Dv{i}.dfdvp{p}*dv(i); end
            end
 
            
            % gradients w.r.t. parameters
            %--------------------------------------------------------------
            dgdp  = D.dgdp;
            dfdp  = D.dfdp;
            for p = 1:np
                Dgdxp     = (D.dgdxp{p} + dgdxp{p})/2;
                Dgdvp     = (D.dgdvp{p} + dgdvp{p})/2;
                Dfdxp     = (D.dfdxp{p} + dfdxp{p})/2;
                Dfdvp     = (D.dfdvp{p} + dfdvp{p})/2;
                dgdp(:,p) = dgdp(:,p) + Dgdxp*dx + Dgdvp*dv;
                dfdp(:,p) = dfdp(:,p) + Dfdxp*dx + Dfdvp*dv;
            end
 
        end
        
    % repeated evaluation of first order derivatives (for Laplace scheme)
    %----------------------------------------------------------------------      
    case{4}
        
        % get derivatives and store expansion point (states)
        %------------------------------------------------------------------
        if isempty(D)
            
            D     = spm_DEM_eval_diff(x,v,qp,M);
            D.x   = x;
            D.v   = v;
            
            % gradients w.r.t. states
            %--------------------------------------------------------------
            dedy  = D.dedy;
            dedc  = D.dedc;
            dfdy  = D.dfdy;
            dfdc  = D.dfdc;
            dgdx  = D.dgdx;
            dgdv  = D.dgdv;
            dfdv  = D.dfdv;
            dfdx  = D.dfdx;
            
            % gradients w.r.t. parameters (state-dependent)
            %--------------------------------------------------------------
            dgdxp = D.dgdxp;
            dfdxp = D.dfdxp;
            dgdvp = D.dgdvp;
            dfdvp = D.dfdvp;
        
            % gradients w.r.t. parameters
            %--------------------------------------------------------------
            dgdp  = D.dgdp;
            dfdp  = D.dfdp;
            
        % re-evaluate first-order derivatives
        %------------------------------------------------------------------
        else
 
            % retain second-order gradients
            %--------------------------------------------------------------
            dgdxp = D.dgdxp;
            dfdxp = D.dfdxp;
            dgdvp = D.dgdvp;
            dfdvp = D.dfdvp;
            
            % re-evaluate first-order gradients
            %--------------------------------------------------------------
            D     = spm_DEM_eval_diff(x,v,qp,M,0);
            dedy  = D.dedy;
            dedc  = D.dedc;
            dfdy  = D.dfdy;
            dfdc  = D.dfdc;
            dgdx  = D.dgdx;
            dgdv  = D.dgdv;
            dfdv  = D.dfdv;
            dfdx  = D.dfdx;
            
            % replace second-order gradients
            %--------------------------------------------------------------
            D.dgdxp = dgdxp;
            D.dfdxp = dfdxp;
            D.dgdvp = dgdvp;
            D.dfdvp = dfdvp;
            
            % gradients w.r.t. parameters
            %--------------------------------------------------------------
            dx    = spm_vec(qu.x{1}) - spm_vec(x);
            dv    = spm_vec(qu.v{1}) - spm_vec(v);
            dgdp  = D.dgdp;
            dfdp  = D.dfdp;
            for p = 1:np
                dgdp(:,p) = D.dgdp(:,p) + D.dgdxp{p}*dx + D.dgdvp{p}*dv;
                if nx
                    dfdp(:,p) = D.dfdp(:,p) + D.dfdxp{p}*dx + D.dfdvp{p}*dv;
                end
            end
 
        end
        
    otherwise
        disp('Unknown method')
 
end
 

% order parameters (d = n = 1 for static models)
%--------------------------------------------------------------------------
d       = M(1).E.d + 1;                    % generalisation order of q(v)
n       = M(1).E.n + 1;                    % embedding order     (n >= d)

% Generalised prediction errors and derivatives
%==========================================================================
Ex      = cell(n,1);
Ev      = cell(n,1);
[Ex{:}] = deal(sparse(nx,1));
[Ev{:}] = deal(sparse(ne,1));
 
% prediction error (E) - causes
%--------------------------------------------------------------------------
for i = 1:n
    qu.y{i} = spm_vec(qu.y{i});
end
Ev{1} = [qu.y{1}; qu.v{1}] - [spm_vec(g); qu.u{1}];
for i = 2:n
    Ev{i} = dedy*qu.y{i} + dedc*qu.u{i} ...      % generalised response
          - dgdx*qu.x{i} - dgdv*qu.v{i};         % and prediction
end
 
% prediction error (E) - states
%--------------------------------------------------------------------------
try
    Ex{1} = qu.x{2} - spm_vec(f);
end
for i = 2:n - 1
    Ex{i} = qu.x{i + 1} ...                      % generalised motion
          - dfdx*qu.x{i} - dfdv*qu.v{i};         % and prediction
end
 
% error
%--------------------------------------------------------------------------
E     = spm_vec({Ev,Ex});
 
 
% Kronecker forms of derivatives for generalised motion
%==========================================================================
if nargout < 2, return, end
 
% dE.dp (parameters)
%--------------------------------------------------------------------------
dgdp  = {dgdp};
dfdp  = {dfdp};
for i = 2:n
    dgdp{i,1} = dgdp{1};
    dfdp{i,1} = dfdp{1};
    for p = 1:np
        dgdp{i,1}(:,p) = dgdxp{p}*qu.x{i} + dgdvp{p}*qu.v{i};
        dfdp{i,1}(:,p) = dfdxp{p}*qu.x{i} + dfdvp{p}*qu.v{i};
    end
end
 
% generalised temporal derivatives: dE.du (states)
%--------------------------------------------------------------------------
dedy  = kron(spm_speye(n,n),dedy);
dedc  = kron(spm_speye(n,d),dedc);
dfdy  = kron(spm_speye(n,n),dfdy);
dfdc  = kron(spm_speye(n,d),dfdc);
dgdx  = kron(spm_speye(n,n),dgdx);
dgdv  = kron(spm_speye(n,d),dgdv);
dfdv  = kron(spm_speye(n,d),dfdv);
dfdx  = kron(spm_speye(n,n),dfdx) - kron(spm_speye(n,n,1),speye(nx,nx));
 
% 1st error derivatives (states)
%--------------------------------------------------------------------------
dE.dy =  spm_cat({dedy; dfdy});
dE.dc =  spm_cat({dedc; dfdc});
dE.dp = -spm_cat({dgdp; dfdp});
dE.du = -spm_cat({dgdx, dgdv  ;
                  dfdx, dfdv});
 
              
% bilinear derivatives
%--------------------------------------------------------------------------
for i = 1:np
    dgdxp{i}  = kron(spm_speye(n,n),dgdxp{i});
    dfdxp{i}  = kron(spm_speye(n,n),dfdxp{i});
    dgdvp{i}  = kron(spm_speye(n,d),dgdvp{i});
    dfdvp{i}  = kron(spm_speye(n,d),dfdvp{i});
    dE.dup{i} = -spm_cat({dgdxp{i}, dgdvp{i};
                          dfdxp{i}, dfdvp{i}});
end
if np
    dE.dpu    =  spm_cell_swap(dE.dup);
else
    dE.dpu    =  {};
end
