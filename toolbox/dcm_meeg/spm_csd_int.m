function [CSD,ERP,csd,mtf,w,pst,x,dP] = spm_csd_int(P,M,U)
% Time frequency response of a neural mass model
% FORMAT [CSD,ERP,csd,mtf,w,pst,x,dP] = spm_csd_int(P,M,U)
%        ERP = spm_csd_int(P,M,U): M.ds = 0
%
% P    - parameters
% M    - neural mass model structure
% M.ds - down-sampling for comutational efficiency (default = 1)
%        if ~M.ds then the ERP is returned
% U    - time-dependent input
%
% ERP  - {E(t,nc}}      - event-related average (sensor space)
% CSD  - {Y(t,w,nc,nc}} - cross-spectral density for nc channels {trials}
%                       - for w frequencies over time t in M.Hz
% csd  - {G(t,w,nc,nc}} - cross spectrum density (before sampling)
% mtf  - {S(t,w,nc,nu}} - transfer functions
% w    - frequencies
% pst  - peristimulus time (sec)
% x    - expectation of hidden (neuronal) states (for last trial)
% dP   - {dP(t,np)}        - parameter fluctuations (plasticity)
%__________________________________________________________________________
%
% This integration routine evaluates the responses of a neural mass model
% to exogenous input - in terms of neuronal states. These are then used as
% expansion point to generate complex cross spectral responses due to
% random neuronal fluctuations. The ensuing spectral (induced) response is
% then convolved (in time) with a window that corresponds to the window of
% a standard wavelet transform. In other words, this routine generates
% predictions of data features based upon a wavelet transform
% characterisation of induced responses.
%
% If M.analysis = 'ERP' then only the ERP is evaluated
%__________________________________________________________________________
% Copyright (C) 2012-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_int.m 6856 2016-08-10 17:55:05Z karl $


% check input - default: one trial (no between-trial effects)
%--------------------------------------------------------------------------
if nargin < 3
    U.dt = 1/256;
    U.u  = sparse(1,M.m);
    U.X  = sparse(1,0);
end


% equations of motion (and OPT for computing CSD vs. ERP)
%--------------------------------------------------------------------------
if ~isfield(M,'g'),   M.g = @(x,u,P,M) x;      end
if ~isfield(M,'h'),   M.h = @(x,u,P,M) 0;      end

M.f = spm_funcheck(M.f);
M.g = spm_funcheck(M.g);
M.h = spm_funcheck(M.h);

% check input function  u = f(t,P,M)
%--------------------------------------------------------------------------
try, ds   = M.ds;   catch, ds    = 1;           end
try, fu   = M.fu;   catch, fu    = 'spm_erp_u'; end
try, ns   = M.ns;   catch, ns    = 128;         end
try, wnum = M.Rft;  catch, wnum  = 8;           end
try, dt   = U.dt;   catch, dt    = 1/256;       end


% within-trial (exogenous) inputs
%==========================================================================
if ~isfield(U,'u')
    u = feval(fu,(1:ns)*dt,P,M)';
else
    u = U.u';
end

% peristimulus time
%--------------------------------------------------------------------------
ns   = size(u,2);
pst  = (1:ns)*dt;
M.dt = dt;

% between-trial (experimental) inputs
%==========================================================================
try
    X = U.X;
    if ~size(X,1)
        X = sparse(1,0);
    end
catch
    X = sparse(1,0);
end

% number of endogenous inputs and hidden states
%==========================================================================
try
    nx = M.n;
catch
    nx = spm_length(M.x);
end
nu    = size(P.A{1},1);
nc    = size(X,1);
dP    = cell(nc,1);
mtf   = cell(nc,1);
csd   = cell(nc,1);
CSD   = cell(nc,1);
ERP   = cell(nc,1);

% cycle over trials or conditions
%--------------------------------------------------------------------------
for c = 1:size(X,1)
    
    % condition-specific parameters
    %----------------------------------------------------------------------
    Q = spm_gen_Q(P,X(c,:));
    
    % initialise hidden states
    %----------------------------------------------------------------------
    x = spm_vec(spm_dcm_neural_x(Q,M));
        
    % remove state (X) and input (Y) dependent parameter from Q
    %----------------------------------------------------------------------
    if isfield(Q,'X'), Q = rmfield(Q,'X'); end
    if isfield(Q,'Y'), Q = rmfield(Q,'Y'); end
    
    
    % get local linear operator LL and delay operator D
    %======================================================================
    if nargout(M.f) >= 3
        [f0,dfdx,D] = M.f(x(:,1),u(:,1),Q,M);
        
    elseif nargout(M.f) == 2
        [f0,dfdx]   = M.f(x(:,1),u(:,1),Q,M);
        D           = 1;
        
    else
        dfdx        = spm_diff(M.f,x(:,1),u(:,1),Q,M,1);
        D           = 1;
    end
    
    % save Delay operator in M to speed computations
    %----------------------------------------------------------------------
    M.D   = D;
    
    % get local linear (Lie) operator L
    %----------------------------------------------------------------------
    p     = max(abs(real(eig(full(dfdx)))));
    N     = ceil(max(1,dt*p*2));
    L     = (spm_expm(dt*D*dfdx/N) - speye(nx,nx))*spm_inv(dfdx);
    
    % cycle over time - expanding around expected states and input
    %======================================================================
    R     = Q;
    Q     = spm_vec(Q);
    dR    = zeros(size(Q));
    dX    = dR;
    dU    = dR;
    
    % clear persistent variables and solve
    %------------------------------------------------------------------
    clear(char(M.h))
    for i = 1:length(pst)
        
        % hidden states
        %------------------------------------------------------------------
        if i > 1, x(:,i) = x(:,i - 1); end
        
        
        % state-dependent parameters (and plasticity)
        %==================================================================
        if isfield(P,'X'), dU  = P.X*u(:,i);             end
        if isfield(P,'Y'), dX  = P.Y*x(:,i);             end
        if isfield(M,'h'), dR  = M.h(x(:,i),u(:,i),P,M); end
        
        % update and save it necessary
        %------------------------------------------------------------------
        dQ = dR + dU + dX;
        R  = spm_unvec(Q + dQ,R);
        
        if nargout > 7
            dP{c}(:,i) = dQ;
        end
        
        % compute complex cross spectral density
        %==================================================================
        if ds
            if ~rem(i + ds - 1,ds)
            
                % add exogenous input and expand around current states
                %----------------------------------------------------------
                M.u     = sparse(nu,1);
                M.x     = spm_unvec(x(:,i),M.x);
                [g,w,s] = spm_csd_mtf(R,M);
                
                % reset exogenous input and expansion point
                %----------------------------------------------------------
                M.x     = spm_unvec(x(:,1),M.x);
                M       = rmfield(M,'u');
            
            end
            
            % place CSD and transfer functions in response
            %--------------------------------------------------------------
            csd{c}(i,:,:,:) = g{:};
            mtf{c}(i,:,:,:) = s{:};
                        
        end
        
        
        % update dx = (expm(dt*J) - I)*inv(J)*f(x,u) = L*f(x,u)
        %==================================================================
        
        % update expected hidden states
        %------------------------------------------------------------------
        for j = 1:N
            x(:,i) = x(:,i) + L*M.f(x(:,i),u(:,i),R,M);
        end
        
        % and get ERP
        %------------------------------------------------------------------
        erp(:,i)    = M.g(x(:,i),u(:,i),R,M);
        
    end
    
    % expected responses (under wavelet transforms)
    %----------------------------------------------------------------------
    ERP{c} = erp';
    if ds
        CSD{c} = spm_morlet_conv(csd{c},w,dt,wnum);
    else
        CSD{c} = ERP{c};
    end
    
end
