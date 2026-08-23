% Integration scheme and graphics
%==========================================================================
function [Q,X,V] = spm_biology_solve(x,u,P,T,dt)
% FORMAT [Q,X,V] = spm_biology_solve(x,u,P,T,dt)
% FORMAT [J]         = spm_biology_solve(Q,X,V,P)
%
% x - hidden states
% u - exogenous input
% P - parameter structure
% 
%   P.d - s.d. of random fluctuations
%
% returns:
%
% Q    - history of microstates (states)
% X    - history of microstates (position)
% V    - history of microstates (velocity)
%
% This auxiliary routine integrates a system - or returns the Jacobian for
% specified states. It deals with the special case of within and between
% particle coupling.
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% equations of motion
%--------------------------------------------------------------------------
fx    = @spm_lorenz_n;

% compute Jacobian for pre-evaluated states
%==========================================================================
if isnumeric(x)
    
    % rearrange inputs
    %----------------------------------------------------------------------
    Q     = x;
    X     = u;
    V     = P;
    P     = T;
    clear x
    x.p   = X(:,:,1);
    x.v   = V(:,:,1);
    x.q   = Q(:,:,1);
    J     = zeros(spm_length(x),spm_length(x),size(Q,3));
    
    %  evaluate Jacobian is for the states provided - and return
    %----------------------------------------------------------------------
    for i = 1:size(Q,3)
        x.p = X(:,:,i);
        x.v = V(:,:,i);
        x.q = Q(:,:,i);
        J(:,:,i) = spm_diff(fx,x,0,P,1);
    end
    Q    = J;
    return
end

% set up
%--------------------------------------------------------------------------
d    = 8;                           % diameter for graphics markers
N    = size(x.p,2);                 % number of microsystems
Q    = zeros(3,N,T);                % history of microstates (states)
X    = zeros(2,N,T);                % history of microstates (position)
V    = zeros(2,N,T);                % history of microstates (velocity)
 
% Integrate
%==========================================================================

% partition indices
%--------------------------------------------------------------------------
ee    = P.p == 1;
ss    = P.p == 2;
ii    = P.p == 3;
aa    = P.p == 4;

dn    = 4;
dt    = dt/dn;
for i = 1:T
    
    % update
    %----------------------------------------------------------------------
    xn    = spm_vec(x);
    for j = 1:dn
        
        % flow
        %------------------------------------------------------------------
        f  = fx(spm_unvec(xn,x),u,P);
        xn = xn + f*dt;

    end
    x    = spm_unvec(xn,x);
    
    % adjacency matrix and states
    %----------------------------------------------------------------------
    X(:,:,i) = x.p;
    V(:,:,i) = x.v;
    Q(:,:,i) = x.q;
    
    % graphics
    %----------------------------------------------------------------------
    plot(x.p(1,ee),x.p(2,ee),'.c','MarkerSize',16), hold on
    plot(x.p(1,ss),x.p(2,ss),'.m','MarkerSize',16)
    plot(x.p(1,ii),x.p(2,ii),'.b','MarkerSize',16)
    plot(x.p(1,aa),x.p(2,aa),'.r','MarkerSize',16)

    px = x.p(1,:) + x.q(2,:)/32;
    py = x.p(2,:) + x.q(3,:)/32;

    %%% plot(px(ee),py(ee),'.c','MarkerSize',4)
    plot(px(ss),py(ss),'.m','MarkerSize',4)
    plot(px(ii),py(ii),'.b','MarkerSize',4)
    plot(px(aa),py(aa),'.r','MarkerSize',4)

    xlabel('Position','FontSize',12), ylabel('Position','FontSize',12)
    title('Dynamics','FontSize',16)
    axis([-1 1 -1 1]*d)
    axis square
    hold off

    % save image
    %------------------------------------------------------------------
%     if i < 128
%         M(i) = getframe(gca);
%     end
    drawnow

end

% set ButtonDownFcn
%--------------------------------------------------------------------------
if false
    h   = findobj(gca);
    set(h(1),'Userdata',{M,16})
    set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
    xlabel('Click for Movie','Color','r')
end

return
 

% Equations of motion
%==========================================================================
function f = spm_lorenz_n(x,u,P)
% Equations of motion for coupled Lorenz attractors
% FORMAT f = spm_lorenz_n(x,u,P)
% x - hidden states (5 x N)
%     p: [2×N double]
%     v: [2×N double]
%     q: [3×N double]
% 
% u - exogenous input
% P - parameters
%
% P.e  - endogenous or energy parameter (depth of potential energy well)
%
% P.p(i) -                % partition
% P.r(i) -                % Rayleigh parameter
% P.d(i) -                % diffusion or random fluctuations
% P.k(i) -                % rate constant
%
% f - flow
%__________________________________________________________________________
 

% orders and flow
%--------------------------------------------------------------------------
[n,N] = size(x.p);
f     = x;

% partition indices
%--------------------------------------------------------------------------
ee    = P.p == 1;
ss    = P.p == 2;
ii    = P.p == 3;
aa    = P.p == 4;

% get distances (Euclidean)
%--------------------------------------------------------------------------
X     = zeros(N,N,n);
for i = 1:n
    d        = ones(N,1)*x.p(i,:);
    X(:,:,i) = (d' - d);
end
D     = sqrt(sum(X.^2,3));
D     = D + eye(N,N)*exp(-8);
X     = X./D;

% Electrochemical dynamics
%==========================================================================

% Lorentz dynamics (Prandtl number = 10)
%--------------------------------------------------------------------------
f.q(1,:) = 10*(x.q(2,:) - x.q(1,:));
f.q(2,:) = P.r.*x.q(1,:) - x.q(2,:) - x.q(1,:).*x.q(3,:);
f.q(3,:) = x.q(2,:).*x.q(1,:) - 8/3*x.q(3,:);


% Electrochemical coupling from x.q(2,:) to x.q(1,:)
%--------------------------------------------------------------------------
f.q(1,ee) = f.q(1,ee) + x.q(2,ee)*exp(-D(ee,ee)*8);
f.q(1,ee) = f.q(1,ee) + x.q(2,ss)*exp(-D(ss,ee)*8);
f.q(1,ee) = f.q(1,ee) + x.q(2,aa)*exp(-D(aa,ee)*8);


f.q(1,ss) = f.q(1,ss) + x.q(2,ss)*exp(-D(ss,ss)*8);
f.q(1,ss) = f.q(1,ss) + x.q(2,aa)*exp(-D(aa,ss)*8);
f.q(1,ss) = f.q(1,ss) + x.q(2,ee)*exp(-D(ee,ss)*8);

f.q(1,ii) = f.q(1,ii) + x.q(2,ii)*exp(-D(ii,ii)*8);
f.q(1,ii) = f.q(1,ii) + x.q(2,aa)*exp(-D(aa,ii)*8);
f.q(1,ii) = f.q(1,ii) + x.q(2,ss)*exp(-D(ss,ii)*8);

f.q(1,aa) = f.q(1,aa) + x.q(2,aa)*exp(-D(aa,aa)*8);
f.q(1,aa) = f.q(1,aa) + x.q(2,ii)*exp(-D(ii,aa)*8);
f.q(1,aa) = f.q(1,aa) + x.q(2,ss)*exp(-D(ss,aa)*8);

% electrochemical response to motion of sensory states
%--------------------------------------------------------------------------
f.q(1,ss) = f.q(1,ss) + x.v(1,ss);

% particle-specific rate constants
%--------------------------------------------------------------------------
f.q = times(f.q,P.k);
 
% Newtonian notion
%==========================================================================

% strong repulsive forces: spatial to spatial coupling
%--------------------------------------------------------------------------
F       = zeros(n,N);
E       = 8*X.*exp(-D*2);

F(:,ee) = F(:,ee) + squeeze(sum(E(ee,ee,:),2))';
F(:,ee) = F(:,ee) + squeeze(sum(E(ee,ss,:),2))';
F(:,ee) = F(:,ee) + squeeze(sum(E(ee,aa,:),2))';

F(:,ss) = F(:,ss) + squeeze(sum(E(ss,ss,:),2))';
F(:,ss) = F(:,ss) + squeeze(sum(E(ss,aa,:),2))';
F(:,ss) = F(:,ss) + squeeze(sum(E(ss,ee,:),2))';

F(:,ii) = F(:,ii) + squeeze(sum(E(ii,ii,:),2))';
F(:,ii) = F(:,ii) + squeeze(sum(E(ii,aa,:),2))';
F(:,ii) = F(:,ii) + squeeze(sum(E(ii,ss,:),2))';

F(:,aa) = F(:,aa) + squeeze(sum(E(aa,aa,:),2))';
F(:,aa) = F(:,aa) + squeeze(sum(E(aa,ii,:),2))';
F(:,aa) = F(:,aa) + squeeze(sum(E(aa,ss,:),2))';
 

% and weak electrochemical coupling (Q)
%--------------------------------------------------------------------------
d       = ones(N,1)*x.q(2,:);
Q       = abs(d' - d)/4;
E       = -X.*(exp(-D).*Q);

F(:,ee) = F(:,ee) + squeeze(sum(E(ee,ee,:),2))';
F(:,ee) = F(:,ee) - squeeze(sum(E(ee,ss,:),2))';
F(:,ee) = F(:,ee) - squeeze(sum(E(ee,aa,:),2))';

F(:,ss) = F(:,ss) + squeeze(sum(E(ss,ss,:),2))';
F(:,ss) = F(:,ss) + squeeze(sum(E(ss,aa,:),2))';
F(:,ss) = F(:,ss) + squeeze(sum(E(ss,ee,:),2))';

F(:,ii) = F(:,ii) + squeeze(sum(E(ii,ii,:),2))';
F(:,ii) = F(:,ii) + squeeze(sum(E(ii,aa,:),2))';
F(:,ii) = F(:,ii) + squeeze(sum(E(ii,ss,:),2))';

F(:,aa) = F(:,aa) + squeeze(sum(E(aa,aa,:),2))';
F(:,aa) = F(:,aa) + squeeze(sum(E(aa,ii,:),2))';
F(:,aa) = F(:,aa) - squeeze(sum(E(aa,ss,:),2))';
 

% Newtonian motion
%--------------------------------------------------------------------------
Q    = 0*[0,-1;1,0] + eye(2,2);
f.p  = x.v;

f.v(:,ee)  = F(:,ee) - x.v(:,ee)   - Q*x.p(:,ee)/4;
f.v(:,ss)  = F(:,ss) - x.v(:,ss)   - Q*x.p(:,ss)/8;
f.v(:,ii)  = F(:,ii) - x.v(:,ii)   - Q*x.p(:,ii)/2;
f.v(:,aa)  = F(:,aa) - x.v(:,aa)   - Q*x.p(:,aa)/8;

% plus random fluctuations
%--------------------------------------------------------------------------
f.v  = f.v  + randn(1,N).*P.d;

% vectorised flow
%--------------------------------------------------------------------------
f    = spm_vec(f);





