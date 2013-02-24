% Integration scheme and graphics
%==========================================================================
function [Q,X,V,A,x] = spm_Manifold_solve(x,u,P,T,dt,PLOT)
% FORMAT [Q,X,V,A,x] = spm_Manifold_solve(x,u,P,T,dt,PLOT)
% PLOT = 0 - no grphics
% PLOT = 1 - quick graphics
% PLOT = 2 - quick graphics with trajectory
% PLOT = 3 - pretty graphics


% defaults
%--------------------------------------------------------------------------
try, T;    catch, T    = 512;  end
try, dt;   catch, dt   = 1/64; end
try, PLOT; catch, PLOT = 0;    end
 
% set up
%--------------------------------------------------------------------------
if PLOT, cla, end
d    = 8;                           % diameter for graphics markers
N    = size(x.p,2);                 % number of microsystems
Q    = zeros(3,N,T);                % history of microstates (states)
X    = zeros(2,N,T);                % history of microstates (position)
V    = zeros(2,N,T);                % history of microstates (dusty)
A    = zeros(N,N,T);                % adjacency matrix
 
% Integrate
%--------------------------------------------------------------------------
fx    = @spm_lorenz_n;
dn    = 16;
dt    = dt/dn;
for i = 1:T
    
    % update
    %----------------------------------------------------------------------
    xn    = spm_vec(x);
    for j = 1:dn
        [f a] = fx(spm_unvec(xn,x),u,P);
        xn    = xn + f*dt;
    end
    x       = spm_unvec(xn,x);
    
    % adjacency matrix and states
    %----------------------------------------------------------------------
    A(:,:,i) = a;
    Q(:,:,i) = x.q;
    X(:,:,i) = x.p;
    V(:,:,i) = x.v;
    
    % graphics
    %----------------------------------------------------------------------
    C     = mean(x.q,2);
    if PLOT == 3
        for j = 1:N
            
            % colour
            %--------------------------------------------------------------
            c  = x.q(:,j) - C;
            c  = spm_softmax(c/std(c));
            
            % plot positions
            %--------------------------------------------------------------
            set(gca,'color','k')
            
            px = x.p(1,j) + x.q([1 2 3],j)/32;
            py = x.p(2,j) + x.q([2 3 1],j)/32;
            plot(px,py,'.c','MarkerSize',24,'color',c), hold on

        end
        xlabel('Position','FontSize',12)
        title('Dynamics','FontSize',16)
        axis([-1 1 -1 1]*d)
        axis square
        hold off
        
        % save image
        %------------------------------------------------------------------
        if i < 128
             M(i) = getframe(gca);
        end
        
    elseif PLOT > 0
        
        % plot positions
        %------------------------------------------------------------------
        set(gca,'color','w')
        
        px = ones(3,1)*x.p(1,:) + x.q([1 2 3],:)/32;
        py = ones(3,1)*x.p(2,:) + x.q([2 3 1],:)/32;
        plot(px,py,'.b','MarkerSize',8), hold on
        px = x.p(1,:);
        py = x.p(2,:);
        plot(px,py,'.c','MarkerSize',24)
        axis([-1 1 -1 1]*d)
        
        % plot trajectories
        %------------------------------------------------------------------
        if PLOT == 1, hold off, end
        xlabel('Position','FontSize',12)
        title('Dynamics','FontSize',16)
        axis square
        
    end
    drawnow
    
end

% set ButtonDownFcn
%--------------------------------------------------------------------------
if PLOT == 3
    h   = findobj(gca);
    set(h(1),'Userdata',{M,16})
    set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
    xlabel('Click for Movie','Color','r')
end
 
 
return
 
 
 
 
% Equations of motion
%==========================================================================
 
function [f,A] = spm_lorenz_n(x,u,P)
% Equations of motion for coupled Lorenz attractors
% FORMAT [f,A] = spm_lorenz_n(x,u,P)
% x - hidden states (3 x N)
% u - exogenous input
% P - parameters
%
% P.e  - endogenous or energy parameter (depth of potential energy well)
%
% P.a  - indices of dynamically inert systems (out)
% P.b  - indices of dynamically inert systems (in)
% P.k  - variations in temporal scale
%
% f - flow
% A - weighted adjacency matrix
%__________________________________________________________________________
 

% defaults
%--------------------------------------------------------------------------
if isfield(P,'b'), P.b = find(P.b); else, P.b = []; end
if isfield(P,'a'), P.a = find(P.a); else, P.a = []; end

% orders and flow
%--------------------------------------------------------------------------
[n N] = size(x.p);
f     = x;

% get distances (Euclidean)
%--------------------------------------------------------------------------
D     = 0;
X     = cell(n);
for i = 1:n
    d    = ones(N,1)*x.p(i,:);
    X{i} = (d' - d);
    D    = D + X{i}.^2;
end
 
% enforce short-range coupling (and constraints)
%--------------------------------------------------------------------------
A        = D < 1;                     % coupling in range
A        = A - diag(diag(A));         % remove self conections
A(P.a,:) = 0;                         % preclude outward edges
A(:,P.b) = 0;                         % preclude inward edges
[i j]    = find(A);
k        = find(A);
C        = sparse(i,j,1./sqrt(D(k)),N,N);

% get distances (Dynamic)
%--------------------------------------------------------------------------
d        = ones(N,1)*x.q(3,:);
D        = abs(d' - d);
Q        = sparse(i,j,8*exp(-D(k)*2) - 2,N,N);
 
% Lorentz dynamics
%==========================================================================
 
% State-dependent coupling
%--------------------------------------------------------------------------
xq  = x.q(:,:)*A/8;
xp  = x.q(1,:)*A;
 
% Lorentz dynamics (Prandtl number = 10)
%--------------------------------------------------------------------------
q        = x.q + xq;
 
f.q(1,:) =        10*(q(2,:) - q(1,:));
f.q(2,:) = (32 + xp).*q(1,:) - q(2,:) - q(1,:).*q(3,:);
f.q(3,:) =    q(2,:).*q(1,:) - 8/3*q(3,:);
 
f.q      = ones(3,1)*P.k.*f.q;
 
 
% Newtonian notion
%==========================================================================
 
% strong (repulsive - C) and weak (electrochemical - Q) forces
%--------------------------------------------------------------------------
CQ    = full((C.^2).*(Q - C));
F     = zeros(n,N);
for i = 1:n
    F(i,:)  = sum(X{i}.*CQ);
end
 
% Newtonian motion
%--------------------------------------------------------------------------
f.p  = x.v;
f.v  = F*32 - x.v*8 - x.p;
 
% vectorised flow
%--------------------------------------------------------------------------
f    = [f.p(:); f.v(:); f.q(:)];

% plus random fluctuations
%--------------------------------------------------------------------------
f    = f + randn(size(f));
