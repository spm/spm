function Markov_blankets_and_NESS
% Meta-modelling of Bayes-optimal responses (Newton's method)
% FORMAT Markov_blankets_and_NESS
%
% This demonstration routine deals with the conditional independence in
% this induced by sparse coupling in a random dynamical systems, where the
% sparse coupling is characterised in terms of the system's Jacobian. At
% nonequilibrium steady-state, this places linear constraints on the
% solenoidal flow (denoted by Q) under a Helmholtz decomposition. The
% resulting curvature of the log density and at nonequilibrium steady-state
% encodes conditional independencies; i.e., when the Hessian is zero. What
% follows are a series of notes illustrating the conditions under which
% conditional independence between internal and external states under a
% Markov blanket partition emerges, either asymptotically as the system
% becomes more dissipative - or under a particular constraints on the
% Jacobian. When invoked symbolic maths is used to illustrate an analytic
% solution for a simple the canonical Markov blanket, using a single
% dimensional state for each subset of a Markovian position. Numerical
% analyses are then used to illustrate how this generalises to high
% dimensional systems. Subsequent notes drill down on particular instances
% in the literature.
%
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% backwards compatibility for symbolic maths
%--------------------------------------------------------------------------
if ~exist('str2sym')
    str2sym = @(x)x;
    
end

%  symbolic maths analysis (four dimensional system)
%==========================================================================
% h s a i:  external, sensory, active and internal states
% syms jhh jhs jha jsh jss jsa jas jaa jai jis jia jii K;
% J = [jhh jhs jha 0; jsh jss jsa 0; 0 jas jaa jai; 0 jis jia jii];

% simple Jacobian (J) with block diagonals parameterised by K
%--------------------------------------------------------------------------
syms  jhs jsh jsa jas K jai jia K;
J = [K jhs 0 0; jsh K jsa 0; 0 jas K jai; 0 0 jia K]


% solenoidal flow (R = -Q') that ensures a symmetric Hessian (see below)
%--------------------------------------------------------------------------
n = size(J,1);
I = eye(n,n);
X = kron(I,J) + kron(J,I);
Y = J - transpose(J);
Y = Y(:);
R = inv(X);
R = R*Y;
R = reshape(R,n,n);
R = simplify(R);

% backwards compatibility for symbolic maths
%--------------------------------------------------------------------------
if ~exist('str2sym')
    
    % exact solution
    %----------------------------------------------------------------------
    str2sym = @(x)x;
    H = -inv(I + R)*J;
    
else
    
    % approximation based upon a Neumann series
    %----------------------------------------------------------------------
    H = -(I - R)*J;

end

% Rational function form of Hessian
%--------------------------------------------------------------------------
N      = cell(n,n);                                  % numerator
D      = cell(n,n);                                  % denominator
C      = cell(n,n);                                  % denominator
for ii = 1:n
    for jj = 1:n
        
        % get the expression for this element of the Hessian
        %------------------------------------------------------------------
        F     = factor(H(ii,jj));
        for i = 1:numel(F)
            
            f     = collect(F(i),K);
            str   = char(f);
            k     = strfind(str,'/');
            
            % numerator or denominator of rational function
            %--------------------------------------------------------------
            if numel(k)
                f = inv(f);
            end
            
            % eliminate all but leading terms in K
            %--------------------------------------------------------------
            for d = 1:16
                dF = diff(f,K,d);
                if eval(dF == 0)
                    f = diff(f,K,d - 1);
                    for j = 1:(d - 1)
                        f = int(f,K);
                    end
                    break
                end
            end
            
            if numel(k)
                f = inv(f);
            end
            F(i) = f;
            
        end
        
        % product of factors and display
        %------------------------------------------------------------------
        C{ii,jj} = char(simplify(collect(prod(F),K)))
        
    end
end

C{end,1}

% numerical analysis of dispersion of the Hessian
%==========================================================================
spm_figure('GetWin','Markov Blankets'); clf;
n     = 8;
for i = 1:64
    
    J = [1 1 0 0;
         1 1 1 0;
         0 1 1 1;
         0 0 1 1];
    J = kron(J,ones(n,n));
    J = J.*randn(size(J));
    J = J - eye(size(J))*8;
    G = diag(1 + rand(size(J,1),1));
    
    
    % Hessian and solenoidal flow at NESS based upon Jacobian
    %----------------------------------------------------------------------
    [H,R]     = spm_ness(J,G);
    RR(:,:,i) = R;
    HH(:,:,i) = H;
    JJ(:,:,i) = J;
    GG(:,:,i) = G;
    
end
subplot(2,2,1),imagesc(std(HH,0,3)), axis image, title('Hessian',   'Fontsize',16)
subplot(2,2,2),imagesc(std(RR,0,3)), axis image, title('Solenoidal','Fontsize',16)
subplot(2,2,3),imagesc(std(JJ,0,3)), axis image, title('Jacobian',  'Fontsize',16)
subplot(2,2,4),imagesc(std(GG,0,3)), axis image, title('Diffusion', 'Fontsize',16)


return

% notes for particular forms (based on Biehl et al)
%==========================================================================

% specify equations of motion in terms of a Jacobian (M)
%--------------------------------------------------------------------------
M = -[5 5 2 0;
     -1 1 2 0;
      0 2 5 5;
      0 -2 -1 1]/3

% solve for solenoidal (R) and divergence (G = I) matrices
%--------------------------------------------------------------------------
G = eye(4);
I = eye(4,4);
X = kron(I,M) + kron(M,I);
Y = spm_vec(M*G - G*M');
R = reshape(X\Y,4,4)

fprintf('%i linear constraints\n',rank(X))
fprintf('solenoidal coupling between internal and active states %2.2f\n\n',R(3,4))


% check consistency M*R + R*M' =  M*G - G*M'
%--------------------------------------------------------------------------
fprintf('constraint norm %2.0f\n',norm((M*R + R*M') - (M*G - G*M')))

% precision (inverse covariance) of steady-state density
%--------------------------------------------------------------------------
U = -(R + G)\M;
C = inv(U);

% check flow constraints (Jacobian M)
%==========================================================================
% where flow   f = (R + G)*d log(p(x))/dx   and
% log(p(x))      = -(1/2)*x'*U*x =>
% d log(p(x))/dx = -U*x =>
% df/dx = M      = -(R + G)*U =>
% U              = -(R + G)\M
%--------------------------------------------------------------------------
dfdx = -(R + G)*U

fprintf('coupling between internal and external states %2.2f\n',dfdx(4,1))


% check for Markov blanket (conditional independence : C|b
%==========================================================================
% conditional covariance C|b = C - C*b*inv(b'*C*b)*b'*C'
% b = contrast of blanket states; i.e., b = [0 1 0 0; 0 0 1 0]';
% or
% C|b = C(i,i) - C(i,b)*pinv(C(b,b))*C(b,i): i = [1,4], b = [2,3]
%--------------------------------------------------------------------------
i   = [1,4]; b = [2,3];
C_b = C(i,i) - C(i,b)*inv(C(b,b))*C(b,i)

fprintf('conditional covariance between internal and external states %2.2f\n',C_b(1,2))


return


% NOTES: nontrivial Markov blankets
%==========================================================================

% architecture
%--------------------------------------------------------------------------
m      = 2/3;                             % off diagonal hessian
N      = 4;                               % number of states
G      = eye(N,N);                        % divergence operator

% specify Hessian; i.e., curvature of steady-state surprisal
%--------------------------------------------------------------------------
U      = [1 m 0 0];
U      = toeplitz(U(1:N));

% direct calculation of combinations of solenoidal flow
%--------------------------------------------------------------------------
clc
JJ    = {};
QQ    = {};
K     = spm_perm_mtx((N^2 - N)/2);
for i = 1:size(K,1)
    
    % block matrix form of solenoidal coupling
    %----------------------------------------------------------------------
    Q    = triu(ones(N,N),1);
    j    = find(Q);
    Q(j) = K(i,:);
    Q    = Q - Q';
    
    % evaluate Jacobian (flow constraints)
    %----------------------------------------------------------------------
    J    = -(G - Q')*U;
    
    % nontrivial Jacobians
    %----------------------------------------------------------------------
    if all(sum(abs(J)<1e-6,2))
        JJ{end + 1} = J;
        QQ{end + 1} = Q;
    end
end

return

% symbolic expressions
%==========================================================================
syms q11 q12 q13 q14 q21 q22 q23 q24 q31 q32 q33 q34 q41 q42 q43 q44
syms j11 j12 j13 j14 j21 j22 j23 j24 j31 j32 j33 j34 j41 j42 j43 j44
syms u11 u12 u13 u14 u21 u22 u23 u24 u31 u32 u33 u34 u41 u42 u43 u44


Q  = [-1 q12 0 0; -q12 -1 0 0; 0 0 -1 q34; 0 0 -q34 -1]
U  = [u11 u12 0 0;u12 u22 u23 0;0 u23 j33 j34;0 0 j34 j44]
J  = Q*U


% NOTES: % sufficient conditions for conditional independence
%==========================================================================

% architecture
%--------------------------------------------------------------------------
m   = 2/3;             % off diagonal hessian
N   = 4;               % external, sensory and active and internal states
G   = eye(N,N);        % divergence operator


% direct calculation of combinations of solenoidal flow
%--------------------------------------------------------------------------
JJ  = {};
QQ  = {};
UU  = {};
K   = spm_perm_mtx((N^2 - N)/2);

for k = 1:size(K,1)
    
    % specify Hessian; i.e., curvature of steady-state surprisal
    %----------------------------------------------------------------------
    U     = triu(ones(N,N),1);
    j     = find(U);
    U(j)  = K(k,:);
    U     = U.*rand(N,N);
    U     = U + U' + G + G;
    
    for i = 1:size(K,1)
        
        % block matrix form of solenoidal coupling
        %------------------------------------------------------------------
        Q    = triu(ones(N,N),1);
        j    = find(Q);
        Q(j) = K(i,:);
        Q    = Q.*rand(N,N);
        Q    = Q - Q';
        
        % evaluate Jacobian (flow constraints)
        %------------------------------------------------------------------
        J    = -(G - Q')*U;
        
        % blanket constraints on Jacobian
        %------------------------------------------------------------------
        if all(sum(~J,2)) && all(all(J^N)) && ...
                ~U(1,4) && ~U(1,3) && ~U(2,4) &&...
                U(1,2) &&  U(2,2) &&  U(3,4)
            
            JJ{end + 1} = J;
            QQ{end + 1} = Q;
            UU{end + 1} = U;
        end
    end
    
end

% remove redundant cases
%--------------------------------------------------------------------------
I      = ones(1,numel(JJ));
for i  = 1:numel(JJ)
    iJ = find(spm_vec([QQ{i},UU{i}]));
    for j = 1:numel(JJ)
        vJ  = spm_vec([QQ{j},UU{j}]);
        if i ~= j && all(vJ(iJ))
            I(i) = 0;
        end
    end
end

% display solutions
%--------------------------------------------------------------------------
I = find(I);
disp('U....'), UU{I}
disp('Q....'), QQ{I}
disp('J....'), JJ{I}
numel(I)


return


% solenoidal coupling and complex analysis
%==========================================================================

% create nearly critical Jacobian
%--------------------------------------------------------------------------
s = 1;
N = 8;
while s > 0 || s < -1/4
    J = (randn(N,N) - eye(N,N))/32;
    e = eig(J);
    s = max(real(e));
end

% simple classical system
%--------------------------------------------------------------------------
% J = [-1 8;
%      -4 -1]/128;

% number of states and diffusion tensor
%--------------------------------------------------------------------------
N   = size(J,1);
G   = eye(N,N)/2048;

% create a model and solve for a realised timeseries
%--------------------------------------------------------------------------
M.f = @(x,u,P,M) P*x + u(:);
M.x = randn(N,1)*4;
U   = randn(1024,N);
U   = U*sqrtm(G*2);
x   = spm_int_J(J,M,U);

% complex analysis
%==========================================================================
eJ  = eig(J,'nobalance','vector');
J   = pinv(eJ)*J*eJ;
G   = pinv(eJ)*G*eJ;

% solve for solenoidal (R) and divergence (G = I) matrices
%--------------------------------------------------------------------------
I   = eye(N,N);
X   = kron(I,J) + kron(conj(J),I);
Y   = spm_vec(J*G - G*J');
R   = reshape(X\Y,N,N);

% check constraint
%--------------------------------------------------------------------------
norm((J*R + R*J') - (J*G - G*J'))

% complex form
%--------------------------------------------------------------------------
R = 1j*real(G)*imag(J)/real(J);

% precision (inverse covariance) of steady-state density
%--------------------------------------------------------------------------
H = -(R + G)\J;
C = inv(H);
J = -(R + G)*H;

% complex forms
%--------------------------------------------------------------------------
H = - real(J)/G;                       % precision
C = - real(J)\G;                       % covariance

% decomposition of kinetic energy
%--------------------------------------------------------------------------
D = real(G*H*G'/G);
Q = real(R*H*R'/G);

% complex form
%--------------------------------------------------------------------------
D = diag(-real(J));                          % dissipative
Q = diag(-real(imag(J)*imag(J)/real(J)));    % solenoidal

% show dynamics
%==========================================================================
[d,i] = max(Q);
pst   = 1:size(x,1);

subplot(2,2,1), bar([D,Q])
xlabel('mode'),ylabel('kinetic energy'), axis square, box off
legend({'dissipative','solendoidal'}), legend(gca,'boxoff')
title('Kinetic energy','FontSize',16)

subplot(2,2,2), plot(x*eJ(:,i))
xlabel('real'),ylabel('imaginary'), axis square, box off
title('Complex motion','FontSize',16)

pst   = 1:size(x,1);
subplot(2,2,4), plot(pst,abs(x*eJ(:,i)),'-',pst,angle(x*eJ(:,i)),':')
xlabel('time'),ylabel('state'), axis square, box off
legend({'dissipative','solendoidal'}), legend(gca,'boxoff')
title('Dissipative dynamics','FontSize',16)


return


% NOTES: predicted variance in terms of the eigenvariates
%==========================================================================
[eJ,vJ] = eig(J,'nobalance','vector');
[eC,vC] = eig(cov(x),'vector');
[vJ,i]  = sort(real(vJ),'descend'); eJ = eJ(:,i);
[vC,i]  = sort(real(vC),'descend'); eC = eC(:,i);

dx    = x*J';
dt    = 2;
Hz    = (1:64)'/64/2;
for i = 1:N
    
    X        = x*eC(:,1:i);
    Y        = X*(pinv(X)*x);
    [psd,Hz] = spm_csd(Y,Hz,1/dt,1);
    VC(:,i)  = sum(psd,2);
    
    X        = x*eJ(:,1:i);
    Y        = X*(pinv(X)*x);
    psd      = spm_csd(Y,Hz,1/dt,1);
    VJ(:,i)  = sum(psd,2);
    
    X        = x*eC(:,1:i);
    Y        = X*(pinv(X)*dx);
    [psd,Hz] = spm_csd(Y,Hz,1/dt,1);
    VC(:,i)  = sum(psd,2);
    
    X        = x*eJ(:,1:i);
    Y        = X*(pinv(X)*dx);
    psd      = spm_csd(Y,Hz,1/dt,1);
    VJ(:,i)  = sum(psd,2);
    
end

% show results as a function of frequency
%--------------------------------------------------------------------------
subplot(2,2,1), plot(Hz,VC)
subplot(2,2,2), plot(Hz,VJ)



