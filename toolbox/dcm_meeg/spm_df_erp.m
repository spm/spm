function [df] = spm_df_erp(x,u,P,M,flag)
% gradients wrt states and parameters for a neural mass model of erps
% FORMAT [df] = spm_df_erp(x,u,P,M,flag)
%
% x      - state vector
%   x(:,1) - voltage (spiny stellate cells)
%   x(:,2) - voltage (pyramidal cells) +ve
%   x(:,3) - voltage (pyramidal cells) -ve
%   x(:,4) - current (spiny stellate cells)    depolarizing
%   x(:,5) - current (pyramidal cells)         depolarizing
%   x(:,6) - current (pyramidal cells)         hyperpolarizing
%   x(:,7) - voltage (inhibitory interneurons)
%   x(:,8) - current (inhibitory interneurons) depolarizing
%   x(:,9) - voltage (pyramidal cells)
% u     - current input to the system
% P     - parameter structure
% M     - model structure
% flag  - 'dfdx' or 'dfdP'
% df    - df(t)/dx(t) or df(t)/dP
%
% Prior fixed parameter scaling [Defaults]
%
%  M.pF.E = [32 16 4];           % extrinsic rates (forward, backward, lateral)
%  M.pF.G = [1 4/5 1/4 1/4]*128; % intrinsic rates (g1, g2 g3, g4)
%  M.pF.D = [2 32];              % propogation delays (intrinsic, extrinsic)
%  M.pF.H = [4 32];              % receptor densities (excitatory, inhibitory)
%  M.pF.T = [8 16];              % synaptic constants (excitatory, inhibitory)
%  M.pF.R = [2 1]/3              % parameters of static nonlinearity
%
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_df_erp.m 2322 2008-10-09 14:54:22Z jean $




% get dimensions and configure state variables
%--------------------------------------------------------------------------
n  = length(P.A{1});             % number of sources
x  = reshape(x,n,9);             % neuronal states

% [default] fixed parameters
%--------------------------------------------------------------------------
try
    E  = M.pF.E;                % extrinsic rates (forward, backward, lateral)
    G  = M.pF.G;                % intrinsic rates (g1, g2 g3, g4)
    D  = M.pF.D;                % propogation delays (intrinsic, extrinsic)
    H  = M.pF.H;                % receptor densities (excitatory, inhibitory)
    T  = M.pF.T;                % synaptic constants (excitatory, inhibitory)
    R  = M.pF.R;                % parameters of static nonlinearity
catch
    E = [32 16 4];              % extrinsic rates (forward, backward, lateral)
    G = [1 4/5 1/4 1/4]*128;    % intrinsic rates (g1, g2 g3, g4)
    D = [2 32];                 % propogation delays (intrinsic, extrinsic)
    H = [4 32];                 % receptor densities (excitatory, inhibitory)
    T = [8 16];                 % synaptic constants (excitatory, inhibitory)
    R = [2 1]/3;                % parameters of static nonlinearity
end

% test for free parameters on intrinsic connections
%--------------------------------------------------------------------------
G     = ones(n,1)*G;
try
    G = G.*exp(P.G);
end

% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
A{1}  = exp(P.A{1})*E(1);
A{2}  = exp(P.A{2})*E(2);
A{3}  = exp(P.A{3})*E(3);
C     = exp(P.C);

% intrinsic connectivity and parameters
%--------------------------------------------------------------------------
Te    = T(1)/1000*exp(P.T);           % excitatory time constants
Ti    = T(2)/1000;                    % inhibitory time constants
Hi    = H(2);                         % inhibitory receptor density
He    = H(1)*exp(P.H);                % excitatory receptor density

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R     = R.*exp(P.S);
S     = 1./(1 + exp(-R(1)*(x - R(2)))) - 1./(1 + exp(R(1)*R(2)));
dSdx  = 1./(1 + exp(-R(1)*(x - R(2)))).^2.*(R(1)*exp(-R(1)*(x - R(2))));

% exogenous input
%==========================================================================
U     = C*u(:);



switch flag

    case 'dfdx'


        % Jacobian: J = df(x)/dx
        %===========================================================================
        I  = speye(n,n);
        J  = kron(sparse(9,9),sparse(n,n));

        % changes in voltage with current
        %--------------------------------------------------------------------------
        S  = sparse([7 2 1 3 9],[8 5 4 6 5],1,9,9);  J = J + kron(S,I);
        S  = sparse(9,6,-1,9,9);                     J = J + kron(S,I);

        % synaptic kernel
        %--------------------------------------------------------------------------
        S  = sparse([8 4 5],[8 4 5],1,9,9); J = J - kron(S,diag(2./Te));
        S  = sparse(6,6,1,9,9);             J = J - kron(S,I)*2/Ti;
        S  = sparse([8 4 5],[7 1 2],1,9,9); J = J - kron(S,diag(1./(Te.*Te)));
        S  = sparse(6,3,1,9,9);             J = J - kron(S,I)/(Ti*Ti);

        % Supragranular layer (inhibitory interneurons)
        %--------------------------------------------------------------------------
        E  = (A{2} + A{3})*diag(dSdx(:,9)) + diag(G(:,3).*dSdx(:,9));
        E  = diag(He./Te)*E;
        S  = sparse(8,9,1,9,9); J = J + kron(S,E);

        % Granular layer (spiny stellate cells)
        %--------------------------------------------------------------------------
        E  = (A{1} + A{3})*diag(dSdx(:,9)) + diag(G(:,1).*dSdx(:,9));
        E  = diag(He./Te)*E;
        S  = sparse(4,9,1,9,9); J = J + kron(S,E);

        % Infra-granular layer (pyramidal cells)
        %--------------------------------------------------------------------------
        E  = (A{2} + A{3})*diag(dSdx(:,9));
        E  = diag(He./Te)*E;
        S  = sparse(5,9,1,9,9); J = J + kron(S,E);

        E  = diag(G(:,2).*dSdx(:,1));
        E  = diag(He./Te)*E;
        S  = sparse(5,1,1,9,9); J = J + kron(S,E);

        % Infra-granular layer (pyramidal cells)
        %--------------------------------------------------------------------------
        E  = diag(G(:,4).*dSdx(:,7));
        E  = Hi*E/Ti;
        S  = sparse(6,7,1,9,9); J = J + kron(S,E);

        df = J;
        return


    case 'dfdP'

        % partial derivative w.r.t parameters (for fast EM)
        %==================================================================

        nu = length(u);
        dfdP = sparse(n*9,(2+4*n)*n+(n+1)*nu+2);

        C = kron(~sparse(9,1),speye(n));
        ind = find(C);


        % Synaptic time constants
        %------------------------------------------------------------------
        df = sparse(n,9);
        df(:,8) = x(:,7)./Te.^2-(He.*((A{2}+A{3})*S(:,9)+G(:,3).*S(:,9))-2*x(:,8)-x(:,7)./Te)./Te;
        df(:,4) = x(:,1)./Te.^2-(He.*((A{1}+A{3})*S(:,9)+G(:,1).*S(:,9)+U)-2*x(:,4)-x(:,1)./Te)./Te;
        df(:,5) = x(:,2)./Te.^2-(He.*((A{2}+A{3})*S(:,9)+G(:,2).*S(:,1))-2*x(:,5)-x(:,2)./Te)./Te;
        C(ind) = df';
        dfdP(:,1:n) = C;


        % Synaptic density
        %------------------------------------------------------------------
        df = sparse(n,9);
        df(:,8) = He.*((A{2} + A{3})*S(:,9)+G(:,3).*S(:,9))./Te;
        df(:,4) = He.*((A{1} + A{3})*S(:,9) + G(:,1).*S(:,9) + U)./Te;
        df(:,5) = He.*((A{2} + A{3})*S(:,9) + G(:,2).*S(:,1))./Te;
        C(ind) = df';
        dfdP(:,n+1:2*n) = C;


        % Activation (sigmoid) parameters
        %------------------------------------------------------------------
        dSdR = zeros([size(S),2]);
        dSdR(:,:,1) = R(1).*(-1./(1+exp(-R(1)*(x-R(2)))).^2.*(-x+R(2)).*exp(-R(1)*(x-R(2)))...
            +1./(1+exp(R(1)*R(2))).^2*R(2)*exp(R(1)*R(2)));
        dSdR(:,:,2) = R(2).*(-1./(1+exp(-R(1)*(x-R(2)))).^2.*R(1).*exp(-R(1)*(x-R(2)))...
            +1./(1+exp(R(1)*R(2))).^2*R(1).*exp(R(1)*R(2)));
        df = sparse(n,9);
        df(:,8) = (He.*((A{2} + A{3})*dSdR(:,9,1) + G(:,3).*dSdR(:,9,1)))./Te;
        df(:,4) = (He.*((A{1} + A{3})*dSdR(:,9,1) + G(:,1).*dSdR(:,9,1)))./Te;
        df(:,5) = (He.*((A{2} + A{3})*dSdR(:,9,1) + G(:,2).*dSdR(:,1,1)))./Te;
        df(:,6) = (Hi*G(:,4).*dSdR(:,7,1))/Ti;
        dfdP(:,2*n+1) = df(:);
        df(:,8) = (He.*((A{2} + A{3})*dSdR(:,9,2) + G(:,3).*dSdR(:,9,2)))./Te;
        df(:,4) = (He.*((A{1} + A{3})*dSdR(:,9,2) + G(:,1).*dSdR(:,9,2)))./Te;
        df(:,5) = (He.*((A{2} + A{3})*dSdR(:,9,2) + G(:,2).*dSdR(:,1,2)))./Te;
        df(:,6) = (Hi*G(:,4).*dSdR(:,7,2))/Ti;
        dfdP(:,2*n+2) = df(:);



        % Connectivity parameters
        %------------------------------------------------------------------
        ze = sparse(n,n.^2);
        C = kron(~sparse(1,n),speye(n));
        ind = find(C);
        nS9 = kron((S(:,9))',diag(He./Te));

        % A{1}
        df = sparse(n.*9,n.^2);
        ze(ind) = A{1};
        df(3*n+1:4*n,:) = ze.*nS9;
        dfdP(:,2*n+3:(2+n)*n+2) = df;

        % A{2}
        df = sparse(n.*9,n.^2);
        ze(ind) = A{2};
        df(4*n+1:5*n,:) = ze.*nS9;
        df(7*n+1:8*n,:) = ze.*nS9;
        dfdP(:,(2+n)*n+3:(2+2*n)*n+2) = df;

        % A{3}
        df = sparse(n.*9,n.^2);
        ze(ind) = A{3};
        df(3*n+1:4*n,:) = ze.*nS9;
        df(4*n+1:5*n,:) = ze.*nS9;
        df(7*n+1:8*n,:) = ze.*nS9;
        dfdP(:,(2+2*n)*n+3:(2+3*n)*n+2) = df;

        % Input gain
        C = kron(~sparse(1,nu),speye(n));
        C(C~=0) = exp(P.C)*diag(u);
        df = sparse(n.*9,n*nu);
        df(3*n+1:4*n,:) = C.*kron(ones(1,nu),diag(He./Te));
        dfdP(:,(2+3*n)*n+3:(2+nu+3*n)*n+2) = df;


        % % Delays (the tricky bit)
        % %----------------------------------------------------------------
        % These are worked out outside this function, since they do not
        % enter the vector field (see [David O. et al., 2005]).


        % Input (df/du))
        %------------------------------------------------------------------
        % NB: to get df/dP_u - ie input parameters-, one has to postmultiply by the
        % input function derivatives du/dP_u...
        df = sparse(n*9,nu);
        df(3*n+1:4*n,:) = diag(He./Te)*exp(P.C);
        dfdP(:,(2+nu+4*n)*n+3:end) = df;

        df = dfdP;
        return

end



