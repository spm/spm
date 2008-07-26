function [f,x] = spm_mountaincar_fun(P,G,Q)
% objective function for mountain car problem
% FORMAT [f x] = spm_mountaincar_fun(P)
%
% P = spm_vec(P)
% P.a - 0th order coefficents of force
% P.b - 1st order coefficents of force
% P.c - 2nd order coefficents of force
%
% G   - world model
% Q   - desired equilibrium density
%
% f   - KL divergence between actual and desired equilibrium densities
% x   - cell of grid points supprting 
%
% see:
% Gaussian Processes in Reinforcement Learning
% Carl Edward Rasmussen and Malte Kuss
% Max Planck Institute for Biological Cybernetics
% Spemannstraße 38, 72076 T¨ubingen, Germany
% {carl,malte.kuss}@tuebingen.mpg.de
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mountaincar_fun.m 1961 2008-07-26 09:38:46Z karl $


% place paramters in model
%--------------------------------------------------------------------------
P(end + 1) = G(1).pE.d;
G(1).pE = spm_unvec(P,G(1).pE);

% state-space
%--------------------------------------------------------------------------
n    = 32;
x{1} = linspace(-1,1,n);
x{2} = linspace(-1,1,n);

% get equilirium density
%--------------------------------------------------------------------------
[M0,q0,X,x] = spm_fp(G,x);

% desired equilibrium density
%--------------------------------------------------------------------------
try
    q = Q;
catch
    qE    = [.5; 0];
    qC    = [(1/32)^2 0; 0 (4/32)^2];
    qP    = inv(qC);
    for i = 1:n
        for j = 1:n
            xij    = [x{1}(i); x{2}(j)] - qE;
            q(i,j) = exp(-xij'*qP*xij/2);
        end
    end
end

  
% KL divergence at target: p(x = [.5; 0]|m)
%--------------------------------------------------------------------------
q     = q/sum(q(:));
D     = q(:)./q0(:);
i     = find(D > exp(-16));
f     = q(i)'*log(D(i)) + P'*P/64;


    
