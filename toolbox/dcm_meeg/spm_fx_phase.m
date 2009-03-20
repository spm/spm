function [f,J] = spm_fx_phase (phi,u,P,M)
% State equation for a phase-coupled oscillator
% FORMAT [f,J] = spm_fx_phase (phi,u,P,M)
%
% phi       state variable
% u         []
% P         model (variable) parameter structure
% M         model (fixed) parameter structure
%
% f         Flow vector, dphi/dt
% J         Jacobian, J(i,j)=df_i/dphi_j
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging
 
% Will Penny
% $Id: spm_fx_phase.m 2908 2009-03-20 14:54:03Z will $

% Sin terms
if isfield(P,'As')
    Nr=size(P.As,1);
    Ns=size(P.As,3);
    
    % negative abs
    Ahs=~(P.As==0);
    As=-abs(P.As).*Ahs;
else
    Ns=0;
end

% Cos terms
if isfield(P,'Ac')
    Nc=size(P.Ac,3);
    Nr=size(P.Ac,1);
    
    % Positive abs
    Ahc=~(P.Ac==0);
    Ac=abs(P.Ac).*Ahc;
else
    Nc=0;
end

% Flow vector
f=zeros(Nr,1);
for i=1:Nr,
    change=0;
    for j=1:Nr,
        phi_diff=phi(i)-phi(j);
        for n=1:Ns,
            change=change+As(i,j,n)*sin(n*phi_diff);
        end
        for n=1:Nc,
            change=change+Ac(i,j,n)*cos(n*phi_diff);
        end
    end
    exc=M.freq+P.df(i)+change;
    f(i)=2*pi*exc;
end

if nargout == 1; return, end

% Jacobian
As=-As;
for i=1:Nr,
    for j=1:Nr,
        jac=0;
        if i==j
            % Diagonal
            for jj=1:Nr,
                for n=1:Ns,
                    jac=jac-n*As(i,jj,n);
                end
            end
        else
            % Off-diagonal
            phi_diff=phi(i)-phi(j);
            for n=1:Ns,
                jac=jac+n*As(i,j,n)*cos(n*phi_diff);
            end
            for n=1:Nc,
                jac=jac+n*Ac(i,j,n)*sin(n*phi_diff);
            end
        end
        J(i,j)=jac;
    end
end
J=J*2*pi;