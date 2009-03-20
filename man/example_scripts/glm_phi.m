function [A,fint] = glm_phi (phi,dt,fb)
% Estimate connectivity parameters using GLM/EMA method
% FORMAT [A,fint] = glm_phi (phi,dt,fb)
%
% phi       [N x Nr] matrix of phase time series
%           (N time points, Nr regions)
% dt        sample period
% fb        bandwidth parameter
%
% A         [Nr x Nr] normalised connectivities
% fint      [Nr x 1] intrinsic frequencies

if iscell(phi)
    Nt=length(phi);
    tmp=[];
    for i=1:Nt,
        tmp=[tmp;phi{i}];
    end
    phi=tmp;
end

[N,Nr]=size(phi);

for i=1:Nr,
    ddphi(:,i)=diff(phi(:,i));
end

% Ignore clocking points in any region
rem=[];
for i=1:Nr,
    rem=[rem;find(abs(ddphi(:,i))>5)];
end
keep=[1:N-1];
keep(rem)=[];

dphi=ddphi(keep,:)/(2*pi*dt);
my_phi=phi(keep,:);
            
A=zeros(Nr,Nr);
for i=1:Nr,
    y=dphi(:,i);
    N=length(y);
    X=[];
    regs=[];
    % Assemble design matrix
    for j=1:Nr,
        if ~(i==j)
            X=[X,sin(my_phi(:,i)-my_phi(:,j))];
            regs=[regs,j];
        end
    end
    X=[X,ones(N,1)];
    beta=pinv(X)*y;
    A(i,regs)=-beta(1:Nr-1);
    fint(i)=beta(Nr);
    yfit(:,i)=X*beta;
end
