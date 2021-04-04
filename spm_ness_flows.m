function spm_ness_flows(P,x,M)
% nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT spm_ness_flows(M,x)
%--------------------------------------------------------------------------
% M   - model specification structure
% Required fields:
%    M.X  - sample points
%    M.W  - (n x n) - precision matrix of random fluctuations
%    M.K  - order of polynomial expansion
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_hd.m 8085 2021-03-21 12:27:26Z karl $


U = spm_ness_U(M,x);                    % get state space and flow

% solve for a trajectory using Laplace form 
%--------------------------------------------------------------------------
for i = 1:size(U.X,1)
    M.X    = U.X(i,:);
    [F,S,Q,L,H,D] = spm_NESS_gen(P,M);
    
    f1(i,:) = 1;
    
end

subplot(3,2,4)
plot3(t(1,:),t(2,:),t(3,:))
title('State-space','Fontsize',16)
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
axis square, box off

subplot(3,2,5)
plot(LF)
title('Potential','Fontsize',16)
xlabel('time'), ylabel('self-information')
axis square xy, box off


return

