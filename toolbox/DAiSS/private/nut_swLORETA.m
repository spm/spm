function [weight]=nut_swLORETA(Lp,data,flags)
% [weight,eta]=nut_swLORETA(Lp,data,flags)
% Lp : lead field ( channels X 3 )
% specify either:
% [1] data.Ryy (sample covariance) normally required (used for data-dependent regularization)
%         currently sets gamma = max(eig(data.Ryy))
%         [probably should be set lower for best compromise between stability and blurriness]
% [2] flags.gamma = regularization constant [optional]

warning('BEWARE: not sure if algorithm is faithfully reproduced here, as paper is a little confusing');

L = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
G = L*L';

clear L;

if(isfield(flags,'gamma'))
    gamma = flags.gamma;
else
    gamma = 1e0*max(eig(data.Ryy)); % data-dependent regularization
end

    
InvG = inv(G+gamma*eye(size(G)));

% specific to swLORETA
%%%%% Lsig is equivalent to norm(Lp(:,:,i)) !!!
% Lsig = sparse(size(Lp,3),size(Lp,3));
% for i=1:size(Lp,3)
%     [u,s,v]=svd(Lp(:,:,i),'econ');
%     Lsig(i,i) = s(1,1);
% end
% sqrtSj=(kron(inv(sqrt(Lsig)),eye(3)));  
%
% weight=zeros(size(Lp));
% for i=1:size(Lp,3)
%     weight(:,:,i) = InvG*Lp(:,1:3,i) * inv(sqrtSj((3*i-2):(3*i),(3*i-2):(3*i)));
% end

weight=zeros(size(Lp));
for i=1:size(Lp,3)
    invsqrtnormL = (norm(Lp(:,:,i))).^-.5;
    weight(:,:,i) = InvG*Lp(:,:,i)*invsqrtnormL;
end



% weight = reshape(w,size(Lp));

disp('done');
% end
