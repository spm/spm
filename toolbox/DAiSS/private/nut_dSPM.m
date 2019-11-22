function [weight]=nut_dSPM(Lp,data, flags) %---------------------------------------------------------
% [weight,eta]=nut_dSPM(Lp,data,flags)
% Lp : lead field
% inputs for regularization contant:
% [1] data.Ryy = sample covariance, for data-dependent regularization
% [2] flags.gamma = user defined regularization constant, or 'auto' for
%     leadfield-based regularization

if nargin<3, flags.gamma=[]; end

L = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
G = L*L'; clear L

if isfield(flags,'snr') 
    if ~isfield(data,'C') || isempty(data.C)
        data.C = eye(size(data.Ryy));
    end
    gamma = trace(G)/(trace(data.C)*flags.snr^2)   
    InvG = inv(G+gamma*data.C);
elseif ~isfield(flags,'gamma') || isempty(flags.gamma)
    gamma = 1e0*max(eig(data.Ryy))
    InvG = inv(G+gamma*eye(size(G)));
elseif isnumeric(flags.gamma)
    gamma = flags.gamma   %* max(eig(data.Ryy))
    InvG = inv(G+gamma*eye(size(G)));
else
    % leadfield-based normalization
    doplot=false;
    x = [-20:20];
    gamma=10.^x; %.* max(eig(data.Ryy));
    numg = length(gamma);
    InvG = zeros([size(G) numg]);
    meanuptr = zeros(numg,1);
    warning('off','MATLAB:nearlySingularMatrix')
    for k=1:numg
        InvG(:,:,k) = inv(G+gamma(k)*eye(size(G)));
        uptr = abs(triu(InvG(:,:,k),1));  
        meanuptr(k) = mean(uptr( find(uptr(:)) ));  % mean of matrix without diagonal
    end
    warning('on','MATLAB:nearlySingularMatrix')
    if doplot, figure; plot(x(1:end-2),-diff(diff(log10(meanuptr)))); end
    [dum,idx]=findpeaks(-diff(diff(log10(meanuptr))),'MINPEAKHEIGHT',.1);  % max of second derivation of InvG across gamma magnitudes
    idx = max(idx);
    fprintf('Optimal gamma: 1e%d\n',x(idx));
    InvG = InvG(:,:,idx);
end

InvG2=InvG^2;

Lp1 = squeeze(Lp(:,1,:));
w1 = zeros(size(Lp1));
if size(Lp,2)>1
    Lp2 = squeeze(Lp(:,2,:));
    w2 = zeros(size(Lp2));
end
if size(Lp,2)>2
    Lp3 = squeeze(Lp(:,3,:));
    w3 = zeros(size(Lp3));
end

for i=1:size(Lp,3)
    InvGLp = InvG*Lp1(:,i);
    InvG2Lp = InvG2*Lp1(:,i);
    J = inv(sqrt(trace(Lp1(:,i)'*InvG2Lp)));
    w1(:,i) = InvGLp * J;
    
    if size(Lp,2)>1
        InvGLp = InvG*Lp2(:,i);
        InvG2Lp = InvG2*Lp2(:,i);
        J = inv(sqrt(trace(Lp2(:,i)'*InvG2Lp)));
        w2(:,i) = InvGLp * J;
    end

    if size(Lp,2)>2
        InvGLp = InvG*Lp3(:,i);
        InvG2Lp = InvG2*Lp3(:,i);
        J = inv(sqrt(trace(Lp3(:,i)'*InvG2Lp)));
        w3(:,i) = InvGLp * J;
    end
end

weight(:,1,:) = w1;
if size(Lp,2)>1
    weight(:,2,:) = w2;
end
if size(Lp,2)>2
    weight(:,3,:) = w3;
end
disp('done');
% end
