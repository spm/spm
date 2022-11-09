function model = fil_prec(model,sett,p)
% Attach matrices for computing priors
% FORMAT model = fil_prec(model,sett)
% model - The learned model from fil_train
% sett  - Settings
%         Uses sett.matname, sett.nu and sett.v0
%
% Takes a fitted model, and converts to a form that allows the
% distributions of latent variables to be estimated by a neural network
% type formulation.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


model(1).P11 = [];
model(1).P12 = [];

dp = [size(model,1), size(model,2), size(model,3)];
latent_names  = cell(dp(3),1);
[pth,nam,~] = fileparts(sett.matname);
for p3=1:dp(3)
    latent_names{p3} = fullfile(pth,sprintf('%s_latent%.3d.mat', nam, p3));
end

fprintf('Converting: ');
latent3    = cell(1,3);
tmp        = load(latent_names{1},'latent');
latent3{3} = tmp.latent;
for p3=1:size(model,3)
    latent3{1} = latent3{2};
    latent3{2} = latent3{3};
    if p3<dp(3)
        tmp        = load(latent_names{p3+1},'latent');
        latent3{3} = tmp.latent;
    else
        latent3{3} = [];
    end

    for p2=1:size(model,2)
        for p1=1:size(model,1)
            patch   = model(p1,p2,p3);

            Z   = latent3{2}(p1,p2).Z;
            V   = latent3{2}(p1,p2).V;
            Z2  = [GetZ(p1  ,p2  ,latent3{3})
                   GetZ(p1  ,p2  ,latent3{1})
                   GetZ(p1  ,p2+1,latent3{2})
                   GetZ(p1  ,p2-1,latent3{2})
                   GetZ(p1+1,p2  ,latent3{2})
                   GetZ(p1-1,p2  ,latent3{2})];
           V2  = blkdiag(GetV(p1  ,p2  ,latent3{3}),...
                         GetV(p1  ,p2  ,latent3{1}),...
                         GetV(p1  ,p2+1,latent3{2}),...
                         GetV(p1  ,p2-1,latent3{2}),...
                         GetV(p1+1,p2  ,latent3{2}),...
                         GetV(p1-1,p2  ,latent3{2}));

            if isempty(Z2), Z2 = zeros(0,size(Z,2)); end

            if ~isempty(patch.mod)
                [patch.P11,patch.P12] = PrecisionParts(Z,Z2,p,V,V2,sett.nu0,sett.v0);
            end
            model(p1,p2,p3) = patch;
        end
    end
    fprintf('.');
end
fprintf('\n');



function Z = GetZ(p1,p2,latent)
if p1>=1 && p1<=size(latent,1) && ...
   p2>=1 && p2<=size(latent,2) && ...
   ~isempty(latent) && ~isempty(latent(p1,p2).Z)
    Z = latent(p1,p2).Z;
else
    Z = [];
end



function V = GetV(p1,p2,latent)
if p1>=1 && p1<=size(latent,1) && ...
   p2>=1 && p2<=size(latent,2) && ...
   ~isempty(latent) && ~isempty(latent(p1,p2).V)
    V = latent(p1,p2).V;
else
    V = [];
end


function [P11,P12] = PrecisionParts(Z,Z2,p,V,V2,nu0,v0)
% FORMAT [P11,P12,W2,V] = PrecisionParts(Z,Z2,p,V,V2,nu0,v0)
% Z   - K  x N
% Z2  - K2 x N
% V   - K  x K
% V2  - K2 x K2
% nu0 - 1  x 1
% v0  - 1  x 1
% p   - N x 1
%
%_______________________________________________________________________
% Copyright (C) 2019-2020 Wellcome Centre for Human Neuroimaging

K    = size(Z ,1);
K2   = size(Z2,1);
Ns   = sum(p); % Effective number of data points
if ~isfinite(nu0), nu0 = size(Z,1)+size(Z2,1)-0.999; end

% E[Z*Z'], weighted by p
EZZ = [Z *bsxfun(@times,p,Z')+V, Z *bsxfun(@times,p,Z2')
       Z2*bsxfun(@times,p,Z'),   Z2*bsxfun(@times,p,Z2')+V2];

% Wishart posterior
% See https://en.wikipedia.org/wiki/Wishart_distribution
% P ~ W(Psi,nu);
%Psi0 = eye(K+K2)/(nu0*v0)
Psi = inv(EZZ + (nu0*v0)*eye(K+K2));
nu  = Ns+nu0;
P   = Psi*nu; % E[P]
P11 = P(1:K,    1:K  );
P12 = P(1:K, K+(1:K2));
