function model = fil_prune(model,sett,p)
% Prune the model
% FORMAT model = fil_prune(model,sett,p)
% model - The learned model from fil_train
%
% Take a fitted model, orthogonalise and remove irrelevent latent
% variables.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


if nargin<2, error('Need settings.'); end
if nargin<3, error('Need weights.'); end

[pth,nam,~] = fileparts(sett.matname);

fprintf('Pruning:    ');
for p3=1:size(model,3)
    latent_name = fullfile(pth,sprintf('%s_latent%.3d.mat', nam, p3));
    tmp         = load(latent_name,'latent');
    latent      = tmp.latent;

    for p2=1:size(model,2)
        for p1=1:size(model,1)
            if ~isempty(latent(p1,p2).Z)
                [model(p1,p2,p3),latent(p1,p2)] = orthogonalise(model(p1,p2,p3),latent(p1,p2),p);
            end
        end
    end

    save(latent_name,'latent');
    fprintf('.');
end
fprintf('\n');


function [patch,latent] = orthogonalise(patch,latent,p)
Z       = latent.Z;
V       = latent.V;
mod     = patch.mod;

ZZ      = Z*bsxfun(@times,p,Z');
[~,~,R] = svd(ZZ); % Rotation to diagonalise ZZ
Z       = R'*Z;    % Rotate the matrices.
V       = R'*V*R;
vw      = 0;
for l=1:numel(mod)
    Nvox     = size(mod(l).W,1);
    M        = size(mod(l).W,2);
    K        = size(mod(l).W,3);
    mod(l).W = reshape(reshape(mod(l).W,[Nvox*M,K])*R,[Nvox,M,K]);
    vw       = vw + squeeze(sum(sum(mod(l).W.^2,1),2));
end
nz  = sqrt(vw.*sum(Z.^2,2)/size(Z,2));
ind = nz>0.001;
Z   = Z(ind,:);
V   = V(ind,ind);
for l=1:numel(mod)
%   mod(l).V = mod(l).V(ind,ind);
    mod(l).W = mod(l).W(:,:,ind);
end
patch.mod = mod;
latent.Z  = Z;
latent.V  = V;
