function out = spm_mb_merge(cfg)
% Combine tissue maps together
% FORMAT out = spm_mb_merge(cfg)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


out      = struct('mu','priors');
odir     = cfg.odir{1};
onam     = cfg.onam;
res_file = cfg.result{1};
ix       = cfg.ix;
bb       = cfg.bb;

[~,~,ext] = fileparts(res_file);
if strcmp(ext,'.mat')
    r = load(res_file);
    if ~isfield(r,'sett'), error('Incorrect format.'); end
    if isfield(r.sett.mu,'create')
        mu_name = r.sett.mu.create.mu;
    else
        mu_name = r.sett.mu.exist.mu;
    end

    if ~exist(mu_name,'file')
        [~,nam,ext] = fileparts(mu_name);
        [pth,~,~]   = fileparts(res_file);
        if ~exist(fullfile(pth,[nam ext]),'file')
            error(sprintf('Cannot find file %s or %s.',mu_name, fullfile(pth,[nam ext])));
        else
            mu_name = fullfile(pth,[nam ext]);
        end
    end

    % Save the reordered priors for each population
    out.priors = cell(numel(r.sett.gmm),1);
    for p=1:numel(r.sett.gmm)
        matname    = fullfile(odir,sprintf('prior_%s_%d.mat',onam,p));
        if max(r.sett.gmm(p).mg_ix) ~= numel(ix)
            error('The indices are incorrectly specified.');
        end
        mg_ix      = ix(r.sett.gmm(p).mg_ix);
        [mg_ix,si] = sort(mg_ix);
        pr         = r.sett.gmm(p).pr;
        pr{1}      = pr{1}(:,si);    % m
        pr{2}      = pr{2}(:,si);    % b
        pr{3}      = pr{3}(:,:,si);  % V
        pr{4}      = pr{4}(:,si);    % n
        save(matname,'mg_ix','pr');
        out.priors{p} = matname;
    end
else
    mu_name = res_file;
end

Nmu       = nifti(mu_name);
[ind,mat] = SubVol(Nmu,bb);
mu        = single(Nmu.dat(ind{:},:));
mu1       = merge_mu(mu,ix);
Nmu.dat.fname = fullfile(odir,['mu_' onam '.nii']);
Nmu.dat.dim   = size(mu1);
Nmu.mat   = mat;
Nmu.mat0  = mat;
create(Nmu);
Nmu.dat(:,:,:,:) = mu1;
out.mu    = {Nmu.dat.fname};
%==========================================================================

%==========================================================================
function mu1 = merge_mu(mu,ix)
K  = size(mu,4)+1;
if numel(ix)~=K
    error(sprintf('Incorrect index dimension.\n The mean image suggests K=%d, whereas the indices assume K=%d.',K, numel(ix)));
end
if any(ix>K) || any(ix<1)
    error('Index out of bounds.');
end

K1  = max(ix);
if numel(unique(ix))~=K1, error('Something wrong with indices.'); end
d   = [size(mu) 1 1];
d   = d(1:3);
mu1 = zeros([d K1-1],'single');
mK  = combine(mu,find(ix==K1));
for k1=1:(K1-1)
    mk  = combine(mu,find(ix==k1));
    mu1(:,:,:,k1) = mk - mK;
end
%==========================================================================

%==========================================================================
function mk = combine(mu,ind)
K = size(mu,4)+1;
d = [size(mu,1) size(mu,2) size(mu,3)];
if numel(ind)==1
    if ind<K
        mk = mu(:,:,:,ind);
    else
        mk = zeros(d,'single');
    end
else
    mx = 0;
    for k=ind
        if k<K
            mx = max(mu(:,:,:,k),mx);
        end
    end
    mk  = 0;
    for k=ind
        if k<K
            mk = mk + exp(mu(:,:,:,k)-mx);
        else
            mk = mk + exp(-mx);
        end
    end
    mk = log(mk)+mx;
end

%==========================================================================

%==========================================================================
function [ind,mat] = SubVol(Nii,bb)
% Information for extracting a subvolume
% FORMAT [ind,mat] = SubVol(Nii,bb)
% V  - mapped file
% bb - bounding box

dm0       = [size(Nii.dat) 1];
dm0       = dm0(1:3);
msk       = ~isfinite(bb(1,:));
bb(1,msk) = 1;
bb(1,:)   = round(max(bb(1,:),1));
msk       = ~isfinite(bb(2,:));
bb(2,msk) = dm0(msk);
bb(2,:)   = round(min(bb(2,:),dm0));
bb        = sort(bb);
ind       = {bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3)};
mat       = Nii.mat*[eye(4,3) [bb(1,:)'-1; 1]];
%==========================================================================

