function spm_mb_merge(cfg)
% Combine tissue maps together

% $Id$

odir     = cfg.odir{1};
onam     = cfg.onam;
res_file = cfg.result{1};
ix       = cfg.ix;

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
            error('Cannot find file %s or %s.',mu_name, fullfile(pth,[nam ext]));
        else
            mu_name = fullfile(pth,[nam ext]);
        end
    end

    % Save the reordered priors for each population
    for p=1:numel(r.sett.gmm)
        matname    = fullfile(odir,['prior_' onam '_' num2str(p) '.mat']);
        mg_ix      = ix(r.sett.gmm(p).mg_ix);
        [mg_ix,si] = sort(mg_ix);
        pr         = r.sett.gmm(p).pr;
        pr{1}      = pr{1}(:,si);    % m
        pr{2}      = pr{2}(:,si);    % b
        pr{3}      = pr{3}(:,:,si);  % W
        pr{4}      = pr{4}(:,si);    % n
        save(matname,'mg_ix','pr');
    end
else
    mu_name = res_file;
end

Nmu = nifti(mu_name);
mu  = single(Nmu.dat(:,:,:,:));
mu1 = MergeMu(mu,ix);
Nmu.dat.fname  = fullfile(odir,['mu_' onam '.nii']);
Nmu.dat.dim(4) = size(mu1,4);
create(Nmu);
Nmu.dat(:,:,:,:) = mu1;




function mu1 = MergeMu(mu,ix)
K  = size(mu,4)+1;
if numel(ix)~=K
    error('Incorrect index dimension.\n The mean image suggests K=%d, whereas the indices assume K=%d.',K, numel(ix));
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
