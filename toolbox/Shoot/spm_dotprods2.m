function out = spm_dotprods2(job)
% Generate a kernel from dot-products of images
% FORMAT spm_dotprods(job)
% job.images  - Images to use
% job.dotprod - Part of filename for results
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dotprods2.m 5260 2013-02-19 16:23:33Z john $

P      = strvcat(job.images);

N = nifti(P);
n = numel(N);
dm= size(N(1).dat);
dat=cell(1,numel(N));
for i=1:numel(N),
    dat{i} = reshape(N(i).dat,[prod(dm),1]);
end

if isfield(job,'images2'),
    P2 = strvcat(job.images2);
    N2 = nifti(P2);
    n2 = numel(N2);
    dm = size(N2(1).dat);
    dat2=cell(1,numel(N2));
    for i=1:numel(N2),
        dat2{i} = reshape(N2(i).dat,[prod(dm),1]);
    end
    K = zeros(n2,n);
else
    K = zeros(n,n);
end

if isfield(job,'weight') && ~isempty(job.weight),
    Pmsk = strvcat(job.weight);
    Nmsk = nifti(Pmsk);
    msk  = Nmsk.dat;
    dmsk = size(msk);
    if any(dmsk(1:3) ~= dm(1:3)),
        error('Wrong sized weighting image.');
    end
    msk = reshape(msk,[prod(dmsk),1]);
    if numel(dmsk)==3,
        msk1 = msk;
        for i=2:prod(dm(4:end)),
            msk = [msk;msk1];
        end
    end
end

mem = 32*1024*1024;  % Mbytes of RAM to use
bs  = ceil(mem/8/n); % Block size
nd  = prod(dm);
nblock = ceil(prod(dm)/bs);
spm_progress_bar('Init',nblock,...
                 'Generating kernel','Blocks complete');
for k=1:nblock,
    o = bs*(k-1)+(1:bs);
    o = o(o<nd);
    if exist('msk','var'),
        wt  = msk(o);
        tmp = wt>0;
        o   = o(tmp);
        wt  = wt(tmp);
    end
    if ~isempty(o),
        X = zeros(numel(o),numel(dat));
        for i=1:n,
            tmp    = dat{i}(o);
            tmp(~isfinite(tmp)) = 0;
            if exist('wt','var'), tmp = tmp.*wt; end
            X(:,i) = tmp;
        end
        if isfield(job,'images2')
            X2 = zeros(numel(o),numel(dat2));
            for i=1:n2,
                tmp    = dat2{i}(o);
                tmp(~isfinite(tmp)) = 0;
                if exist('wt','var'), tmp = tmp.*wt; end
                X2(:,i) = tmp;
            end
            K    = K + X2'*X;
            clear X2
        else
            K    = K + X'*X;
        end
        clear X
    end
    spm_progress_bar('Set',k);
end
spm_progress_bar('Clear');

if isfield(job,'dotprod')
    [pth,nam,ext] = fileparts(job.dotprod);
    ofname = fullfile(pwd,['dp_' nam '.mat']);
    input  = job;
    typ    = 'images';
    save(ofname,'K','input','typ');
    out.fname = {ofname};
else
    out = K;
end

