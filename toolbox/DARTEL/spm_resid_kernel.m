function spm_resid_kernel(job)
% Generate a kernel from dot-products of images
% FORMAT spm_resid_kernel(job)
% job.images  - Images to use
% job.dotprod - Part of filename for results
%_______________________________________________________________________
% Copyright (C) 2007 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$

P      = strvcat(job.images);
[pth,nam,ext] = fileparts(job.dotprod);
ofname = fullfile(pwd,['dp_' nam '.mat']);

N = nifti(P);
n = numel(N);
for i=1:numel(N),
    dat{i} = reshape(N(i).dat,[prod(size(N(i).dat)) 1]);
end
Phi = zeros(n,n);

mem = 32*1024*1024;  % Mbytes of RAM to use
bs  = ceil(mem/8/n); % Block size
nd  = prod(size(dat{1}));
nblock = ceil(prod(size(dat{1}))/bs);
spm_progress_bar('Init',nblock,...
                 'Generating kernel','Blocks complete');
for k=1:nblock,
    o = bs*(k-1)+(1:bs);
    o = o(o<nd);
    X = zeros(numel(o),numel(dat));

    for i=1:n,
        X(:,i) = dat{i}(o);
    end
    Phi    = Phi + X'*X;
    spm_progress_bar('Set',k);
end
spm_progress_bar('Clear');
input = job;
typ   = 'images';
save(ofname,'Phi','input','typ');

