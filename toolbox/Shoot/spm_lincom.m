function out = spm_lincom(job)
% Generate linear combinations of images.
% FORMAT spm_lincom(job)
% job.images   - Images to use
% job.weights  - Matrix of weights
% job.basename - Part of filename for results
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_lincom.m 4573 2011-11-25 23:01:01Z john $

P      = strvcat(job.images);
[pth,nam,ext] = fileparts(job.basename);
ofname = fullfile(pwd,['dp_' nam '.mat']);

N = nifti(P);
n = numel(N);
dm= size(N(1).dat);
dat=cell(1,numel(N));
for i=1:numel(N),
    dat{i} = reshape(N(i).dat,[prod(dm),1]);
end
W = job.weights;
if size(W,1)~=numel(N),
    error('Incompatible dimensions.');
end

out = cell(size(W,2),1);
for i=1:size(W,2)
    [pth,nam,ext]=fileparts(job.basename);
    Po = fullfile(pwd,sprintf('%s%.4d.nii', nam, i));
    out{i} = Po;
    f  = file_array(Po,dm,'float32',352,1,0);
    No = nifti;
    No.mat = N(1).mat;
    No.mat_intent = N(1).mat_intent;
    No.descrip = 'Linear combination';
    No.dat     = f;
    create(No);
    f  = No.dat;
    f  = reshape(f,[prod(size(f)),1]);
    f(end)=0;
    if i==1,
        odat = f;
    else
        odat = [odat f];
    end
end

mem = 128*1024*1024;  % Mbytes of RAM to use
bs  = ceil(mem/8/n); % Block size
nd  = prod(dm);
nblock = ceil(prod(dm)/bs);
spm_progress_bar('Init',100,...
                 'Generating Linear Combs.','Percent complete');
for k=1:nblock,
    o = bs*(k-1)+(1:bs);
    o = o(o<nd);
    if ~isempty(o),
        X = zeros(numel(o),numel(dat));
        for i=1:n,
            tmp    = dat{i}(o);
            tmp(~isfinite(tmp)) = 0;
            X(:,i) = tmp;
        end
        odat(o,:) = X*W;
        clear X
    end
    spm_progress_bar('Set',k/nblock*100);
end
spm_progress_bar('Clear');

