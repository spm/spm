function varargout = fil_io(varargin)
% Function handles for I/O, used by fil_train
% FORMAT io = fil_io
% io - Structure of function handles
%      io.init  - initialise
%      io.block - read in a block of data
%      io.patch - extract a patch from the read in block
%
% FORMAT dat = io.block(dat,x0,y0,z0)
% dat - data structure
% x0  - x indices
% y0  - y indices
% z0  - z indices
%
% FORMAT [X,J,C] = io.patch(dat,x0,y0,z0, r)
% dat  - data structure
% x0   - x indices
% y0   - y indices
% z0   - z indices
% r    - search radius
%
% FORMAT dat = io.init(varargin)
% varargin - arrays of input filenames containing output from
%            fil_push_train_data
% dat      - data structure
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


[varargout{1:nargout}] = spm_subfun(localfunctions,varargin{:});


function dat = init(varargin)
chan = struct('fname','','type','nii','nii',[],'block',cell(1,nargin));
dat  = struct('dw',[],'Mw',[],'dd',0,'x',[],'y',[],'z',[],'chan',chan);
for i=1:nargin
    fname = varargin{i};
    dat.chan(i).fname = fname;
    [~,~,ext] = fileparts(fname(1,:));
    if strcmp(ext,'.mat')
        dat.chan(i).type = 1;
        t0 = load(deblank(fname(1,:)),'dw','Mw');
        for ii=2:size(fname,1)
            t = load(deblank(fname(ii,:)),'dw','Mw');
            if any(t0.dw(:) ~= t.dw(:)) || any(t0.Mw(:) ~= t.Mw(:))
                error('Incompatible dimensions etc (sparse).');
            end
        end
    else
        dat.chan(i).type = 2;
        Nii = nifti(fname);
        dat.chan(i).nii = Nii;
        t0  = struct('dw',size(Nii(1).dat),'Mw',Nii(1).mat);
        for ii=2:numel(Nii)
            t = struct('dw',size(Nii(ii).dat),'Mw',Nii(ii).mat);
            if any(t0.dw(:) ~= t.dw(:)) || any(t0.Mw(:) ~= t.Mw(:))
                error('Incompatible dimensions etc (nii).');
            end
        end
    end
    if isempty(dat.dw)
        dat.dw = t0.dw(1:3);
        dat.Mw = t0.Mw;
    else
        if any(t0.dw(1:3) ~= dat.dw(1:3)) || any(t0.Mw(:) ~= dat.Mw(:))
            error('Incompatible dimensions etc (nifti).');
        end
    end
end


function dat = block(dat,x0,y0,z0)
msk1 = x0>=1 & x0<= dat.dw(1);
msk2 = y0>=1 & y0<= dat.dw(2);
msk3 = z0>=1 & z0<= dat.dw(3);

[x,y,z] = ndgrid(x0(msk1),y0(msk2),z0(msk3));
ind     = x(:)+dat.dw(1)*(y(:)-1+dat.dw(2)*(z(:)-1));
ind1    = find(kron(msk3',kron(msk2',msk1')));
dat.x   = x0;
dat.y   = y0;
dat.z   = z0;
dat.dd  = [numel(x0) numel(y0) numel(z0)];

for m=1:numel(dat.chan)
    chan = dat.chan(m);
    chan.block = cell(size(chan.fname,1),1);
    if chan.type == 1 % Sparse matrix representation
        for n=1:size(chan.fname,1)
            t = load(deblank(chan.fname(n,:)));
            chan.block{n}         = sparse(256,prod(dat.dd));
            chan.block{n}(:,ind1) = t.S(:,ind);
        end
    else              % NIfTI representation
        for n=1:numel(chan.nii)
            fa = chan.nii(n).dat;
            fa = reshape(fa,[prod(size(fa,1:3)) size(fa,4)]);
            chan.block{n}         = zeros([size(fa,2),prod(dat.dd)],'single');
            chan.block{n}(:,ind1) = single(fa(ind,:))';
        end
    end
    dat.chan(m) = chan;
end


function [X,J,C] = patch(dat,x0,y0,z0, r)
if nargin<5, r = 0; end

X  = cell(1,numel(dat.chan));
J  = cell(1,numel(dat.chan));
C  = cell(1,numel(dat.chan));

% Get jitter info
rr  = -ceil(r):ceil(r);
[gx,gy,gz] = ndgrid(rr,rr,rr);
r2  = gx.^2+gy.^2+gz.^2;
msk = r2<=r.^2;
gx  = [gx(msk) gy(msk) gz(msk)];
nj  = size(gx,1);

for m=1:numel(dat.chan)
    cc = cell(nj,1);
    Xc = cell(nj,1);
    for i=1:nj
        [Xc{i},cc{i}] = read_patch(dat,m,x0+gx(i,1),y0+gx(i,2),z0+gx(i,3));
    end
    c = cc{1};
    for i=2:numel(cc)
        c = union(c,cc{i});
    end
    Jc = cell(nj,1);
    for i=1:numel(cc)
        if numel(cc{i})<numel(c)
            xi = Xc{i};
            x  = zeros([size(xi,1),numel(c),size(xi,3)]);
            x(:,index(cc{i},c),:) = xi;
            Xc{i} = x;
        end
        Jc{i} = single(reshape(sum(Xc{i},2),[size(Xc{i},1),size(Xc{i},3)]));
        Xc{i} = single(Xc{i}(:,1:(end-1),:));
    end
    J{m} = cat(2,Jc{:});
    X{m} = cat(3,Xc{:});
    C{m} = c;
end


function [X,c] = read_patch(dat,m,x0,y0,z0)
ind  = indices(dat,x0,y0,z0);
blck = dat.chan(m).block;
S    = cell(numel(blck),1);
for n=1:numel(S)
    S{n} = blck{n}(:,ind);
end

if dat.chan(m).type == 1
    c    = false;
    for n=1:numel(S)
        c    = c | any(S{n},2);
    end
    c = find(c)';
else
    c = 1:size(blck{1},1);
end

X = zeros([numel(ind) numel(c) numel(S)],'single');
for n=1:numel(S)
    X(:,:,n) = full(S{n}(c,:))';
end


function ind = indices(dat,x0,y0,z0)
x = index(x0,dat.x);
y = index(y0,dat.y)-1;
z = index(z0,dat.z)-1;
[x1,y1,z1] = ndgrid(x,y,z);
ind  = x1(:)+dat.dd(1)*(y1(:)+dat.dd(2)*(z1(:)));


function ind = index(A,B)
assert(all(ismember(A,B)));
ind = zeros(size(A));
for i=1:numel(A)
    ind(i) = find(A(i)==B);
end
