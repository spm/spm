function varargout = spm_mb_shape(varargin)
%__________________________________________________________________________
%
% Functions for shape model related.
%
% FORMAT psi0          = spm_mb_shape('Affine',d,Mat)
% FORMAT B             = spm_mb_shape('AffineBases',code)
% FORMAT psi           = spm_mb_shape('Compose',psi1,psi0)
% FORMAT id            = spm_mb_shape('Identity',d)
% FORMAT dat           = spm_mb_shape('InitDef',dat,sett.ms)
% FORMAT l             = spm_mb_shape('LSE',mu,ax)
% FORMAT sett          = spm_mb_shape('MuValOutsideFOV',mu,sett);
% FORMAT a1            = spm_mb_shape('Pull1',a0,psi,r)
% FORMAT [f1,w1]       = spm_mb_shape('Push1',f,psi,d,r)
% FORMAT sd            = spm_mb_shape('SampDens',Mmu,Mn)
% FORMAT varargout     = spm_mb_shape('Shoot',v0,kernel,args)
% FORMAT mu1           = spm_mb_shape('ShrinkTemplate',mu,oMmu,sett)
% FORMAT P             = spm_mb_shape('Softmax',mu,ax)
% FORMAT E             = spm_mb_shape('TemplateEnergy',mu,sett)
% FORMAT mun           = spm_mb_shape('TemplateK1',mun)
% FORMAT dat           = spm_mb_shape('UpdateAffines',dat,mu,sett)
% FORMAT [mu,dat]      = spm_mb_shape('UpdateMean',dat, mu, sett)
% FORMAT dat           = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett)
% FORMAT [mu,dat]      = spm_mb_shape('UpdateSimpleMean',dat, mu, sett)
% FORMAT dat           = spm_mb_shape('UpdateVelocities',dat,mu,sett)
% FORMAT dat           = spm_mb_shape('UpdateWarps',dat,sett)
% FORMAT [dat,mu]      = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu)
% FORMAT sz            = spm_mb_shape('ZoomSettings', v_settings, mu, n)
% FORMAT [P,datn]      = spm_mb_io('GetClasses',datn,mu,sett)
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% $Id$

if nargin == 0
    help spm_mb_shape
    error('Not enough argument. Type ''help spm_mb_shape'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch id
    case 'Affine'
        [varargout{1:nargout}] = Affine(varargin{:});
    case 'AffineBases'
        [varargout{1:nargout}] = AffineBases(varargin{:});
    case 'Compose'
        [varargout{1:nargout}] = Compose(varargin{:});
    case 'Identity'
        [varargout{1:nargout}] = Identity(varargin{:});
    case 'InitDef'
        [varargout{1:nargout}] = InitDef(varargin{:});
    case 'LSE'
        [varargout{1:nargout}] = LSE(varargin{:});
    case 'MuValOutsideFOV'
        [varargout{1:nargout}] = MuValOutsideFOV(varargin{:});
    case 'Pull1'
        [varargout{1:nargout}] = Pull1(varargin{:});
    case 'Push1'
        [varargout{1:nargout}] = Push1(varargin{:});
    case 'SampDens'
        [varargout{1:nargout}] = SampDens(varargin{:});
    case 'Shoot'
        [varargout{1:nargout}] = Shoot(varargin{:});
    case 'ShrinkTemplate'
        [varargout{1:nargout}] = ShrinkTemplate(varargin{:});
    case 'Softmax'
        [varargout{1:nargout}] = Softmax(varargin{:});
    case 'TemplateEnergy'
        [varargout{1:nargout}] = TemplateEnergy(varargin{:});
    case 'TemplateK1'
        [varargout{1:nargout}] = TemplateK1(varargin{:});
    case 'UpdateAffines'
        [varargout{1:nargout}] = UpdateAffines(varargin{:});
    case 'UpdateMean'
        [varargout{1:nargout}] = UpdateMean(varargin{:});
    case 'UpdateSimpleAffines'
        [varargout{1:nargout}] = UpdateSimpleAffines(varargin{:});
    case 'UpdateSimpleMean'
        [varargout{1:nargout}] = UpdateSimpleMean(varargin{:});
    case 'UpdateVelocities'
        [varargout{1:nargout}] = UpdateVelocities(varargin{:});
    case 'UpdateWarps'
        [varargout{1:nargout}] = UpdateWarps(varargin{:});
    case 'ZoomVolumes'
        [varargout{1:nargout}] = ZoomVolumes(varargin{:});
    case 'ZoomSettings'
        [varargout{1:nargout}] = ZoomSettings(varargin{:});
    case 'GetClasses'
        [varargout{1:nargout}] = GetClasses(varargin{:});
    otherwise
        help spm_mb_shape
        error('Unknown function %s. Type ''help spm_mb_shape'' for help.', id)
end
end
%==========================================================================

%==========================================================================
function psi0 = Affine(d,Mat,samp)
% Generate an affine deformation field
if nargin<3 || numel(samp)~=3
    samp = [1 1 1];
end
d     = [d 1 1];
d     = d(1:3);
id    = Identity(d,samp);
%psi0 = vec2vol(vol2vec(id)*Mat(1:3,1:3)' + Mat(1:3,4)',size(id));
psi0  = vec2vol(bsxfun(@plus,vol2vec(id)*Mat(1:3,1:3)',Mat(1:3,4)'),size(id));
if d(3) == 1, psi0(:,:,:,3) = 1; end
end
%==========================================================================

%==========================================================================
function B = AffineBases(code)
g     = regexpi(code,'(?<code>\w*)\((?<dim>\d*)\)','names');
g.dim = str2num(g.dim);
if numel(g.dim)~=1 || (g.dim ~=0 && g.dim~=2 && g.dim~=3)
    error('Can not use size');
end
if g.dim==0
    B        = zeros(4,4,0);
elseif g.dim==2
    switch g.code
    case 'T'
        % T(2) - translations
        B        = zeros(4,4,2);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
    case 'SO'
        % SO(2) - Special Orthogonal Group - rotations
        B        = zeros(4,4,1);
        B(1,2,1) =  1;
        B(2,1,1) = -1;
    case 'SE'
        % SE(2) - Special Euclidean Group - rigid
        B        = zeros(4,4,3);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
        B(1,2,3) =  1;
        B(2,1,3) = -1;
    otherwise
        error('Unknown group.');
    end
elseif g.dim==3
    switch g.code
    case 'T'
        % T(3)
        B        = zeros(4,4,3);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
        B(3,4,3) =  1;
    case 'SO'
        % SO(3) - Special Orthogonal
        B        = zeros(4,4,3);
        B(1,2,1) =  1;
        B(2,1,1) = -1;
        B(1,3,2) =  1;
        B(3,1,2) = -1;
        B(2,3,3) =  1;
        B(3,2,3) = -1;
    case 'SE'
        % SE(3) - Special Euclidean
        B        = zeros(4,4,6);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
        B(3,4,3) =  1;
        B(1,2,4) =  1;
        B(2,1,4) = -1;
        B(1,3,5) =  1;
        B(3,1,5) = -1;
        B(2,3,6) =  1;
        B(3,2,6) = -1;
    otherwise
        error('Unknown group.');
    end
end
end
%==========================================================================

%==========================================================================
function psi = Compose(psi1,psi0)
% Function composition
% \psi \gets \psi^1 \circ \psi^0
if isempty(psi1)
    psi = psi0;
else
    psi = spm_diffeo('comp',psi1,psi0);
end
if size(psi,3) == 1, psi(:,:,:,3) = 1; end % 2D
end
%==========================================================================

%==========================================================================
function id = Identity(d,samp)
% Identity transform of dimensions d, but with subsampling (samp)
if nargin<2 || numel(samp)~=3
    samp = [1 1 1];
end
d  = [d 1 1];
d  = d(1:3);
ds = floor((d(:)'-1)./samp(:)')+1;
id = zeros([ds(:)',3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:samp(1):d(1)),...
                                               single(1:samp(2):d(2)),...
                                               single(1:samp(3):d(3)));
end
%==========================================================================

%==========================================================================
function dat = InitDef(dat,ms)
% Possibly break up this function and put parts of it in the io file.

% Parse function settings
d       = ms.d;
Mmu     = ms.Mmu;

v    = zeros([d,3],'single');
psi1 = Identity(d);
for n=1:numel(dat)
    if ~isnumeric(dat(n).v) || ~isnumeric(dat(n).psi)
        fa       = file_array('placeholder.nii',[d(1:3) 1 3],'float32',0);
        nii      = nifti;
        nii.dat  = fa;
        nii.mat  = Mmu;
        nii.mat0 = Mmu;
    end
    if isnumeric(dat(n).v)
        dat(n).v   = v;
    else
        if ischar(dat(n).v)
            nii.dat.fname = dat(n).v;
        else
            nii.dat.fname = dat(n).v.dat.fname;
        end
        nii.descrip = 'Velocity';
        create(nii);
        nii.dat(:,:,:,:) = v;
        dat(n).v         = nii;
    end

    if isnumeric(dat(n).psi)
        dat(n).psi   = psi1;
    else
        if ischar(dat(n).psi)
            nii.dat.fname = dat(n).psi;
        else
            nii.dat.fname = dat(n).psi.dat.fname;
        end
        nii.descrip = 'Deformation (WIP)';
        create(nii);
        nii.dat(:,:,:,:) = psi1;
        dat(n).psi       = nii;
    end
end
end
%==========================================================================

%==========================================================================
function l = LSE(mu,ax)
% log-sum-exp function

if nargin<2, ax = 4; end
mx = max(mu,[],ax);
%l = log(exp(-mx) + sum(exp(mu - mx),ax)) + mx;
l  = log(exp(-mx) + sum(exp(bsxfun(@minus,mu,mx)),ax)) + mx;
end
%==========================================================================

%==========================================================================
function a1 = Pull1(a0,psi,r)
% Resample an image or set of images
% FORMAT a1 = Pull1(a0,psi,r)
%
% a0  - Input image(s)
% psi - Deformation
% r   - subsampling density in each dimension (default: [1 1 1])
%
% a1  - Output image(s)
%

if nargin<3, r=[1 1 1]; end

if isempty(a0)
    a1 = a0;
elseif isempty(psi)
    a1 = a0;
else
    if r==1
        a1 = spm_diffeo('pull',a0,psi);
        return
    end
    d  = [size(a0) 1 1];
    if d(3)>1, zrange = range(r(3)); else zrange = 0; end
    if d(2)>1, yrange = range(r(2)); else yrange = 0; end
    if d(1)>1, xrange = range(r(1)); else xrange = 0; end
    id = Identity(size(psi));
    a1 = zeros([size(psi,1),size(psi,2),size(psi,3),size(a0,4)],'single');
    for l=1:d(4)
        tmp = single(0);
        al  = single(a0(:,:,:,l));
        for dz=zrange
            for dy=yrange
                for dx=xrange
                    ids  = id  + cat(4,dx,dy,dz);
                    psi1 = spm_diffeo('pull',psi-id,    ids)+ids;
                    as   = spm_diffeo('pull',al,psi1);
                   %ids  = id  - cat(4,dx,dy,dz);
                    tmp  = tmp + spm_diffeo('push',as,  ids);
                end
            end
        end
        a1(:,:,:,l) = tmp/(numel(zrange)*numel(yrange)*numel(xrange));
    end
end
end
%==========================================================================

%==========================================================================
function [f1,w1] = Push1(f,psi,d,r)
% Push an image (or set of images) accorging to a spatial transform
% FORMAT [f1,w1] = Push1(f,psi,d,r)
%
% f   - Image (3D or 4D)
% psi - Spatial transform
% d   - dimensions of output (default: size of f)
% r   - subsampling density in each dimension (default: [1 1 1])
%
% f1  - "Pushed" image
%

if nargin<4, r = [1 1 1]; end
if nargin<3, d = [size(f,1) size(f,2) size(f,3)]; end
if numel(r)==1, r = repmat(r,[1 3]); end

%msk    = isfinite(f);
%f(msk) = 0;
if ~isempty(psi)
    if all(r==1)
        if nargout==1
            f1      = spm_diffeo('push',single(f),psi,d);
        else
            [f1,w1] = spm_diffeo('push',single(f),psi,d);
        end
        return
    end

    if d(3)>1, zrange = range(r(3)); else zrange = 0; end
    if d(2)>1, yrange = range(r(2)); else yrange = 0; end
    if d(1)>1, xrange = range(r(1)); else xrange = 0; end

    id    = Identity(size(psi));
    f1    = single(0);
    w1    = single(0);
    for dz=zrange
        for dy=yrange
            for dx=xrange
               %ids       = id + cat(4,dx,dy,dz);
                ids       = bsxfun(@plus, id, cat(4,dx,dy,dz));
                psi1      = spm_diffeo('pullc',psi-id,    ids)+ids;
                fs        = spm_diffeo('pull',single(f), ids);
               %fs=single(f);
                if nargout==1
                    fs        = spm_diffeo('push',fs,        psi1,d);
                    f1        = f1  + fs;
                else
                    [fs,ws]   = spm_diffeo('push',fs,        psi1,d);
                    f1        = f1  + fs;
                    w1        = w1  + ws;
                end
            end
        end
    end
    scale = 1/(numel(zrange)*numel(yrange)*numel(xrange));
    f1    = f1*scale;
    w1    = w1*scale;
else
    msk      = isfinite(f);
    f1       = f;
    f1(~msk) = 0;
    w1       = single(all(msk,4));
end
end
%==========================================================================

%==========================================================================
function r = range(n)
r = (-floor((n-1)/2):ceil((n-1)/2))/n;
end
%==========================================================================

%==========================================================================
function sd = SampDens(Mmu,Mn)
vx_mu = sqrt(sum(Mmu(1:3,1:3).^2,1));
vx_f  = sqrt(sum( Mn(1:3,1:3).^2,1));
sd    = max(round(2.0*vx_f./vx_mu),1);
end
%==========================================================================

%==========================================================================
function mu = ShrinkTemplate(mu,oMmu,sett,smo_wt)

% Parse function settings
if nargin<4, smo_wt = 1; end
d     = sett.ms.d;
Mmu   = sett.ms.Mmu;
d0    = [size(mu,1) size(mu,2) size(mu,3)];
Mzoom = Mmu\oMmu;
if any(d0~=d) || norm(Mzoom-eye(4))>1e-4
   %y      = reshape(reshape(Identity(d0),[prod(d0),3])*Mzoom(1:3,1:3)'...
   %                 +Mzoom(1:3,4)',[d0 3]);
    y      = reshape(bsxfun(@plus,...
             reshape(Identity(d0),[prod(d0),3])*Mzoom(1:3,1:3)',...
             Mzoom(1:3,4)'),[d0 3]);
   %mu     = exp(mu-LSE(mu,4));
    mu     = exp(bsxfun(@minus,mu,LSE(mu,4)));
    [mu,c] = Push1(mu,y,d,1);
    e      = eps('single');
   %mu     = mu./max(c,e);
    mu     = bsxfun(@rdivide,mu,max(c,e));
   %mu     = log(max(mu,e))-log(max(1-sum(mu,4),e));
    mu     = bsxfun(@minus,log(max(mu,e)),log(max(1-sum(mu,4),e)));
end
smo_wt = min(max(smo_wt,0),1);
if smo_wt~=0
    smo = single([0.25 0.5 0.25]);
    if d(3)>1
        smo = reshape(kron(kron(smo,smo),smo),[3 3 3]);
        del = zeros([3 3 3],'single');
        del(2,2,2) = 1;
        smo = smo*smo_wt+del*(1-smo_wt);
    else
        smo = smo'*smo;
        del = zeros([3 3],'single');
        del(2,2) = 1;
        smo = smo*smo_wt+del*(1-smo_wt);
    end
    for k=1:size(mu,4)
        mu(:,:,:,k) = convn(mu(:,:,:,k),smo,'same');
    end
end
end
%==========================================================================

%==========================================================================
function P = Softmax(mu,ax)
% safe softmax function

if nargin<2, ax = 4; end
mx  = max(mu,[],ax);
%E  = exp(mu-mx);
E   = exp(bsxfun(@minus,mu,mx));
den = sum(E,ax)+exp(-mx);
%P  = E./den;
P   = bsxfun(@rdivide,E,den);
end
%==========================================================================

%==========================================================================
function E = TemplateEnergy(mu,mu_settings)
if ~isempty(mu_settings)
    spm_field('bound',1);
    m0 = sum(mu,4)/(size(mu,4)+1);
   %mu = mu - m0;
    mu = bsxfun(@minus,mu,m0);
    g  = spm_field('vel2mom', mu, mu_settings);
    E  = 0.5*mu(:)'*g(:);
    E = E + 0.5*sum(sum(sum(m0.*spm_field('vel2mom', m0, mu_settings))));
else
    E  = 0;
end
end
%==========================================================================

%==========================================================================
% TemplateK1()
function mun = TemplateK1(mun,ax)
if nargin<2, ax = 4; end
mx   = max(max(mun,[],ax),0);
%lse = mx + log(sum(exp(mun - mx),ax) + exp(-mx));
lse  = bsxfun(@plus,mx,log(sum(exp(bsxfun(@minus,mun,mx)),ax) + exp(-mx)));
%mun = cat(ax,mun - lse, -lse);
mun  = cat(ax,bsxfun(@minus,mun,lse), -lse);
end
%==========================================================================

%==========================================================================
function dat = UpdateAffines(dat,mu,sett)

% Parse function settings
B         = sett.B;
groupwise = isa(sett.mu,'struct') && isfield(sett.mu,'create');
nw        = GetNumWork(sett,max(27,sett.K*5+17));

% Update the affine parameters
if ~isempty(B)
    if nw > 1 && numel(dat) > 1 % PARFOR
        parfor(n=1:numel(dat),nw) dat(n) = UpdateAffinesSub(dat(n),mu,sett); end
    else % FOR
        for n=1:numel(dat), dat(n) = UpdateAffinesSub(dat(n),mu,sett); end
    end

    if groupwise
        % Zero-mean the affine parameters
        mq = sum(cat(2,dat(:).q),2)/numel(dat);
        for n=1:numel(dat)
            dat(n).q = dat(n).q - mq;
        end
    end
end
end
%==========================================================================

%==========================================================================
function datn = UpdateAffinesSub(datn,mu,sett)
% This could be made more efficient.

% Parse function settings
accel = sett.accel;
B     = sett.B;
d     = sett.ms.d;
Mmu   = sett.ms.Mmu;
scal  = 0.9;

q    = double(datn.q);
Mn   = datn.Mat;
samp = datn.samp;
[Mr,dM3] = spm_dexpm(q,B);
dM   = zeros(12,size(B,3));
for m=1:size(B,3)
    tmp     = Mmu\dM3(:,:,m)*Mn;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

df   = datn.dm;
psi0 = Affine(df,Mmu\Mr*Mn,samp);
ds   = [size(psi0,1),size(psi0,2),size(psi0,3)];
psi1 = spm_mb_io('GetData',datn.psi);
if ~isempty(psi1)
    J    = spm_diffeo('jacobian',psi1);
    J    = reshape(Pull1(reshape(J,[d 3*3]),psi0),[ds 3 3]);
    psi  = Compose(psi1,psi0);
else
    J    = [];
    psi  = psi0;
end
clear psi0  psi1

mu1  = Pull1(mu,psi);
[f,datn] = GetClasses(datn,mu1,sett);
M    = size(mu,4);
G    = zeros([ds M 3],'single');
for m=1:M
    [~,Gm{1},Gm{2},Gm{3}] = spm_diffeo('bsplins',mu(:,:,:,m),psi,[1 1 1  0 0 0]);
    for i1=1:3
        if ~isempty(J)
            tmp = single(0);
            for j1=1:3
                tmp = tmp + J(:,:,:,j1,i1).*Gm{j1};
            end
        else
            tmp = Gm{j1};
        end
        tmp(~isfinite(tmp)) = 0;
        G(:,:,:,m,i1) = tmp;
    end
    clear Gm
end
clear J mu

msk       = all(isfinite(f),4) & all(isfinite(mu1),4);
mu1(~isfinite(mu1)) = 0;
a         = Mask(f - Softmax(mu1,4),msk);
[H,g]     = AffineHessian(mu1,G,a,single(msk),accel,samp);
g         = double(dM'*g);
H         = dM'*H*dM;
H         = H + eye(numel(q))*(norm(H)*1e-5 + 0.01);
q         = q + scal*(H\g);
datn.q    = q;
end
%==========================================================================

%==========================================================================
function [H,g] = AffineHessian(mu,G,a,w,accel,samp)
if nargin<6, samp = [1 1 1]; end
d  = [size(mu,1),size(mu,2),size(mu,3)];
I  = Horder(3);
H  = zeros(12,12);
g  = zeros(12, 1);
[x{1:4}] = ndgrid(((1:d(1))-1)*samp(1)+1,((1:d(2))-1)*samp(2)+1,1,1);
for i=1:d(3)
    x{3} = x{3}*0+(i-1)*samp(3)+1;
   %gv   = reshape(sum(a(:,:,i,:).*G(:,:,i,:,:),4),[d(1:2) 1 3]);
    gv   = reshape(sum(bsxfun(@times,a(:,:,i,:),G(:,:,i,:,:)),4),[d(1:2) 1 3]);
   %Hv   = w(:,:,i).*VelocityHessian(mu(:,:,i,:),G(:,:,i,:,:),accel);
    Hv   = bsxfun(@times,w(:,:,i),VelocityHessian(mu(:,:,i,:),G(:,:,i,:,:),accel));
    for i1=1:12
        k1g   = rem(i1-1,3)+1;
        k1x   = floor((i1-1)/3)+1;
        g(i1) = g(i1) + sum(sum(sum(x{k1x}.*gv(:,:,:,k1g))));
        for i2=1:12
            k2g      = rem(i2-1,3)+1;
            k2x      = floor((i2-1)/3)+1;
            H(i1,i2) = H(i1,i2) + sum(sum(sum(x{k1x}.*Hv(:,:,:,I(k1g,k2g)).*x{k2x})));
        end
    end
end
end
%==========================================================================

%==========================================================================
function [mu,dat] = UpdateMean(dat, mu, sett)

% Parse function settings
accel       = sett.accel;
mu_settings = sett.ms.mu_settings;
s_settings  = [2 2];
nw          = GetNumWork(sett,4*sett.K+4);

spm_field('bound',1);
m0 = sum(mu,4)/(size(mu,4)+1);
%g = spm_field('vel2mom', mu-m0, mu_settings);
g  = spm_field('vel2mom', bsxfun(@minus,mu,m0), mu_settings);
w  = zeros(sett.ms.d,'single');
if nw > 1 && numel(dat) > 1 % PARFOR
    parfor(n=1:numel(dat),nw)
        [gn,wn,dat(n)] = UpdateMeanSub(dat(n),mu,sett);
        g              = g + gn;
        w              = w + wn;
    end
else
    for n=1:numel(dat) % FOR
        [gn,wn,dat(n)] = UpdateMeanSub(dat(n),mu,sett);
        g              = g + gn;
        w              = w + wn;
    end
end
clear gn wn
H  = AppearanceHessian(mu,accel,w);
% Note that spm_field could be re-written to make these updates
% converge more effectively.
mu = mu - spm_field(H, g, [mu_settings s_settings]);
end
%==========================================================================

%==========================================================================
function [g,w,datn] = UpdateMeanSub(datn,mu,sett)

% Parse function settings
B        = sett.B;
d        = sett.ms.d;
Mmu      = sett.ms.Mmu;
df       = datn.dm;
q        = double(datn.q);
Mn       = datn.Mat;
samp     = datn.samp;

psi      = Compose(spm_mb_io('GetData',datn.psi),Affine(df, Mmu\spm_dexpm(q,B)*Mn,samp));
mu       = Pull1(mu,psi);
[f,datn] = GetClasses(datn,mu,sett);
[g,w]    = Push1(Softmax(mu,4) - f,psi,d,1);
end
%==========================================================================

%==========================================================================
function H = AppearanceHessian(mu,accel,w)
M  = size(mu,4);
d  = [size(mu,1) size(mu,2) size(mu,3)];
if accel>0, s  = Softmax(mu,4); end
Ab = 0.5*(eye(M)-1/(M+1)); % See Bohning's paper
I  = Horder(M);
H  = zeros([d (M*(M+1))/2],'single');
for m1=1:M
    for m2=m1:M
        if accel==0
            tmp = Ab(m1,m2)*ones(d,'single');
        else
            if m2~=m1
                tmp = accel*(-s(:,:,:,m1).*s(:,:,:,m2))           + (1-accel)*Ab(m1,m2);
            else
                tmp = accel*(max(s(:,:,:,m1).*(1-s(:,:,:,m1)),0)) + (1-accel)*Ab(m1,m2);
            end
        end
        if nargin>=3
            H(:,:,:,I(m1,m2)) = tmp.*w;
        else
            H(:,:,:,I(m1,m2)) = tmp;
        end
    end
end
end
%==========================================================================

%==========================================================================
function dat = UpdateSimpleAffines(dat,mu,sett)

% Parse function settings
accel     = sett.accel;
B         = sett.B;
groupwise = isa(sett.mu,'struct') && isfield(sett.mu,'create');

% Update the affine parameters
G  = spm_diffeo('grad',mu);
H0 = VelocityHessian(mu,G,accel);

if ~isempty(B)
    for n=1:numel(dat)
        dat(n) = UpdateSimpleAffinesSub(dat(n),mu,G,H0,sett);
    end

    if groupwise
        % Zero-mean the affine parameters
        mq = sum(cat(2,dat(:).q),2)/numel(dat);
        for n=1:numel(dat)
            dat(n).q = dat(n).q - mq;
        end
    end
end
end
%==========================================================================

%==========================================================================
function datn = UpdateSimpleAffinesSub(datn,mu,G,H0,sett)

% Parse function settings
B    = sett.B;
d    = sett.ms.d;
Mmu  = sett.ms.Mmu;
scal = 0.9;

df   = datn.dm;
q    = double(datn.q);
Mn   = datn.Mat;
samp = datn.samp;
[Mr,dM3] = spm_dexpm(q,B);
dM   = zeros(12,size(B,3));
for m=1:size(B,3)
    tmp     = (Mr*Mmu)\dM3(:,:,m)*Mmu;
    dM(:,m) = reshape(tmp(1:3,:),12,1);
end

psi      = Affine(df,Mmu\Mr*Mn,samp);
mu1      = Pull1(mu,psi);

[f,datn] = GetClasses(datn,mu1,sett);
[a,w]    = Push1(f - Softmax(mu1,4),psi,d,1);
clear mu1 psi f

[H,g]     = SimpleAffineHessian(mu,G,H0,a,w);
g         = double(dM'*g);
H         = dM'*H*dM;
H         = H + eye(numel(q))*(norm(H)*1e-5 + 0.1);
q         = q + scal*(H\g);
datn.q    = q;
end
%==========================================================================

%==========================================================================
function [H,g] = SimpleAffineHessian(mu,G,H0,a,w)
samp = [1 1 1];
d  = [size(mu,1),size(mu,2),size(mu,3)];
I  = Horder(3);
H  = zeros(12,12);
g  = zeros(12, 1);
[x{1:4}] = ndgrid(1:d(1),1:d(2),1,1);
for i=1:d(3)
    x{3} = x{3}*0+i;
   %gv   = reshape(sum(a(:,:,i,:).*G(:,:,i,:,:),4),[d(1:2) 1 3]);
    gv   = reshape(sum(bsxfun(@times,a(:,:,i,:),G(:,:,i,:,:)),4),[d(1:2) 1 3]);

   %Hv   = w(:,:,i).*H0(:,:,i,:);
    Hv   = bsxfun(@times,w(:,:,i),H0(:,:,i,:));
    for i1=1:12
        k1g   = rem(i1-1,3)+1;
        k1x   = floor((i1-1)/3)+1;
        g(i1) = g(i1) + sum(sum(sum(x{k1x}.*gv(:,:,:,k1g))));
        for i2=1:12
            k2g      = rem(i2-1,3)+1;
            k2x      = floor((i2-1)/3)+1;
            H(i1,i2) = H(i1,i2) + sum(sum(sum(x{k1x}.*Hv(:,:,:,I(k1g,k2g)).*x{k2x})));
        end
    end
end
end
%==========================================================================

%==========================================================================
function [mu,dat] = UpdateSimpleMean(dat, mu, sett)

% Parse function settings
accel       = sett.accel;
mu_settings = sett.ms.mu_settings;
s_settings  = [2 2];
nw          = GetNumWork(sett,4*sett.K+4);

w  = zeros(sett.ms.d,'single');
gf = zeros(size(mu),'single');
if nw > 1 && numel(dat) > 1 % PARFOR
    parfor(n=1:numel(dat),nw)
        [gn,wn,dat(n)] = UpdateSimpleMeanSub(dat(n),mu,sett);
        gf             = gf + gn;
        w              = w  + wn;
    end
else
    for n=1:numel(dat) % FOR
        [gn,wn,dat(n)] = UpdateSimpleMeanSub(dat(n),mu,sett);
        gf             = gf + gn;
        w              = w  + wn;
    end
end
clear gn wn
for it=1:ceil(4+2*log2(numel(dat)))
    H  = AppearanceHessian(mu,accel,w);
    g  = w.*Softmax(mu,4) - gf;
    spm_field('bound',1);
    m0 = sum(mu,4)/(size(mu,4)+1);
    g  = g  + spm_field('vel2mom', mu-m0, mu_settings);
    % Note that spm_field could be re-written to make the updates slightly
    % more effective.
    mu = mu - spm_field(H, g, [mu_settings s_settings]);
end
end
%==========================================================================

%==========================================================================
% UpdateSimpleMeanSub()
function [g,w,datn] = UpdateSimpleMeanSub(datn,mu,sett)

% Parse function settings
B     = sett.B;
d     = sett.ms.d;
Mmu   = sett.ms.Mmu;
df    = datn.dm;
q     = double(datn.q);
Mn    = datn.Mat;
samp  = datn.samp;

psi      = Compose(spm_mb_io('GetData',datn.psi),Affine(df,Mmu\spm_dexpm(q,B)*Mn,samp));
mu       = Pull1(mu,psi);
[f,datn] = GetClasses(datn,mu,sett);
[g,w]    = Push1(f,psi,d,1);
end
%==========================================================================

%==========================================================================
function dat = UpdateVelocities(dat,mu,sett)

% Parse function settings
accel = sett.accel;
nw    = GetNumWork(sett,4*sett.K+4*3);

G  = spm_diffeo('grad',mu);
H0 = VelocityHessian(mu,G,accel);
if size(G,3) == 1
    % Data is 2D -> add some regularisation
    H0(:,:,:,3) = H0(:,:,:,3) + mean(reshape(H0(:,:,:,[1 2]),[],1));
end
if nw > 1 && numel(dat) > 1 % PARFOR
    parfor(n=1:numel(dat),nw) dat(n) = UpdateVelocitiesSub(dat(n),mu,G,H0,sett); end
else % FOR
    for n=1:numel(dat), dat(n) = UpdateVelocitiesSub(dat(n),mu,G,H0,sett); end
end
end
%==========================================================================

%==========================================================================
function datn = UpdateVelocitiesSub(datn,mu,G,H0,sett)

% Parse function settings
B          = sett.B;
d          = sett.ms.d;
Mmu        = sett.ms.Mmu;
s_settings = [2 2];
scal       = 0.9;
v_settings = sett.ms.v_settings;

v         = spm_mb_io('GetData',datn.v);
q         = datn.q;
Mn        = datn.Mat;
samp      = datn.samp;
Mr        = spm_dexpm(q,B);
Mat       = Mmu\Mr*Mn;
df        = datn.dm;
psi       = Compose(spm_mb_io('GetData',datn.psi),Affine(df,Mat,samp));
mu        = Pull1(mu,psi);
[f,datn]  = GetClasses(datn,mu,sett);
[a,w]     = Push1(f - Softmax(mu,4),psi,d,1);
clear psi f mu

%g        = reshape(sum(a.*G,4),[d 3]);
g         = reshape(sum(bsxfun(@times,a,G),4),[d 3]);
%H        = w.*H0;
H         = bsxfun(@times,w,H0);
clear a w

u0        = spm_diffeo('vel2mom', v, v_settings);                          % Initial momentum
datn.E(2) = 0.5*sum(u0(:).*v(:));                                          % Prior term
v         = v - scal*spm_diffeo('fmg',H, g + u0, [v_settings s_settings]); % Gauss-Newton update

if d(3)==1, v(:,:,:,3) = 0; end % If 2D
if v_settings(1)==0             % Mean displacement should be 0
    avg = mean(mean(mean(v,1),2),3);
    v   = v - avg;
end
datn.v = spm_mb_io('SetData',datn.v,v);
end
%==========================================================================

%==========================================================================
function [P,datn] = GetClasses(datn,mun,sett)

if false
% Put a lower level on how confident the model is when template
% probabilities gecome extremely small. This may improve the
% apparent convergence because it tries less hard to make the
% probabilities in the template approach 0 or 1.
t    = log(eps('single')*128);
lse  = LSE(mun,4);
%mun = max(mun-lse,t)+min(lse,-t);
mun  = bsxfun(@plus,max(bsxfun(@minus,mun,lse),t),min(lse,-t));
end


if isfield(datn.model,'cat')
    % Categorical model
    P  = spm_mb_io('GetData',datn.model.cat.f);
    sk = datn.samp;
    P  = P(1:sk:end,1:sk:end,1:sk:end,:);
    if nargout > 1
        % Compute subject-specific categorical cross-entropy loss between
        % segmentation and template
        msk       = all(isfinite(P),4) & all(isfinite(mun),4);
        tmp       = sum(P.*mun,4) - LSE(mun,4);
        datn.E(1) = -sum(tmp(msk(:)));
    end
else
    % Update appearance model
    [datn,P] = spm_mb_appearance('Update',datn,mun,sett);
end
end
%==========================================================================

%==========================================================================
function H = VelocityHessian(mu,G,accel)
d  = [size(mu,1),size(mu,2),size(mu,3)];
M  = size(mu,4);
if accel>0, s  = Softmax(mu,4); end
Ab = 0.5*(eye(M)-1/(M+1)); % See Bohning's paper
H1 = zeros(d,'single');
H2 = H1;
H3 = H1;
H4 = H1;
H5 = H1;
H6 = H1;
for m1=1:M
    Gm11 = G(:,:,:,m1,1);
    Gm12 = G(:,:,:,m1,2);
    Gm13 = G(:,:,:,m1,3);
    if accel==0
        tmp = Ab(m1,m1);
    else
        sm1 = s(:,:,:,m1);
        tmp = (max(sm1.*(1-sm1),0))*accel + (1-accel)*Ab(m1,m1);
    end
    H1 = H1 + tmp.*Gm11.*Gm11;
    H2 = H2 + tmp.*Gm12.*Gm12;
    H3 = H3 + tmp.*Gm13.*Gm13;
    H4 = H4 + tmp.*Gm11.*Gm12;
    H5 = H5 + tmp.*Gm11.*Gm13;
    H6 = H6 + tmp.*Gm12.*Gm13;
    for m2=(m1+1):M
        if accel==0
            tmp = Ab(m1,m2);
        else
            sm2 = s(:,:,:,m2);
            tmp = (-sm1.*sm2)*accel + (1-accel)*Ab(m1,m2);
        end
        Gm21 = G(:,:,:,m2,1);
        Gm22 = G(:,:,:,m2,2);
        Gm23 = G(:,:,:,m2,3);
        H1 = H1 + 2*tmp.*Gm11.*Gm21;
        H2 = H2 + 2*tmp.*Gm12.*Gm22;
        H3 = H3 + 2*tmp.*Gm13.*Gm23;
        H4 = H4 + tmp.*(Gm11.*Gm22 + Gm21.*Gm12);
        H5 = H5 + tmp.*(Gm11.*Gm23 + Gm21.*Gm13);
        H6 = H6 + tmp.*(Gm12.*Gm23 + Gm22.*Gm13);
    end
end
clear Gm11 Gm12 Gm13 Gm21 Gm22 Gm23 sm1 sm2 tmp
H = cat(4, H1, H2, H3, H4, H5, H6);
end
%==========================================================================

%==========================================================================
function dat = UpdateWarps(dat,sett)

% Parse function settings
groupwise  = isa(sett.mu,'struct') && isfield(sett.mu,'create');
v_settings = sett.ms.v_settings;
d          = sett.ms.d;

if groupwise
    nw    = GetNumWork(sett,9);
    % Total initial velocity should be zero (Khan & Beg), so mean correct
    avg_v = single(0);
    if nw > 1 && numel(dat) > 1 % PARFOR
        parfor(n=1:numel(dat),nw) avg_v = avg_v + spm_mb_io('GetData',dat(n).v); end
    else % FOR
        for n=1:numel(dat), avg_v = avg_v + spm_mb_io('GetData',dat(n).v); end
    end
    avg_v = avg_v/numel(dat);
    d     = [size(avg_v,1) size(avg_v,2) size(avg_v,3)];
else
    avg_v = [];
end

nw     = GetNumWork(sett,33);
kernel = Shoot(d,v_settings);
if nw > 1 && numel(dat) > 1 % PARFOR
    parfor(n=1:numel(dat),nw) dat(n) = UpdateWarpsSub(dat(n),avg_v,kernel); end
else % FOR
    for n=1:numel(dat), dat(n) = UpdateWarpsSub(dat(n),avg_v,kernel); end
end
end
%==========================================================================

%==========================================================================
function datn = UpdateWarpsSub(datn,avg_v,kernel)
v        = spm_mb_io('GetData',datn.v);
if ~isempty(avg_v)
    v      = v - avg_v;
    datn.v = spm_mb_io('SetData',datn.v,v);
end
u0         = spm_diffeo('vel2mom', v, kernel.v_settings); % Initial momentum
datn.E(2)  = 0.5*sum(u0(:).*v(:));                        % Prior term
psi1       = Shoot(v, kernel, 8);                     % Geodesic shooting
datn.psi   = spm_mb_io('SetData',datn.psi,psi1);
end
%==========================================================================

%==========================================================================
function [dat,mu] = ZoomVolumes(dat,mu,sett,oMmu)

% Parse function settings
d     = sett.ms.d;
Mmu   = sett.ms.Mmu;
nw    = GetNumWork(sett,6+0.5*6);

d0    = [size(mu,1) size(mu,2) size(mu,3)];
z     = single(reshape(d./d0,[1 1 1 3]));
Mzoom = oMmu\Mmu;
%y    = reshape(reshape(Identity(d),[prod(d),3])*Mzoom(1:3,1:3)' + Mzoom(1:3,4)',[d 3]);
y     = reshape(bsxfun(@plus,reshape(Identity(d),[prod(d),3])*Mzoom(1:3,1:3)',Mzoom(1:3,4)'),[d 3]);
if nargout > 1 || ~isempty(mu), mu = spm_diffeo('pullc',mu,y); end % only resize template if updating it

if ~isempty(dat)
    if nw > 1 && numel(dat) > 1 % PARFOR
        parfor(n=1:numel(dat),nw)
            v          = spm_mb_io('GetData',dat(n).v);
           %v          = spm_diffeo('pullc',v.*z,y);
            v          = spm_diffeo('pullc',bsxfun(@times,v,z),y);
            dat(n).v   = ResizeFile(dat(n).v  ,d,Mmu);
            dat(n).psi = ResizeFile(dat(n).psi,d,Mmu);
            dat(n).v   = spm_mb_io('SetData',dat(n).v,v);
        end
    else % FOR
        for n=1:numel(dat)
            v          = spm_mb_io('GetData',dat(n).v);
           %v          = spm_diffeo('pullc',v.*z,y);
            v          = spm_diffeo('pullc',bsxfun(@times,v,z),y);
            dat(n).v   = ResizeFile(dat(n).v  ,d,Mmu);
            dat(n).psi = ResizeFile(dat(n).psi,d,Mmu);
            dat(n).v   = spm_mb_io('SetData',dat(n).v,v);
        end
    end
end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
function I = Horder(d)
I = diag(1:d);
l = d;
for i1=1:d
    for i2=(i1+1):d
        l = l + 1;
        I(i1,i2) = l;
        I(i2,i1) = l;
    end
end
end
%==========================================================================

%==========================================================================
function f = Mask(f,msk)
f(~isfinite(f)) = 0;
%f = f.*msk;
f  = bsxfun(@times,f,msk);
end
%==========================================================================

%==========================================================================
function fin = ResizeFile(fin,d,Mat)
if isa(fin,'char')
    fin = nifti(fin);
end
if isa(fin,'nifti')
    for m=1:numel(fin)
        fin(m).dat.dim(1:3) = d(1:3);
        fin(m).mat  = Mat;
        fin(m).mat0 = Mat;
        create(fin(m));
    end
end
end
%==========================================================================

%==========================================================================
function varargout = Shoot(v0,kernel,args)
% Geodesic shooting
% FORMAT psi = Shoot(v0,kernel,args)
%
% v0       - Initial velocity field n1*n2*n3*3 (single prec. float)
% kernel   - structure created previously
% args     - Integration parameters
%            - [1] Num time steps
%
% psi      - Inverse deformation field n1*n2*n3*3 (single prec. float)
%
% FORMAT kernel = Shoot(d,v_settings)
% d          - dimensions of velocity fields
% v_settings - 8 settings
%              - [1][2][3] Voxel sizes
%              - [4][5][6][7][8] Regularisation settings.
%              Regularisation uses the sum of
%              - [4] - absolute displacements
%              - [5] - laplacian
%              - [6] - bending energy
%              - [7] - linear elasticity mu
%              - [8] - linear elasticity lambda
%
% kernel     - structure encoding Greens function
%
% This code generates inverse deformations from
% initial velocity fields by gedesic shooting.  See the work of Miller,
% Younes and others.
%
% LDDMM (Beg et al) uses the following evolution equation:
%     d\phi/dt = v_t(\phi_t)
% where a variational procedure is used to find the stationary solution
% for the time varying velocity field.
% In principle though, once the initial velocity is known, then the
% velocity at subsequent time points can be computed.  This requires
% initial momentum (u_0), computed (using differential operator L) by:
%     u_0 = L v_0
% Then (Ad_{\phi_t})^* m_0 is computed:
%     u_t = |d \phi_t| (d\phi_t)^T u_0(\phi_t)
% The velocity field at this time point is then obtained by using
% multigrid to solve:
%     v_t = L^{-1} u_t
%
% These equations can be found in:
% Younes (2007). "Jacobi fields in groups of diffeomorphisms and
% applications". Quarterly of Applied Mathematics, vol LXV,
% number 1, pages 113-134 (2007).

if nargin==2
    if numel(v0)>5
        d = [size(v0) 1];
        d = d(1:3);
    else
        d = v0;
    end
    v_settings = kernel;
    spm_diffeo('boundary',0);
    F   = spm_shoot_greens('kernel',d,v_settings);
    varargout{1} = struct('d',d, 'v_settings',v_settings, 'F', F);
    return;
end

if isempty(v0)
    varargout{1} = [];
    varargout{2} = [];
    return;
end

args0 = 8;
if nargin<3
    args = args0;
else
    if numel(args)<numel(args0)
        args = [args args0((numel(args)+1):end)];
    end
end

T     = args(1);   % # Time steps
d     = size(v0);
d     = d(1:3);
id    = Identity(d);

if sum(v0(:).^2)==0
    varargout{1} = id;
    varargout{2} = v0;
end

if ~isfinite(T)
    % Number of time steps from an educated guess about how far to move
    T = double(floor(sqrt(max(max(max(v0(:,:,:,1).^2+v0(:,:,:,2).^2+v0(:,:,:,3).^2)))))+1);
end

spm_diffeo('boundary',0);
v   = v0;
u   = spm_diffeo('vel2mom',v,kernel.v_settings); % Initial momentum (u_0 = L v_0)
psi = id - v/T;

for t=2:abs(T)
    % The update of u_t is not exactly as described in the paper, but describing this might be a bit
    % tricky. The approach here was the most stable one I could find - although it does lose some
    % energy as < v_t, u_t> decreases over time steps.
    Jdp         = spm_diffeo('jacobian',id-v/T);
    u1          = zeros(size(u),'single');
    u1(:,:,:,1) = Jdp(:,:,:,1,1).*u(:,:,:,1) + Jdp(:,:,:,2,1).*u(:,:,:,2) + Jdp(:,:,:,3,1).*u(:,:,:,3);
    u1(:,:,:,2) = Jdp(:,:,:,1,2).*u(:,:,:,1) + Jdp(:,:,:,2,2).*u(:,:,:,2) + Jdp(:,:,:,3,2).*u(:,:,:,3);
    u1(:,:,:,3) = Jdp(:,:,:,1,3).*u(:,:,:,1) + Jdp(:,:,:,2,3).*u(:,:,:,2) + Jdp(:,:,:,3,3).*u(:,:,:,3);
    clear Jdp

    u = spm_diffeo('pushc',u1,id+v/T);

    % v_t \gets L^g u_t
    v = spm_shoot_greens(u,kernel.F,kernel.v_settings); % Convolve with Greens function of L

    if size(v,3)==1, v(:,:,:,3) = 0; end

    % $\psi \gets \psi \circ (id - \tfrac{1}{T} v)$
    % I found that simply using $\psi \gets \psi - \tfrac{1}{T} (D \psi) v$ was not so stable.
    psi          = spm_diffeo('comp',psi,id-v/T);

    if size(v,3)==1, psi(:,:,:,3) = 1; end
end
varargout{1} = psi;
varargout{2} = v;
end
%==========================================================================

%==========================================================================
function sz = ZoomSettings(v_settings, mu, n)
d   = mu.d;
Mmu = mu.Mmu;
[dz{1:n}] = deal(d);
if isfield(mu,'create')
    sz        = struct('Mmu',Mmu,'d',dz, 'v_settings', [0 0 0 0 0],...
                   'mu_settings',[0 0 0]);
else
    sz        = struct('Mmu',Mmu,'d',dz, 'v_settings', [0 0 0 0 0]);
end

scale = 1/abs(det(Mmu(1:3,1:3)));
for i=1:n
    sz(i).d           = ceil(d/(2^(i-1)));
    z                 = d./sz(i).d;
    sz(i).Mmu         = Mmu*[diag(z), (1-z(:))*0.5; 0 0 0 1];
    vx                = sqrt(sum(sz(i).Mmu(1:3,1:3).^2));
    scale_i           = scale*abs(det(sz(i).Mmu(1:3,1:3)));
    % This adhoc value should really be 1, but this gives less extreme
    % warps in the early iterations, which might help the clustering
    % associate the right priors to each tissue class - without
    % warping the priors to the wrong tissue.
    adhoc             = 1.1; 
    sz(i).v_settings  = [vx v_settings*scale_i.^adhoc]; 
    if isfield(mu,'create')
        mu_settings       = mu.create.mu_settings;
        sz(i).mu_settings = [vx mu_settings*scale_i];
    end
end
end
%==========================================================================

%==========================================================================
function X2d = vol2vec(X4d)
d   = [size(X4d) 1 1];
X2d = reshape(X4d,[prod(d(1:3)) d(4)]);
end
%==========================================================================

%==========================================================================
function X4d = vec2vol(X2d,dm)
dm = [dm(:)' 1 1];
if size(X2d,1)~=prod(dm(1:3))
    error('Incompatible dimensions.');
end
X4d = reshape(X2d,[dm(1:3) size(X2d,2)]);
end
%==========================================================================

%==========================================================================
function nw = GetNumWork(sett,NumVol)
% Estimate number of parfor workers from available system RAM
% (if sett.gen.num_workers = -1)

% Parse function settings
MemMax       = 0;          % max memory usage (in MB)
NumWork      = sett.nworker;
dm           = sett.ms.d;  % current template dimensions
K            = sett.K;     % template classes

if MemMax == 0 % default
    try
        % Get memory info automatically (in MB)
        if ispc
            % Windows
            [~,meminfo] = memory;
            MemMax      = meminfo.PhysicalMemory.Available;
        else
            % UNIX
            [~,meminfo] = system('free --mega');
            meminfo     = string(meminfo);
            meminfo     = strsplit(meminfo,' ');
            MemMax      = str2double(meminfo{14}); % field that holds available RAM (MATLAB 2018a)
        end
    catch
        MemMax = 1024;
    end
end

% Get memory requirement (with current template size)
if nargin<2
    NumVol     = (K*(K+1)/2+4*K);
end
NumFloats      = NumVol*prod(dm(1:3));            % times two..we also keep images, etc in memory (rough)
FloatSizeBytes = 4;                               % One float is four bytes (single precision)
MemReq         = (NumFloats*FloatSizeBytes)/1e6;  % to MB
nw             = max(floor(MemMax/MemReq) - 1,0); % Number of parfor workers to use (minus one..for main thread)
if NumWork >= 0
    nw = min(NumWork,nw);
end
end
%==========================================================================

