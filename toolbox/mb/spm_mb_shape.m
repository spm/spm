function varargout = spm_mb_shape(action,varargin)
% Shape model
%
% FORMAT psi0      = spm_mb_shape('affine',d,Mat)
% FORMAT B         = spm_mb_shape('affine_bases',code)
% FORMAT psi       = spm_mb_shape('compose',psi1,psi0)
% FORMAT id        = spm_mb_shape('identity',d)
% FORMAT dat       = spm_mb_shape('init_def',dat,sett.ms)
% FORMAT l         = spm_mb_shape('LSE',mu,ax)
% FORMAT a1        = spm_mb_shape('pull1',a0,psi,r)
% FORMAT [f1,w1]   = spm_mb_shape('push1',f,psi,d,r)
% FORMAT sd        = spm_mb_shape('samp_dens',Mmu,Mn)
% FORMAT varargout = spm_mb_shape('shoot',v0,kernel,args)
% FORMAT mu1       = spm_mb_shape('shrink_template',mu,oMmu,sett)
% FORMAT P         = spm_mb_shape('softmax',mu,ax)
% FORMAT E         = spm_mb_shape('template_energy',mu,sett)
% FORMAT dat       = spm_mb_shape('update_affines',dat,mu,sett)
% FORMAT [mu,dat]  = spm_mb_shape('update_mean',dat, mu, sett)
% FORMAT dat       = spm_mb_shape('update_simple_affines',dat,mu,sett)
% FORMAT [mu,dat]  = spm_mb_shape('update_simple_mean',dat, mu, sett)
% FORMAT dat       = spm_mb_shape('update_velocities',dat,mu,sett)
% FORMAT dat       = spm_mb_shape('update_warps',dat,sett)
% FORMAT [dat,mu]  = spm_mb_shape('zoom_volumes',dat,mu,sett,oMmu)
% FORMAT sz        = spm_mb_shape('zoom_settings', v_settings, mu, n)
% FORMAT psi       = get_def(dat,sett.ms.Mmu)
%__________________________________________________________________________
% Copyright (C) 2019-2020 Wellcome Centre for Human Neuroimaging

% $Id: spm_mb_shape.m 7928 2020-08-14 16:33:20Z john $
[varargout{1:nargout}] = spm_subfun(localfunctions,action,varargin{:});
%==========================================================================

%==========================================================================
function psi0 = affine(d,Mat,samp)
% Generate an affine deformation field
if nargin<3 || numel(samp)~=3
    samp = [1 1 1];
end
d     = [d 1 1];
d     = d(1:3);
id    = identity(d,samp);
psi0  = vec2vol(bsxfun(@plus,vol2vec(id)*Mat(1:3,1:3)',Mat(1:3,4)'),size(id));
if d(3) == 1, psi0(:,:,:,3) = 1; end
%==========================================================================

%==========================================================================
function B = affine_bases(code)
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
%==========================================================================

%==========================================================================
function psi = compose(psi1,psi0)
% Function composition
% \psi \gets \psi^1 \circ \psi^0
if isempty(psi1)
    psi = psi0;
else
    bc = spm_diffeo('bound');
    spm_diffeo('bound',1);
    psi = spm_diffeo('comp',psi1,psi0);
    spm_diffeo('bound',bc);
end
if size(psi,3) == 1, psi(:,:,:,3) = 1; end % 2D
%==========================================================================

%==========================================================================
function id = identity(d,samp)
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
%==========================================================================

%==========================================================================
function l = LSE(mu,ax)
% log-sum-exp function

if nargin<2, ax = 4; end
mx = max(mu,[],ax);
l  = log(exp(-mx) + sum(exp(bsxfun(@minus,mu,mx)),ax)) + mx;
%==========================================================================

%==========================================================================
function a1 = pull1(a0,psi,r)
% Resample an image or set of images
% FORMAT a1 = pull1(a0,psi,r)
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
    id = identity(size(psi));
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
%==========================================================================

%==========================================================================
function [f1,w1] = push1(f,psi,d,r)
% push an image (or set of images) accorging to a spatial transform
% FORMAT [f1,w1] = push1(f,psi,d,r)
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

    id    = identity(size(psi));
    f1    = single(0);
    w1    = single(0);
    for dz=zrange
        for dy=yrange
            for dx=xrange
                ids       = bsxfun(@plus, id, cat(4,dx,dy,dz));
                psi1      = spm_diffeo('pullc',psi-id,    ids)+ids;
                fs        = spm_diffeo('pull',single(f), ids);
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
%==========================================================================

%==========================================================================
function r = range(n)
r = (-floor((n-1)/2):ceil((n-1)/2))/n;
%==========================================================================

%==========================================================================
function sd = samp_dens(Mmu,Mn)
vx_mu = sqrt(sum(Mmu(1:3,1:3).^2,1));
vx_f  = sqrt(sum( Mn(1:3,1:3).^2,1));
sd    = max(round(2.0*vx_f./vx_mu),1);
%==========================================================================

%==========================================================================
function mu = shrink_template(mu,oMmu,sett,smo_wt)

% Parse function settings
if nargin<4, smo_wt = 1; end
d     = sett.ms.d;
Mmu   = sett.ms.Mmu;
d0    = [size(mu,1) size(mu,2) size(mu,3)];
Mzoom = Mmu\oMmu;
if any(d0~=d) || norm(Mzoom-eye(4))>1e-4
    y      = reshape(bsxfun(@plus,...
             reshape(identity(d0),[prod(d0),3])*Mzoom(1:3,1:3)',...
             Mzoom(1:3,4)'),[d0 3]);
    mu     = exp(bsxfun(@minus,mu,LSE(mu,4)));
    [mu,c] = push1(mu,y,d,1);
    e      = eps('single');
    mu     = bsxfun(@rdivide,mu,max(c,e));
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
%==========================================================================

%==========================================================================
function P = softmax(mu,ax)
% safe softmax function

if nargin<2, ax = 4; end
mx  = max(mu,[],ax);
E   = exp(bsxfun(@minus,mu,mx));
den = sum(E,ax)+exp(-mx);
P   = bsxfun(@rdivide,E,den);
%==========================================================================

%==========================================================================
function E = template_energy(mu,mu_settings)
% mu(:)'*kron(eye(K)-1/(K+1),L)*mu(:), where L is the vel2mom regulariser
if ~isempty(mu_settings)
    spm_field('bound',1);
    g  = reg_mu(mu, mu_settings);
    E  = 0.5*mu(:)'*g(:);
else
    E  = 0;
end
%==========================================================================

%==========================================================================
function dat = update_affines(dat,mu,sett)

% Parse function settings
B         = sett.B;
groupwise = isa(sett.mu,'struct') && isfield(sett.mu,'create');
nw        = get_num_workers(sett,max(27,sett.K*5+17));

% Update the affine parameters
if ~isempty(B)
    if nw > 1 && numel(dat) > 1 % PARFOR
        parfor(n=1:numel(dat),nw) dat(n) = update_affines_sub(dat(n),mu,sett); end
    else % FOR
        for n=1:numel(dat), dat(n) = update_affines_sub(dat(n),mu,sett); end
    end

    if groupwise
        % Zero-mean the affine parameters
        mq = sum(cat(2,dat(:).q),2)/numel(dat);
        for n=1:numel(dat)
            dat(n).q = dat(n).q - mq;
        end
    end

    % Update orientations in deformation headers when appropriate
    for n=1:numel(dat)
        Mdef       = spm_dexpm(dat(n).q,B)\sett.ms.Mmu;
        dat(n).psi = spm_mb_io('save_mat',dat(n).psi,Mdef);
    end
end
%==========================================================================

%==========================================================================
function datn = update_affines_sub(datn,mu,sett)
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
psi0 = affine(df,Mmu\Mr*Mn,samp);
ds   = [size(psi0,1),size(psi0,2),size(psi0,3)];
psi1 = get_def(datn,sett.ms.Mmu);
if ~isempty(psi1)
    J    = spm_diffeo('jacobian',psi1);
    J    = reshape(pull1(reshape(J,[d 3*3]),psi0),[ds 3 3]);
    psi  = compose(psi1,psi0);
else
    J    = [];
    psi  = psi0;
end
clear psi0  psi1

mu1  = pull1(mu,psi);
[f,datn] = spm_mb_classes(datn,mu1,sett);
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
            tmp = Gm{i1};
        end
        tmp(~isfinite(tmp)) = 0;
        G(:,:,:,m,i1) = tmp;
    end
    clear Gm
end
clear J mu

msk    = all(isfinite(f),4) & all(isfinite(mu1),4);
mu1(~isfinite(mu1)) = 0;
a      = mask(f - softmax(mu1,4),msk);
[H,g]  = affine_hessian(mu1,G,a,single(msk),accel,samp);
g      = double(dM'*g);
H      = dM'*H*dM;
H      = H + eye(numel(q))*(norm(H)*1e-5 + 0.01);
q      = q + scal*(H\g);
datn.q = q;
%==========================================================================

%==========================================================================
function [H,g] = affine_hessian(mu,G,a,w,accel,samp)
if nargin<6, samp = [1 1 1]; end
d  = [size(mu,1),size(mu,2),size(mu,3)];
I  = horder(3);
H  = zeros(12,12);
g  = zeros(12, 1);
[x{1:4}] = ndgrid(((1:d(1))-1)*samp(1)+1,((1:d(2))-1)*samp(2)+1,1,1);
for i=1:d(3)
    x{3} = x{3}*0+(i-1)*samp(3)+1;
    gv   = reshape(sum(bsxfun(@times,a(:,:,i,:),G(:,:,i,:,:)),4),[d(1:2) 1 3]);
    Hv   = bsxfun(@times,w(:,:,i),velocity_hessian(mu(:,:,i,:),G(:,:,i,:,:),accel));
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
%==========================================================================

%==========================================================================
function [mu,dat] = update_mean(dat, mu, sett)

% Parse function settings
accel       = sett.accel;
mu_settings = sett.ms.mu_settings;
nw          = get_num_workers(sett,4*sett.K+4);

spm_field('bound',1);
g  = reg_mu(mu,mu_settings);
w  = zeros(sett.ms.d,'single');
if nw > 1 && numel(dat) > 1 % PARFOR
    parfor(n=1:numel(dat),nw)
        [gn,wn,dat(n)] = update_mean_sub(dat(n),mu,sett);
        g              = g + gn;
        w              = w + wn;
    end
else
    for n=1:numel(dat) % FOR
        [gn,wn,dat(n)] = update_mean_sub(dat(n),mu,sett);
        g              = g + gn;
        w              = w + wn;
    end
end
clear gn wn
mu = gn_mu_update(mu,g,w,mu_settings,accel);
%==========================================================================

%==========================================================================
function mu = gn_mu_update(mu,g,w,mu_settings,accel)
% Use Gauss-Seidel method for computing updates
% https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
spm_field('bound',1);
K   = size(mu,4);
update_settings = [mu_settings(1:3) mu_settings(4:end)*(1-1/(K+1)) 2 2];
if accel>0, s   = softmax(mu); end
dmu = zeros(size(mu),'like',mu);
for it=1:16
    for k=1:K

        % Diagonal elements of Hessian
        if accel==0
            h_kk = (0.5*(1-1/(K+1)))*w;
        else
            h_kk = (accel*(s(:,:,:,k)-s(:,:,:,k).^2) + (1-accel)*0.5*(1-1/(K+1))).*w;
        end

        % Dot product betwen off-diagonals of likelihood Hessian and dmu
        g_k = zeros([size(mu,1),size(mu,2),size(mu,3)],'like',g);
        for k1=1:K
            if k1~=k
                % Off-diaginal elements of Hessian
                if accel==0
                    g_k = g_k - dmu(:,:,:,k1).*(0.5/(K+1));
                else
                    g_k = g_k - dmu(:,:,:,k1).*(accel*(s(:,:,:,k).*s(:,:,:,k1)) + ((1-accel)*0.5/(K+1)));
                end
            end
        end
        g_k          = g_k.*w;

        % Dot product between off-diagonals of regularisation Hessian and dmu
        g_k          = g_k - spm_field('vel2mom', (sum(dmu,4)-dmu(:,:,:,k))/(K+1), mu_settings);

        % Gauss-Seidel update
        dmu(:,:,:,k) = spm_field(h_kk, g(:,:,:,k) - g_k, update_settings);
    end
end
mu = mu - dmu;
%==========================================================================

%==========================================================================
function g = reg_mu(mu,mu_settings)
K = size(mu,4);
g = spm_field('vel2mom', bsxfun(@minus,mu, sum(mu,4)/(K+1)), mu_settings);
%==========================================================================

%==========================================================================
function [g,w,datn] = update_mean_sub(datn,mu,sett)

% Parse function settings
B        = sett.B;
d        = sett.ms.d;
Mmu      = sett.ms.Mmu;
df       = datn.dm;
q        = double(datn.q);
Mn       = datn.Mat;
samp     = datn.samp;

psi      = compose(get_def(datn,sett.ms.Mmu),affine(df, Mmu\spm_dexpm(q,B)*Mn,samp));
mu       = pull1(mu,psi);
[f,datn] = spm_mb_classes(datn,mu,sett);
[g,w]    = push1(softmax(mu,4) - f,psi,d,1);
%==========================================================================

%==========================================================================
function H = appearance_hessian(mu,accel,w)
M  = size(mu,4);
d  = [size(mu,1) size(mu,2) size(mu,3)];
if accel>0, s  = softmax(mu,4); end
Ab = 0.5*(eye(M)-1/(M+1)); % See Bohning's paper
I  = horder(M);
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
%==========================================================================

%==========================================================================
function dat = update_simple_affines(dat,mu,sett)

% Parse function settings
accel     = sett.accel;
B         = sett.B;
if ~isempty(B)
    groupwise = isa(sett.mu,'struct') && isfield(sett.mu,'create');
    nw        = get_num_workers(sett,4*sett.K+4);

    % Update the affine parameters
    G  = spm_diffeo('grad',mu);
    H0 = velocity_hessian(mu,G,accel);

    if nw > 1 && numel(dat) > 1 % PARFOR
        parfor(n=1:numel(dat),nw)
            dat(n) = update_simple_affines_sub(dat(n),mu,G,H0,sett);
        end
    else
        for n=1:numel(dat)
            dat(n) = update_simple_affines_sub(dat(n),mu,G,H0,sett);
        end
    end

    if groupwise
        % Zero-mean the affine parameters
        mq = sum(cat(2,dat(:).q),2)/numel(dat);
        for n=1:numel(dat)
            dat(n).q = dat(n).q - mq;
        end
    end

    % Update orientations in deformation headers when appropriate
    for n=1:numel(dat)
        Mdef       = spm_dexpm(dat(n).q,B)\sett.ms.Mmu;
        dat(n).psi = spm_mb_io('save_mat',dat(n).psi,Mdef);
    end
end
%==========================================================================

%==========================================================================
function datn = update_simple_affines_sub(datn,mu,G,H0,sett)

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

psi      = affine(df,Mmu\Mr*Mn,samp);
mu1      = pull1(mu,psi);

[f,datn] = spm_mb_classes(datn,mu1,sett);
[a,w]    = push1(f - softmax(mu1,4),psi,d,1);
clear mu1 psi f

[H,g]    = simple_affine_hessian(mu,G,H0,a,w);
g        = double(dM'*g);
H        = dM'*H*dM;
H        = H + eye(numel(q))*(norm(H)*1e-5 + 0.1);
q        = q + scal*(H\g);
datn.q   = q;
%==========================================================================

%==========================================================================
function [H,g] = simple_affine_hessian(mu,G,H0,a,w)
d  = [size(mu,1),size(mu,2),size(mu,3)];
I  = horder(3);
H  = zeros(12,12);
g  = zeros(12, 1);
[x{1:4}] = ndgrid(1:d(1),1:d(2),1,1);
for i=1:d(3)
    x{3} = x{3}*0+i;
    gv   = reshape(sum(bsxfun(@times,a(:,:,i,:),G(:,:,i,:,:)),4),[d(1:2) 1 3]);
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
%==========================================================================

%==========================================================================
function [mu,dat] = update_simple_mean(dat, mu, sett)
% Unused functionality

% Parse function settings
accel       = sett.accel;
mu_settings = sett.ms.mu_settings;
s_settings  = [2 2];
nw          = get_num_workers(sett,4*sett.K+4);

w  = zeros(sett.ms.d,'single');
gf = zeros(size(mu),'single');
if nw > 1 && numel(dat) > 1 % PARFOR
    parfor(n=1:numel(dat),nw)
        [gn,wn,dat(n)] = update_simple_mean_sub(dat(n),mu,sett);
        gf             = gf + gn;
        w              = w  + wn;
    end
else
    for n=1:numel(dat) % FOR
        [gn,wn,dat(n)] = update_simple_mean_sub(dat(n),mu,sett);
        gf             = gf + gn;
        w              = w  + wn;
    end
end
clear gn wn
for it=1:ceil(4+2*log2(numel(dat)))
    H  = appearance_hessian(mu,accel,w);
    g  = w.*softmax(mu,4) - gf;
    spm_field('bound',1);
    g  = g  + reg_mu(mu, mu_settings);
    % Note that spm_field could be re-written to make the updates slightly
    % more effective.
    mu = mu - spm_field(H, g, [mu_settings s_settings]);
end
%==========================================================================

%==========================================================================
function [g,w,datn] = update_simple_mean_sub(datn,mu,sett)
% Unused functionality

% Parse function settings
B     = sett.B;
d     = sett.ms.d;
Mmu   = sett.ms.Mmu;
df    = datn.dm;
q     = double(datn.q);
Mn    = datn.Mat;
samp  = datn.samp;

psi      = compose(get_def(datn,sett.ms.Mmu),affine(df,Mmu\spm_dexpm(q,B)*Mn,samp));
mu       = pull1(mu,psi);
[f,datn] = spm_mb_classes(datn,mu,sett);
[g,w]    = push1(f,psi,d,1);
%==========================================================================

%==========================================================================
function dat = update_velocities(dat,mu,sett)

% Parse function settings
accel = sett.accel;
nw    = get_num_workers(sett,4*sett.K+4*3);

G  = spm_diffeo('grad',mu);
H0 = velocity_hessian(mu,G,accel);
if size(G,3) == 1
    % Data is 2D -> add some regularisation
    H0(:,:,:,3) = H0(:,:,:,3) + mean(reshape(H0(:,:,:,[1 2]),[],1));
end
if nw > 1 && numel(dat) > 1 % PARFOR
    parfor(n=1:numel(dat),nw) dat(n) = update_velocities_sub(dat(n),mu,G,H0,sett); end
else % FOR
    for n=1:numel(dat), dat(n) = update_velocities_sub(dat(n),mu,G,H0,sett); end
end
%==========================================================================

%==========================================================================
function datn = update_velocities_sub(datn,mu,G,H0,sett)

% Parse function settings
B          = sett.B;
d          = sett.ms.d;
Mmu        = sett.ms.Mmu;
s_settings = [2 2];
scal       = 0.9;
v_settings = sett.ms.v_settings;

v         = spm_mb_io('get_data',datn.v);
q         = datn.q;
Mn        = datn.Mat;
samp      = datn.samp;
Mr        = spm_dexpm(q,B);
Mat       = Mmu\Mr*Mn;
df        = datn.dm;
psi       = compose(get_def(datn,sett.ms.Mmu),affine(df,Mat,samp));
mu        = pull1(mu,psi);
[f,datn]  = spm_mb_classes(datn,mu,sett);
[a,w]     = push1(f - softmax(mu,4),psi,d,1);
clear psi f mu

g         = reshape(sum(bsxfun(@times,a,G),4),[d 3]);
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
datn.v = spm_mb_io('set_data',datn.v,v);
%==========================================================================

%==========================================================================
function H = velocity_hessian(mu,G,accel)
d  = [size(mu,1),size(mu,2),size(mu,3)];
M  = size(mu,4);
if accel>0, s  = softmax(mu,4); end
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
        H1 = H1 + 2*tmp.* Gm11.*Gm21;
        H2 = H2 + 2*tmp.* Gm12.*Gm22;
        H3 = H3 + 2*tmp.* Gm13.*Gm23;
        H4 = H4 +   tmp.*(Gm11.*Gm22 + Gm21.*Gm12);
        H5 = H5 +   tmp.*(Gm11.*Gm23 + Gm21.*Gm13);
        H6 = H6 +   tmp.*(Gm12.*Gm23 + Gm22.*Gm13);
    end
end
clear Gm11 Gm12 Gm13 Gm21 Gm22 Gm23 sm1 sm2 tmp
H = cat(4, H1, H2, H3, H4, H5, H6);
%==========================================================================

%==========================================================================
function dat = update_warps(dat,sett)

% Parse function settings
groupwise  = isa(sett.mu,'struct') && isfield(sett.mu,'create');
v_settings = sett.ms.v_settings;
d          = sett.ms.d;

if groupwise
    nw    = get_num_workers(sett,9);
    % Total initial velocity should be zero (Khan & Beg), so mean correct
    avg_v = single(0);
    if nw > 1 && numel(dat) > 1 % PARFOR
        parfor(n=1:numel(dat),nw) avg_v = avg_v + spm_mb_io('get_data',dat(n).v); end
    else % FOR
        for n=1:numel(dat), avg_v = avg_v + spm_mb_io('get_data',dat(n).v); end
    end
    avg_v = avg_v/numel(dat);
    d     = [size(avg_v,1) size(avg_v,2) size(avg_v,3)];
else
    avg_v = [];
end

nw     = get_num_workers(sett,33);
kernel = shoot(d,v_settings);
if nw > 1 && numel(dat) > 1 % PARFOR
    parfor(n=1:numel(dat),nw) dat(n) = update_warps_sub(dat(n),avg_v,kernel,sett); end
else % FOR
    for n=1:numel(dat), dat(n) = update_warps_sub(dat(n),avg_v,kernel,sett); end
end
%==========================================================================

%==========================================================================
function datn = update_warps_sub(datn,avg_v,kernel,sett)
v          = spm_mb_io('get_data',datn.v);
if ~isempty(avg_v)
    v      = v - avg_v;
    datn.v = spm_mb_io('set_data',datn.v,v);
end
u0         = spm_diffeo('vel2mom', v, kernel.v_settings); % Initial momentum
datn.E(2)  = 0.5*sum(u0(:).*v(:));                        % Prior term
psi1       = shoot(v, kernel, 8);                         % Geodesic shooting
datn       = set_def(datn,sett.ms.Mmu,psi1);
%==========================================================================

%==========================================================================
function [dat,mu] = zoom_volumes(dat,mu,sett,oMmu)

% Parse function settings
d     = sett.ms.d;
Mmu   = sett.ms.Mmu;
B     = sett.B;
nw    = get_num_workers(sett,6+0.5*6);

d0    = [size(mu,1) size(mu,2) size(mu,3)];
z     = single(reshape(d./d0,[1 1 1 3]));
Mzoom = oMmu\Mmu;
y     = reshape(bsxfun(@plus,reshape(identity(d),[prod(d),3])*Mzoom(1:3,1:3)',Mzoom(1:3,4)'),[d 3]);
if nargout > 1 || ~isempty(mu), mu = spm_diffeo('pullc',mu,y); end % only resize template if updating it

if ~isempty(dat)
    if nw > 1 && numel(dat) > 1 % PARFOR
        parfor(n=1:numel(dat),nw)
            v          = spm_mb_io('get_data',dat(n).v);
            v          = spm_diffeo('pullc',bsxfun(@times,v,z),y);
            dat(n).v   = resize_file(dat(n).v  ,d,Mmu);
            dat(n).v   = spm_mb_io('set_data',dat(n).v,v);
            if ~isempty(B)
                Mdef   = spm_dexpm(dat(n).q,B)\Mmu;
            else
                Mdef   = Mmu;
            end
            dat(n).psi = resize_file(dat(n).psi,d,Mdef);
        end
    else % FOR
        for n=1:numel(dat)
            v          = spm_mb_io('get_data',dat(n).v);
            v          = spm_diffeo('pullc',bsxfun(@times,v,z),y);
            dat(n).v   = resize_file(dat(n).v  ,d,Mmu);
            dat(n).v   = spm_mb_io('set_data',dat(n).v,v);
            if ~isempty(B)
                Mdef   = spm_dexpm(dat(n).q,B)\Mmu;
            else
                Mdef   = Mmu;
            end
            dat(n).psi = resize_file(dat(n).psi,d,Mdef);
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
function I = horder(d)
I = diag(1:d);
l = d;
for i1=1:d
    for i2=(i1+1):d
        l = l + 1;
        I(i1,i2) = l;
        I(i2,i1) = l;
    end
end
%==========================================================================

%==========================================================================
function f = mask(f,msk)
f(~isfinite(f)) = 0;
f  = bsxfun(@times,f,msk);
%==========================================================================

%==========================================================================
function varargout = shoot(v0,kernel,args)
% Geodesic shooting
% FORMAT psi = shoot(v0,kernel,args)
%
% v0       - Initial velocity field n1*n2*n3*3 (single prec. float)
% kernel   - structure created previously
% args     - Integration parameters
%            - [1] Num time steps
%
% psi      - Inverse deformation field n1*n2*n3*3 (single prec. float)
%
% FORMAT kernel = shoot(d,v_settings)
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
id    = identity(d);

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
%==========================================================================

%==========================================================================
function sz = zoom_settings(v_settings, mu, n)
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
%==========================================================================

%==========================================================================
function X2d = vol2vec(X4d)
d   = [size(X4d) 1 1];
X2d = reshape(X4d,[prod(d(1:3)) d(4)]);
%==========================================================================

%==========================================================================
function X4d = vec2vol(X2d,dm)
dm = [dm(:)' 1 1];
if size(X2d,1)~=prod(dm(1:3))
    error('Incompatible dimensions.');
end
X4d = reshape(X2d,[dm(1:3) size(X2d,2)]);
%==========================================================================

%==========================================================================
function nw = get_num_workers(sett,NumVol)
% Estimate number of parfor workers from available system RAM
% (if sett.gen.num_workers = -1)

% Parse function settings
MemMax         = 0;          % max memory usage (in MB)
NumWork        = sett.nworker;
dm             = sett.ms.d;  % current template dimensions
K              = sett.K;     % template classes

if NumWork <= 1
    nw = 0;
    return
end

if MemMax == 0
    MemMax     = spm_platform('memory','available') / 1024 / 1024;
    if isnan(MemMax)
        MemMax = 1024;
    end
end

% Get memory requirement (with current template size)
if nargin<2
    NumVol     = (K*(K+1)/2+4*K);
end
NumFloats      = NumVol*prod(dm(1:3));            % float size of mean Hessian + padding (we also keep images, etc in memory (rough))
FloatSizeBytes = 4;                               % One float is four bytes (single precision)
MemReq         = (NumFloats*FloatSizeBytes)/1e6;  % to MB
nw             = max(floor(MemMax/MemReq) - 1,0); % Number of parfor workers to use (minus one..for main thread)
if NumWork >= 0
    nw = min(NumWork,nw);
end
%==========================================================================

%==========================================================================
function phi = MatDefMul(phi,M)
d   = size(phi);
phi = reshape(bsxfun(@plus,reshape(phi,[prod(d(1:3)),3])*M(1:3,1:3)',M(1:3,4)'),d);
%==========================================================================

%==========================================================================
% Functions that should be integrated better with spm_mb_io.m.
% This should include anything that relies on the nifti data
% structure.
%==========================================================================

%==========================================================================
function psi = get_def(datn,Mmu)
psi    = spm_mb_io('get_data',datn.psi);
if isa(datn.psi,'nifti')
    psi = MatDefMul(psi,inv(Mmu));
end
%==========================================================================

%==========================================================================
function datn = set_def(datn,Mmu,psi)
if isa(datn.psi,'nifti')
    psi = MatDefMul(psi,Mmu);
end
datn.psi = spm_mb_io('set_data',datn.psi,psi);
%==========================================================================

%==========================================================================
function fin = resize_file(fin,d,Mat)
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
%==========================================================================

%==========================================================================
function dat = init_def(dat,ms)
% Possibly break up this function and put parts of it in the io file.

% Parse function settings
d       = ms.d;
Mmu     = ms.Mmu;

v    = zeros([d,3],'single');
psi1 = identity(d);
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
        nii.descrip = 'Deformation';
        nii.mat     = Mmu;
        create(nii);
        dat(n).psi  = nii;
        dat(n)      = set_def(dat(n),Mmu,psi1);
    end
end
%==========================================================================
