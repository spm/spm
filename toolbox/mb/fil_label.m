function files = fil_label(fil,mbsett,mbdat,iterations,vsett_scale,odir,df,Mf)
% Label image(s)
% FORMAT files = fil_label(fil,mbsett,mbdat,iterations,vsett_scale,odir,df,Mf)
% fil         - a trained model (see fil_train) loaded with
%                   fil    = load('fil_blah.mat');
% mbsett      - global parameters from mb toolbox registration
% mbdat       - subject data from mb toolbox registration
%                   mb     = load('mb_blah.mat');
%                   mbsett = mb.sett;
%                   mbdat  = mb.dat;
% iterations  - three elements containing
%                   Number of registration Gauss-Newton updates
%                   Number of outer iterations to update the latent vars
%                   Number of inner iterations to update the latent vars
%               (defaults to [6 10 10])
% vsett_scale - scaling of the regularisation, relative to what was used
%               originally by the mb toolbox (defaults to 0.25)
% odir        - output directory name (defaults to '.')
% df          - dimensions of label image (optional)
% Mf          - voxel-to-world matrix of label image (optional)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


if nargin<4, iterations  = [6 10 10]; end
if nargin<5, vsett_scale = 1/4; end
if nargin<6, odir = '.'; end
if nargin<7, df = []; end
if nargin<8, Mf = []; end
v_settings = mbsett.ms.v_settings;
v_settings(4:end) = v_settings(4:end)*vsett_scale;

files = '';
for n=1:numel(mbdat)
    onam  = fullfile(odir,['label_' mbdat(n).onam '.nii']);
    label_one(fil, mbdat(n), mbsett, iterations, v_settings,1:3,df,Mf,onam);
    files = strvcat(files,onam);
end


function write_image(Y,Mat,pth)
No         = nifti;
No.dat     = file_array(pth,size(Y),'UINT8',0,1,0);
No.mat     = Mat;
No.mat0    = Mat;
No.descrip = 'Neuromorphometrics labels';
create(No);
No.dat(:,:,:) = Y;


function Y = label_one(fil, datn, mbsett, niter, v_settings,sel, df,Mf, onam)
[Z,phi1,R1] = estimate_params(fil, datn, mbsett, niter, v_settings,sel);
sh       = spm_mb_shape;
if isempty(df), df = datn.dm;  end
if isempty(Mf), Mf = datn.Mat; end
phi      = sh.compose(phi1,sh.affine(df,R1*Mf));
Y        = project(fil.model,Z,phi);
write_image(Y,Mf,onam);


function [Z,phi1,R1] = estimate_params(fil, datn, mbsett, niter, v_settings,sel)
io     = spm_mb_io;
sh     = spm_mb_shape;
ind    = 1;
Ml     = fil.Mw;         % Label space voxel-to-world
dl     = fil.dw;         % Label space dimensions
dmu    = mbsett.ms.d;    % Template space dimensions
Mmu    = mbsett.ms.Mmu;  % Template space voxel-to-world
kernel = sh.shoot(dmu,v_settings); % Greens kernel
df     = datn.dm;       % Native space dimensions

% Should really account for the Jacobians here, in case there are rotations/zooms/etc
if norm(Mmu(1:3,1:3)\Ml(1:3,1:3) - eye(3))>1e-3
    error('Does not handle Mmu and Ml of this form. See code for more info.');
end

% Map from native voxels to rotated versions in coordinate system of template
R1     = Mmu\spm_dexpm(double(datn.q),mbsett.B);
R      = R1*datn.Mat;
v      = io.get_data(datn.v);

phi1   = sh.get_def(datn,mbsett.ms.Mmu);   % Load diffeomorphism
%mu    = io.get_data(mbsett.mu.exist.mu);  % Load mb-template
%phi   = sh.compose(phi1,sh.affine(df,R)); % Mapping from native voxels to template voxels
%mu1   = sh.pull1(mu,phi);                 % Template warped to individual
%fl    = likelihoods(datn,mbsett,mu1);
%f     = responsibilities(datn,mbsett,mu1,sel); % Native-space tissue probabilities
f      = responsibilities(datn,mbsett,sh.pull1(io.get_data(mbsett.mu.exist.mu),sh.compose(phi1,sh.affine(df,R))),sel);

phi1   = sh.matdefmul(phi1,Ml\Mmu);        % Now maps to voxels in label space
phi    = sh.compose(phi1,sh.affine(df,R)); % Mapping from native voxels to label space voxels
%f     = classes(fl, sh.pull1(reconmu(fil.model,ind),phi));
[f1,c] = sh.push1(f,phi,dl);               % Push native-space responsibilities to label-space
Z      = latent(fil.model,f1,c,ind,niter(2:3));    % Compute latent variables

for it=1:niter(1)
    mu    = reconmu(fil.model,ind,Z);       % Reconstruct label-space mu
    spm_diffeo('boundary',1);               % Zero gradient at edges
    G0    = spm_diffeo('grad',mu);          % Gradients of label-space mu
    H0    = sh.velocity_hessian(mu, G0, datn.delta, 0.8); % Label-space Hessian
   %f     = classes(fl,mu1);


    % Push differences from native-space to label-space
    % This approach was used in the validations
   %[a,c]  = sh.push1(f - sh.softmax0(sh.pull1(mu,phi),4),phi,dl,1);

    % Compute gradients that match the objective function of the labelling
    [f1,c] = sh.push1(f,phi,dl,1);
    a      = f1 - bsxfun(@times,c, sh.softmax0(mu,4));

    % Note that the affine mapping should ideally be translation by integer numbers of voxels
    g     = sh.push1(reshape(sum(bsxfun(@times,a,G0),4),[dl 3]), sh.affine(dl,Mmu\Ml),dmu); % Registration gradients
    H     = sh.push1(bsxfun(@times,c,H0),                        sh.affine(dl,Mmu\Ml),dmu); % Registration Hessian

    % Update velocity field in template-space
    spm_diffeo('boundary',0);
    u0    = spm_diffeo('vel2mom', v, v_settings);               % Initial momentum
    v     = v - spm_diffeo('fmg', H, g + u0, [v_settings 2 2]); % Gauss-Newton update

    % Reconstruct mapping from native-space voxels to label-space voxels
   %phi1  = sh.shoot(v, kernel, 8);           % Mapping to voxels of mu
   %phi1  = sh.matdefmul(phi1,Ml\Mmu);        % Mapping to voxels in labelling space 
   %phi   = sh.compose(phi1,sh.affine(df,R)); % Account for rigid transform
    phi1  = sh.matdefmul(sh.shoot(v, kernel, 8),Ml\Mmu);
    phi   = sh.compose(phi1,sh.affine(df,R));

   %f      = classes(fl,sh.pull1(mu,phi));
    [f1,c] = sh.push1(f,phi,dl);                    % Push responsibilities to label-space
    Z      = latent(fil.model,f1,c,ind,niter(2:3)); % Update the latent variables
end



function Ynative = project(model,Z,phi) %C,P) % later
% Reconstruct native-space data
% FORMAT Y = PatchCCAnative(model,Z,phi)
% model - The learned model
% Z     - Cell array of latent variables
% phi   - Mapping from native to template space
%
% Y     - Native-space patch image
%

ind = 2; % Use the 2nd view from the model data
dm_temp = cellfun(@max,model(end,end,end).pos); % Size of template-space data

% Find what label values are used
labels = [];
for p=1:numel(model)
    labels = unique([labels(:); model(p).c{ind}(:)]);
end
if max(labels)>255, error('Maximum label value too big.'); end

spm_progress_bar('Init',numel(labels),sprintf('Labelling %dx%dx%d',dm_temp),'Labels complete');

% Compute a template-space image (L) of log-sum-exp (lse) results
L = zeros(dm_temp,'single'); % Map of lse values
for p=1:numel(model)         % Loop over patches
    patch     = model(p);                          % Model for current patch
    dm_patch  = cellfun(@numel,patch.pos);         % Dimensions of patch
    L(patch.pos{:}) = getLSE0(Z{p},patch.mod(ind),dm_patch); % log-sum-exp from variables
end

% Native-space data
dm_native = [size(phi,1) size(phi,2) size(phi,3)]; % Native-space dimensions
Ynative   = zeros(dm_native,'uint8');  % Native-space labels
Pnative   = zeros(dm_native,'single'); % Native-space maximum label probability
%Cend     = 1-sum(C,4); % later

ii = 0;
for label=labels' % Loop over labels

    % Generate template-space probability map for this label
    Pl = zeros(dm_temp,'single');
    for p=1:numel(model)                 % Loop over patches
        patch     = model(p);            % Model for current patch
        selection = patch.c{ind}==label; % Where current label is encoded in the patch
        if any(selection)
            m        = find(selection);
            dm_patch = cellfun(@numel,patch.pos);              % Dimensions of patch
            psi      = getPsi(Z{p},patch.mod(ind),dm_patch,m); % Linear combination of bases
            Pl(patch.pos{:}) = exp(psi-L(patch.pos{:}));       % Softmax values for this label
        end
    end

    % Warp label probabilities to native space
    Pl = spm_diffeo('pull',Pl,phi);

   %% Code for if/when training uses the model as priors for tissue classification
   %Pc = Cend*P(label+1,end);
   %for i=1:size(C,4)
   %    Pc = Pc + C(:,:,:,i)*P(label,i); %%%% +1 ?
   %end
   %Pl = Pl.*Pc;

    % Update native space data
    replace          = find(Pl > Pnative); % Find native voxels to replace
    Pnative(replace) = Pl(replace);        % Update maximum probabilities
    Ynative(replace) = label-1;            % Update native space labels

    ii = ii+1;
    spm_progress_bar('Set',ii);
end
spm_progress_bar('Clear');



function psi = getPsi(z,mod,dm,m)
K   = size(mod.W,3);
M   = size(mod.W,2);
if nargin>=4
    % Reconstruct psi using only mod.W(:,m,:) and mod.mu(:,m) and reshape to dm
    if m==M+1
        psi = zeros(dm);
    else
        psi = reshape(reshape(mod.W(:,m,:),[prod(dm),K])*z,dm)+reshape(mod.mu(:,m),dm);
    end
else
    psi = reshape(reshape(mod.W,[prod(dm)*M,K])*z+mod.mu(:),[dm M]);
end


function lse = getLSE0(z,mod,dm)
% Compute log-sum-exp
K   = size(mod.W,3);
M   = size(mod.W,2);
if M>0
    psi = reshape(reshape(mod.W,[prod(dm)*M,K])*z,[dm M])+reshape(mod.mu,[dm M]);
    mx  = max(max(psi,[],4),0);
    lse = log(sum(exp(bsxfun(@minus,psi,mx)),4)+exp(-mx))+mx;
else
    lse = zeros(dm);
end


function fl = responsibilities(datn,mbsett,mu,select)
app    = spm_mb_appearance;
io     = spm_mb_io;
gmm    = datn.model.gmm;
fn     = io.get_image(gmm);
gmms   = mbsett.gmm(gmm.pop);
C      = gmms.C;
mg_ix  = gmms.mg_ix;
gam    = gmm.gam;
chan   = app.inu_basis(gmm.T,size(fn),datn.Mat,ones(1,C));
[~,mf,vf] = app.inu_recon(fn,[],chan,gmm.T,gmm.Sig);

mu  = spm_mb_classes('template_k1',mu,datn.delta);
fl  = spm_gmmlib('resp',gmm.m,gmm.b,gmm.W,gmm.nu,gam,...
                 uint64(mg_ix), mu,mf,vf, uint64(gmm.samp));
fl  = fl(:,:,:,select);


function [Z,model] = latent(model,F,J,ind,nit,Z)
% Obtain latent variables for a new image
% Z = latent(model,F,J,ind,nit,Z)
% model - The model (computed by PatchCCAtrain and Patch2NN)
% F     - Pushed data to fit the latent variables to (3D/4D)
% J     - Pushed voxel counts
% ind   - Index of the data channel that F corresponds with
% nit   - Iteration settings
% Z     - Optional starting estimates for latent variables

if ~isfield(model,'P11'), error('Missing precisions.'); end

% Needs some error checking for dimensions etc
if nargin<6
    % Initial allocation of cell array of latent variables
    Z     = cell(size(model));
    for p=1:numel(model)
        if ~isempty(model(p).mod)
            Z{p} = zeros([size(model(p).mod(1).W,3),1],'single');
        end
    end
end

% Compute V once and attach to model
for p=1:numel(model)
    jp  = J(model(p).pos{:},ind);
    model(p).mod(ind).V  = computeV(model(p).P11, model(p).mod(ind).W, jp(:));
end

% Estimate most probable latent variables
for it=1:(2*nit(1)) % Black & White chessboard updates
    for p3=1:size(model,3) % Loop over z
        for p2=1:size(model,2) % Loop over y
            for p1=1:size(model,1) % Loop over x
                if ((rem(p1,2)==rem(p2,2))==rem(p3,2))==rem(it,2) % Is it a black/white patch
                    patch = model(p1,p2,p3);             % Get the patch
                    if ~isempty(patch.mod)               % Has it been modelled
                        dm  = cellfun(@numel,patch.pos); % Patch dimensions
                        Fp  = F(patch.pos{:},patch.c{ind}(1:(end-1))); % Get the image patch
                        Fp  = reshape(Fp,[prod(dm) size(Fp,4)]); % Convert 3D patches to column vectors

                        Jp  = J(patch.pos{:},:);         % Get the image patch
                        Jp  = reshape(Jp,[prod(dm) 1]);  % Convert 3D patches to column vectors

                        % Neighbouring latent variables (of other colour)
                        z2  = [getZ(p1  ,p2  ,p3+1, Z); getZ(p1  ,p2  ,p3-1, Z)
                               getZ(p1  ,p2+1,p3  , Z); getZ(p1  ,p2-1,p3  , Z)
                               getZ(p1+1,p2  ,p3  , Z); getZ(p1-1,p2  ,p3  , Z)];

                        % Estimate latent variables, conditional on neighbouring latent variables
                        Z{p1,p2,p3} = updateZ(Fp,Jp,patch.mod(ind).mu,patch.mod(ind).W,...
                                               patch.mod(ind).V,patch.P12,z2,Z{p1,p2,p3},nit(2));
                    end
                end
            end
        end
    end
end


function z = getZ(p1,p2,p3,Z)
if p1>=1 && p1<=size(Z,1) &&  p2>=1 && p2<=size(Z,2) && p3>=1 && p3<=size(Z,3) && ~isempty(Z{p1,p2,p3})
    z = Z{p1,p2,p3};
else
    z = zeros(0,1);
end

function V = computeV(P11, W, jp)
V  =  inv(P11 + hessZ(W,jp));


function H = hessZ(W,wt)
[Nvox,M,K] = size(W);
A    = Abohning(M);
H    = 0;
MK   = [M,K];
for i=1:Nvox
    Wi = reshape(W(i,:,:),MK);
    H  = H + Wi'*(wt(i)*A)*Wi;
end


function z = updateZ(Fp,jp,mu,W,V,P12,z2,z,nit)
g0 = -P12*z2;
for iter=1:nit
    z   = V*(g0 + numeratorZ(Fp,z,mu,W,jp));
end


function g = numeratorZ(F,z,mu,W,wt)
if isempty(W), g = zeros([size(W,3),1],'single'); return; end
[Nvox,M,K] = size(W);
A    = Abohning(M);
Psi0 = reshape(  reshape(W,[Nvox*M,K])*z,[Nvox,M]);
P    = softmax0(Psi0+mu,2);
r    = reshape(F+(Psi0*A-P).*wt,[1,Nvox*M]);
g    = reshape(r*reshape(W,[Nvox*M,K]),[K,1]);


function A = Abohning(M)
A  = 0.5*(eye(M)-1/(M+1));


function P = softmax0(Psi,d)
mx  = max(Psi,[],d);
E   = exp(Psi-mx);
P   = E./(sum(E,d)+exp(-mx));


%function P = softmax0Old(Psi,d)
%mx  = max(Psi,[],d);
%E   = exp(bsxfun(@minus,Psi,mx));
%P   = bsxfun(@rdivide, E, sum(E,d)+exp(-mx));


function r = classes(fl,mu)
% Responsibilities from likelihoods (fl) and priors (mu)
mu   = exp(bsxfun(@minus,mu,spm_mb_shape('lse',mu)));
sl   = sum(fl,4);
sp   = sum(mu,4);
r    = fl.*mu;
r    = bsxfun(@rdivide,r,sum(r,4)+(1-sl).*(1-sp));
r(~isfinite(r)) = 0;


function fl = likelihoods(datn,mbsett,mu1,select)
% Labelling using likelihoods, rather than responsibilities,
% is not yet very effective. If training also made use of the
% likelihood and priors reconstructed from the factorisation,
% then I would expect this to be rather more effective.
% 
if nargin<4, select=1:3; end
app    = spm_mb_appearance;
io     = spm_mb_io;
gmm    = datn.model.gmm;
fn     = io.get_image(gmm);
gmms   = mbsett.gmm(gmm.pop);
C      = gmms.C;
mg_ix  = gmms.mg_ix;
gam    = gmm.gam;
chan   = app.inu_basis(gmm.T,size(fn),datn.Mat,ones(1,C));
[~,mf,vf] = app.inu_recon(fn,[],chan,gmm.T,gmm.Sig);

% Adjust priors
mu1     = spm_mb_classes('template_k1',mu1,datn.delta); % Add implicit class
mu1(:,:,:,select) = 0;                       % Compute only likelihoods for foreground tissues
other  = setxor(select,1:size(mu,4));        % Tissues to be merged into the background class
mu2    = mu1(:,:,:,other);                   % Background tissue probabilities
mu1(:,:,:,other)  = bsxfun(@minus,mu2,lse1(mu2,4));

% Compute likelihoods
fl  = spm_gmmlib('resp',gmm.m,gmm.b,gmm.W,gmm.nu,gam,...
                 uint64(mg_ix), mu1,mf,vf, uint64(gmm.samp));
fl  = fl(:,:,:,select); % Select the foreground tissues


function l = lse1(mu,d)
mx = max(mu,[],d);
l  = log(sum(exp(bsxfun(@minus,mu,mx)),d)) + mx;


function P = reconmu(model,ind,Z)
% Reconstruct native-space data
% FORMAT Y = reconmu(model,ind,Z)
% model - The learned model
% ind   - Which view to reconstruct
% Z     - Cell array of latent variables

dm_temp = cellfun(@max,model(end,end,end).pos); % Size of template-space data
if nargin<2, ind = 1; end
c       = model(1).c{ind}(1:(end-1));
P       = zeros([dm_temp numel(c)],'single');
if nargin<3
    for p=1:numel(model) % Loop over patches
        patch             = model(p);
        dm_patch          = cellfun(@numel,patch.pos);            % Dimensions of patch
        psi               = reshape(patch.mod(ind).mu(:,c),[dm_patch numel(c)]);
        P(patch.pos{:},c) = psi;
    end
else
    for p=1:numel(model) % Loop over patches
        patch             = model(p);
        dm_patch          = cellfun(@numel,patch.pos);            % Dimensions of patch
        psi               = getPsi(Z{p},patch.mod(ind),dm_patch); % Linear combination of bases
        P(patch.pos{:},c) = psi;
    end
end
