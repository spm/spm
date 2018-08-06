function out = spm_shoot_warp(job)
% Register images with template
% format spm_shoot_warp(job)
% Fields of job:
%     job.images{1} first set of images (eg rc1*.nii)
%     job.images{2} second set of images (eg rc2*.nii)
%     etc
%     job.templates template files
% Other settings are defined in spm_shoot_defaults.m
%
% The outputs are flow fields (v*.nii), deformation fields (y*.nii) and
% Jacobian determinants (j*.nii)
%_______________________________________________________________________
% Copyright (C) Wellcome Trust Centre for Neuroimaging (2009)

% John Ashburner
% $Id: spm_shoot_warp.m 7389 2018-08-06 13:35:48Z john $

%_______________________________________________________________________
d       = spm_shoot_defaults;
tname   = d.tname;   % Base file name for templates

cyc_its = d.cyc_its; % No. multigrid cycles and inerations
sched   = d.sched;   % Schedule for coarse to fine
nits    = numel(sched)-1;
rparam  = d.rparam;  % Regularisation parameters for deformation
eul_its = d.eul_its; % Start with fewer steps
scale   = d.scale;   % Fraction of Gauss-Newton update step to use

bs_args = d.bs_args; % B-spline settings for interpolation
%_______________________________________________________________________

spm_diffeo('boundary',0); % Set boundary condition

% Sort out handles to images
n1 = numel(job.images);
n2 = numel(job.images{1});
NF = struct('NI',[],'vn',[1 1]);
NF(n1,n2) = struct('NI',[],'vn',[1 1]);

% Pick out individual volumes within NIfTI files
for i=1:n1
    if numel(job.images{i}) ~= n2
        error('Incompatible number of images');
    end
    for j=1:n2
        [pth,nam,ext,num] = spm_fileparts(job.images{i}{j});
        NF(i,j).NI        = nifti(fullfile(pth,[nam ext]));
        num               = [str2num(num) 1 1];
        NF(i,j).vn        = num(1:2);
    end
end

dm = [size(NF(1,1).NI.dat) 1];
dm = dm(1:3);

% Sort out which template for each iteration
tmpl_no = round(((1:nits)-1)/(nits-1)*(numel(job.templates)-0.51))+1;

ok = true(n2,1);

NU     = nifti;
NU(n2) = nifti;
NY     = nifti;
NY(n2) = nifti;
NJ     = nifti;
NJ(n2) = nifti;

for i=1:n2 % Loop over subjects

    % Generate files for flow fields, deformations and Jacobian determinants.
    [pth,nam,ext]   = fileparts(NF(1,i).NI.dat.fname);
    if ~isempty(tname), nam = [nam '_' tname]; end
    offs  = 352;

    NU(i) = nifti;
    NU(i).dat = file_array(fullfile(pth,['v_' nam '.nii']),[dm 1 3], 'float32-le', offs, 1, 0);
    NU(i).descrip = sprintf('Velocity (%d %.4g %.4g %.4g)', rparam(1), rparam(2), rparam(3), rparam(4));
    NU(i).mat     = NF(1,i).NI.mat;
    NU(i).mat0    = NF(1,i).NI.mat0;
    create(NU(i)); NU(i).dat(:,:,:,:,:) = 0;

    NY(i) = nifti;
    NY(i).dat = file_array(fullfile(pth,['y_' nam '.nii']),[dm 1 3], 'float32-le', offs, 1, 0);
    NY(i).descrip = 'Deformation (templ. to. ind.)';
    NY(i).mat     = NF(1,i).NI.mat;
    create(NY(i)); NY(i).dat(:,:,:,:,:) = reshape(affind(spm_diffeo('Exp',zeros([dm,3],'single'),[0 1]),NU(i).mat0),[dm,1,3]);

    NJ(i) = nifti;
    NJ(i).dat = file_array(fullfile(pth,['j_' nam '.nii']),[dm 1 1], 'float32-le', offs, 1, 0);
    NJ(i).descrip = 'Jacobian det (templ. to. ind.)';
    NJ(i).mat     = NF(1,i).NI.mat;
    create(NJ(i)); NJ(i).dat(:,:,:)     = 1;

    drawnow
end

for i=1:n2 % Loop over subjects. Can replace FOR with PARFOR.

    % Load image data for this subject
    f  = loadimage(NF(:,i));
    u  = squeeze(single(NU(i).dat(:,:,:,:,:)));
    y  = affind(squeeze(single(NY(i).dat(:,:,:,:,:))),inv(NU(i).mat0));
    dt = squeeze(single(NJ(i).dat(:,:,:)));

    % Re-load first template
    [g,vx] = load_template(job.templates{tmpl_no(1)}, n1, bs_args);

    % The actual work
    for it=1:nits

        % Load template appropriate for this iteration
        if (it>1) && (tmpl_no(it)~=tmpl_no(it-1))
            [g,vx] = load_template(job.templates{tmpl_no(it)}, n1, bs_args);
        end

        % More regularisation in the early iterations, as well as a
        % a less accurate approximation in the integration.
        prm      = [vx, rparam*sched(it+1)*prod(vx)];
        int_args = [eul_its(it), cyc_its];
        drawnow

        fprintf(' %-5d %-3d\t| ',i,it);

        % Gauss-Newton iteration to re-estimate deformations for this subject
        u      = spm_shoot_update(g,f,u,y,dt,prm,bs_args,scale); drawnow
        [y,dt] = defdet(u,prm,int_args);

        if any(~isfinite(dt(:)) | dt(:)>100 | dt(:)<1/100)
            ok(i) = false;
            fprintf('Problem with %s (dets: %g .. %g)\n', NU(i).dat.fname, min(dt(:)), max(dt(:)));
           %clear dt
            break
        end

        drawnow
        NU(i).dat(:,:,:,:,:) = reshape(u,[dm 1 3]);
        NY(i).dat(:,:,:,:,:) = reshape(affind(y,NU(i).mat0),[dm 1 3]);
        NJ(i).dat(:,:,:)     = dt;

    end
    fprintf('\n');

    drawnow
end

if any(~ok)
    fprintf('Problems with:\n');
    for i=find(~ok)'
        fprintf('\t%s\n', NU(i).dat.fname);
    end
end

% Finish off
out.vel = cell(n2,1);
out.def = cell(n2,1);
out.jac = cell(n2,1);
for i=1:n2
    out.vel{i} = NU(i).dat.fname;
    out.def{i} = NY(i).dat.fname;
    out.jac{i} = NJ(i).dat.fname;
end
%=======================================================================

%=======================================================================
function y1 = affind(y0,M)
% Affine transform of deformation
y1 = zeros(size(y0),'single');
for d=1:3
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%=======================================================================

%=======================================================================
function f = loadimage(NF)
n1      = size(NF,1);
f       = cell(n1+1,1);
dm      = [NF(1).NI.dat.dim 1 1 1];
dm      = dm(1:3);
f{n1+1} = ones(dm,'single');
for j=1:n1
    vn      = NF(j,1).vn;
    f{j}    = single(NF(j,1).NI.dat(:,:,:,vn(1),vn(2)));
    msk     = ~isfinite(f{j});
    f{j}(msk) = 0;
    f{n1+1} = f{n1+1} - f{j};
    drawnow
end
f{n1+1}(msk) = 0.00001;
%=======================================================================

%=======================================================================
function [g,vx] = load_template(template, n1, bs_args)
g  = cell(n1+1,1);
NG = nifti(template);
if size(NG.dat,4) < n1+1
    error('Not enough tissues in template (%d < %d+1).', size(NG.dat,4),n1);
end

bg = ones([size(NG.dat,1), size(NG.dat,2), size(NG.dat,3)]);
for j=1:n1
    tmp  = NG.dat(:,:,:,j);
    g{j} = spm_bsplinc(log(tmp), bs_args);
    bg   = bg - tmp;
    clear tmp;
end
g{n1+1}  = log(max(bg,eps));
clear bg

vx = sqrt(sum(NG.mat(1:3,1:3).^2));
%=======================================================================

%=======================================================================
function [y,dt] = defdet(u,prm,int_args)
% Generate deformation
[y,J] = spm_shoot3d(u,prm,int_args);
dt    = spm_diffeo('det',J);
%=======================================================================

%=======================================================================

