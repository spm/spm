function out = spm_shoot_template(job)
% Iteratively compute a template with mean shape and intensities
% format spm_shoot_template(job)
% Fields of job:
%     job.images{1} first set of images (eg rc1*.nii)
%     job.images{2} second set of images (eg rc2*.nii)
%     etc
%
% Other settings are defined in spm_shoot_defaults.m
%
% The outputs are flow fields (v*.nii), deformation fields (y*.nii)
% and a series of Template images.
%_______________________________________________________________________
% Copyright (C) Wellcome Trust Centre for Neuroimaging (2009)

% John Ashburner
% $Id: spm_shoot_template.m 4324 2011-05-06 11:29:30Z john $

%_______________________________________________________________________
d       = spm_shoot_defaults;
tname   = d.tname;   % Base file name for templates
issym   = d.issym;   % Use a symmetric template

cyc_its = d.cyc_its; % No. multigrid cycles and inerations
smits   = d.smits;   % No. smoothing iterations
sched   = d.sched;   % Schedule for coarse to fine
nits    = numel(sched)-1;
rparam  = d.rparam;  % Regularisation parameters for deformation
sparam  = d.sparam;  % Regularisation parameters for blurring
eul_its = d.eul_its; % Start with fewer steps

bs_args = d.bs_args; % B-spline settings for interpolation
%_______________________________________________________________________


% Sort out handles to images
n1 = numel(job.images);
n2 = numel(job.images{1});
NF = struct('NI',[],'vn',[1 1]);
NF(n1,n2) = struct('NI',[],'vn',[1 1]);

% Pick out individual volumes within NIfTI files
for i=1:n1,
    if numel(job.images{i}) ~= n2,
        error('Incompatible number of images');
    end;
    for j=1:n2,
        [pth,nam,ext,num] = spm_fileparts(job.images{i}{j});
        NF(i,j).NI        = nifti(fullfile(pth,[nam ext]));
        num               = [str2num(num) 1 1];
        NF(i,j).vn        = num(1:2);
    end;
end;

spm_progress_bar('Init',n2,'Initial mean','Subjects done');
dm = [size(NF(1,1).NI.dat) 1];
dm = dm(1:3);

NU     = nifti;
NU(n2) = nifti;
NY     = nifti;
NY(n2) = nifti;
NJ     = nifti;
NJ(n2) = nifti;

t  = zeros([dm n1+1],'single');
for i=1:n2,
    % Generate files for flow fields, deformations and Jacobian determinants.
    [pth,nam,ext]   = fileparts(NF(1,i).NI.dat.fname);

    NU(i) = nifti;
    NY(i) = nifti;
    NJ(i) = nifti;

    offs  = 352;
    if ~isempty(tname),
        NU(i).dat = file_array(fullfile(pth,['v_' nam '_' tname '.nii']),...
                               [dm 1 3], 'float32-le', 352, 1, 0);
        NY(i).dat = file_array(fullfile(pth,['y_' nam '_' tname '.nii']),...
                               [dm 1 3], 'float32-le', 352, 1, 0);
        NJ(i).dat = file_array(fullfile(pth,['j_' nam '_' tname '.nii']),...
                               [dm 1 1], 'float32-le', 352, 1, 0);
    else
        NU(i).dat = file_array(fullfile(pth,['v_' nam '.nii']),...
                               [dm 1 3], 'float32-le', 352, 1, 0);
        NY(i).dat = file_array(fullfile(pth,['y_' nam '.nii']),...
                               [dm 1 3], 'float32-le', 352, 1, 0);
        NJ(i).dat = file_array(fullfile(pth,['y_' nam '.nii']),...
                               [dm 1 1], 'float32-le', 352, 1, 0);
    end

    NU(i).descrip = sprintf('Velocity (%d %.4g %.4g %.4g)', rparam(1), rparam(2), rparam(3), rparam(4));
    NU(i).mat     = NF(1,i).NI.mat;
    NU(i).mat0    = NF(1,i).NI.mat0;

    NY(i).descrip = 'Deformation (templ. to. ind.)';
    NY(i).mat     = NF(1,i).NI.mat;

    NJ(i).descrip = 'Jacobian det (templ. to. ind.)';
    NJ(i).mat     = NF(1,i).NI.mat;


    create(NU(i)); NU(i).dat(:,:,:,:,:) = 0;
    create(NY(i)); NY(i).dat(:,:,:,:,:) = reshape(affind(dartel3('Exp',zeros([dm,3],'single'),[0 1]),NU(i).mat0),[dm,1,3]);
    create(NJ(i)); NJ(i).dat(:,:,:)     = 1;

    % Add to sufficient statistics for generating initial template
    for j=1:n1,
        vn         = NF(j,i).vn;
        dat        = NF(j,i).NI.dat(:,:,:,vn(1),vn(2));
        msk        = isfinite(dat);
        dat(~msk)  = 0;
        t(:,:,:,j) = t(:,:,:,j) + dat;
    end;
    t(:,:,:,end) = t(:,:,:,end) + msk;
    clear tmp msk
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');

% Make symmetric (if necessary)
if issym, t = t + t(end:-1:1,:,:,:); end

% Generate template from sufficient statistics
tmp = t(:,:,:,end);
msk = tmp<=0;
for j=1:n1,
    tmp = tmp - t(:,:,:,j);
end
tmp(msk) = 0.01; % Most NaNs are likely to be background
t(:,:,:,end) = tmp;
clear tmp msk
M  = NF(1,1).NI.mat;
vx = sqrt(sum(M(1:3,1:3).^2));
t  = max(t,0);
g  = cell(n1+1,1);

%%%%%%%%%%%%%%%
% save SuffStats0.mat t sparam vx sched smits
%%%%%%%%%%%%%%%

if ~isempty(sparam),
    g0 = spm_shoot_blur(t,[sparam(1), vx, sched(1)*sparam(2) sparam(3:4)],smits);
    for j=1:n1+1,
        g{j} = max(g0(:,:,:,j),1e-4);
    end
    clear g0
else
    sumt = max(sum(t,4),0)+eps;
    for j=1:n1+1,
        g{j} = (t(:,:,:,j)+0.1)./(sumt+0.1*(n1+1));
    end
    clear sumt
end

if ~isempty(tname),
    % Write template
    NG = NF(1,1).NI;
    NG.descrip       = sprintf('Avg of %d', n2);
    [tdir,nam,ext]   = fileparts(job.images{1}{1});
    NG.dat.fname     = fullfile(tdir,[tname, '_0.nii']);
    NG.dat.dim       = [dm n1+1];
    NG.dat.dtype     = 'float32-le';
    NG.dat.scl_slope = 1;
    NG.dat.scl_inter = 0;
    NG.mat0          = NG.mat;
    create(NG);
    for j=1:n1+1,
        NG.dat(:,:,:,j)  = g{j};
    end
end
for j=1:n1+1, g{j} = spm_bsplinc(log(g{j}), bs_args); end

vx = sqrt(sum(NG.mat(1:3,1:3).^2));

ok = true(n2,1);

% The actual work
for it=1:nits,

    % More regularisation in the early iterations, as well as a
    % a less accurate approximation in the integration.
    prm      = [rparam(1), vx, rparam(2:3)*sched(it+1), rparam(4)];
    int_args = [eul_its(it), cyc_its];
    drawnow

    t   = zeros([dm n1+1],'single');
    su  = zeros([dm 3]);

    % Update velocities
    spm_progress_bar('Init',n2,sprintf('Update velocities (%d)',it),'Subjects done');
    for i=1:n2, % Loop over subjects

        if ok(i)
            fprintf(' %-3d %-5d\t| ',it,i);

            % Load image data for this subject
            f = loadimage(NF(:,i));

            % Load this subject's flow field and deformation
            u  = squeeze(single(NU(i).dat(:,:,:,:,:)));
            y  = affind(squeeze(single(NY(i).dat(:,:,:,:,:))),inv(NU(i).mat0));
            dt = squeeze(single(NJ(i).dat(:,:,:)));
            drawnow

            % Gauss-Newton iteration to re-estimate deformations for this subject
            u = spm_shoot_update(g,f,u,y,dt,prm,int_args,bs_args);
            clear f y

            drawnow
            NU(i).dat(:,:,:,:,:) = reshape(u,[dm 1 3]);
            su = su + u;
            clear u
            spm_progress_bar('Set',i);
        end

    end
    spm_progress_bar('Clear');

    if issym
        su(:,:,:,1)   = (su(:,:,:,1)   - su(end:-1:1,:,:,1)  )/(sum(ok)*2);
        su(:,:,:,2:3) = (su(:,:,:,2:3) + su(end:-1:1,:,:,2:3))/(sum(ok)*2);
    else
        su = su/sum(ok);
    end
    
    % Update template sufficient statistics
    spm_progress_bar('Init',n2,sprintf('Update deformations and template (%d)',it),'Subjects done');
    for i=1:n2, % Loop over subjects

        if ok(i)
            % Load velocity, mean adjust and re-save
            u  = squeeze(single(NU(i).dat(:,:,:,:,:)));
            u  = u - su; % Subtract mean
            NU(i).dat(:,:,:,:,:) = reshape(u,[dm 1 3]);

            % Generate inverse deformation and save
            [y,J] = spm_shoot3d(u,prm,int_args);
            clear u
            dt    = dartel3('det',J); clear J

            if any(~isfinite(dt(:)) | dt(:)>100 | dt(:)<1/100)
                ok(i) = false;
                fprintf('Problem with %s (dets: %g .. %g)\n', NU(i).dat.fname, min(dt(:)), max(dt(:)));
                clear dt
            end
        end

        if ok(i)
            NY(i).dat(:,:,:,:,:) = reshape(affind(y,NU(i).mat0),[dm 1 3]);
            NJ(i).dat(:,:,:)     = dt;
            drawnow;

            % Load image data for this subject
            f = loadimage(NF(:,i));

            % Increment sufficient statistic for template
            for j=1:n1+1,
                tmp        = dartel3('samp',f{j},y);
                t(:,:,:,j) = t(:,:,:,j) + tmp.*dt;
                drawnow
            end;

            clear f y dt
            fprintf('.');
            spm_progress_bar('Set',i);
        end

    end
    clear su
    fprintf('\n');
    spm_progress_bar('Clear');

    % Make symmetric (if necessary)
    if issym, t = t + t(end:-1:1,:,:,:); end

    %%%%%%%%%%%%%%%
    % save(['SuffStats' num2str(it) '.mat'],'t', 'sparam', 'vx', 'sched', 'smits', 'g');
    %%%%%%%%%%%%%%%

    % Re-generate template data from sufficient statistics
    if ~isempty(sparam),
        g0 = reconv(g,bs_args);
        g0 = spm_shoot_blur(t,[sparam(1), vx, sched(it+1)*sparam(2) sparam(3:4)],smits,g0);
        g  = cell(n1+1,1);
        for j=1:n1+1,
            g{j} = max(g0(:,:,:,j),1e-4);
        end
        clear g0
    else
        sumt = max(sum(t,4),0)+eps;
        for j=1:n1+1,
            g{j} = (t(:,:,:,j)+0.1)./(sumt+0.1*(n1+1));
        end
        clear sumt
    end
    clear t

    % Write template
    if ~isempty(tname),
        NG.dat.fname    = fullfile(tdir,[tname '_' num2str(ceil(it/6)) '.nii']);
        create(NG);
        for j=1:n1+1,
            NG.dat(:,:,:,j) = g{j};
        end
    end

    % Compute template's B-spline coefficients
    for j=1:n1+1, g{j} = spm_bsplinc(log(g{j}), bs_args); end
    drawnow
end

if any(~ok)
    fprintf('Problems with:\n');
    for i=find(~ok)',
        fprintf('\t%s\n', NU(i).dat.fname);
    end
end

% Finish off
[tdir,nam,ext] = fileparts(job.images{1}{1});
out.template = cell(floor(nits/3)+1,1);
if ~isempty(tname),
    for it=0:floor(nits/3),
        fname    = fullfile(tdir,[tname '_' num2str(it) '.nii']);
        out.template{it+1} = fname;
    end
end
out.vel = cell(n2,1);
out.def = cell(n2,1);
out.jac = cell(n2,1);
for i=1:n2,
    out.vel{i} = NU(i).dat.fname;
    out.def{i} = NY(i).dat.fname;
    out.jac{i} = NJ(i).dat.fname;
end
%=======================================================================

%=======================================================================
function y1 = affind(y0,M)
% Affine transform of deformation
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%=======================================================================

%=======================================================================
function g0 = reconv(g,bs_args)
d = [size(g{1}), 1];
[i1,i2,i3]=ndgrid(1:d(1),1:d(2),1:d(3));
g0 = zeros([d,numel(g)],'single');
for k=1:numel(g),
    g0(:,:,:,k) = max(exp(spm_bsplins(g{k},i1,i2,i3,bs_args)),1e-4);
end
%=======================================================================

%=======================================================================
function f = loadimage(NF)
n1      = size(NF,1);
f       = cell(n1+1,1);
dm      = [NF(1).NI.dat.dim 1 1 1];
dm      = dm(1:3);
f{n1+1} = ones(dm,'single');
for j=1:n1,
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

