function out = spm_dartel_template(job)
% Iteratively compute a template with mean shape and intensities
% format spm_dartel_template(job)
%
% The outputs are flow fields, and a series of Template images.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_template2.m 2231 2008-09-29 15:01:05Z guillaume $

Ntime = 4;

code = 2;
st = job.settings;
K  = st.param(1).K;
n1 = numel(job.images);
n2 = numel(job.images{1});
NF = struct('NI',[],'vn',[1 1]);
NF(n1,n2) = struct('NI',[],'vn',[1 1]);

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

spm_progress_bar('Init',n2,'Initial mean','Images done');
dm = [size(NF(1,1).NI.dat) 1];
dm = dm(1:3);
NU = cat(2,NF(1,:).NI);
t  = zeros([dm n1+1],'single');
tname     = deblank(job.settings.template);
out.files = cell(n2,1);
for i=1:n2,
    [pth,nam,ext]   = fileparts(NU(i).dat.fname);
    if ~isempty(tname),
        NU(i).dat.fname = fullfile(pth,['u_' nam '_' tname '.nii']);
    else
        NU(i).dat.fname = fullfile(pth,['u_' nam '.nii']);
    end
    NU(i).dat.dim   = [dm Ntime 3];
    NU(i).dat.dtype = 'float32-le';
    NU(i).dat.scl_slope = 1;
    NU(i).dat.scl_inter = 0;
    NU(i).descrip = 'Flow Field';

    if exist(NU(i).dat.fname,'file'),
        ntmp = nifti(NU(i).dat.fname);
        if size(ntmp.dat,4) ~= Ntime,
            nt0  = size(ntmp.dat,4);
            nt1  = size(NU(i).dat,4);
            scal = nt1/nt0;
            ind  = max(min(round(((1:nt1)-0.5)/nt1*nt0+0.5),nt0),1);
            uo   = single(ntmp.dat(:,:,:,:,:));
            for j=1:nt1,
                NU(i).dat(:,:,:,j,:) = uo(:,:,:,ind(j),:)*scal;
            end
        end
        y   = spm_dartel_integrate(NU(i).dat,[0 1],K);
        tmp = cell(1);
        for j=1:n1,
            vn = NF(j,i).vn;
            if j==n1, tmp=cell(1,2); end
            [tmp{:}] = dartel3('push',...
                single(NF(j,i).NI.dat(:,:,:,vn(1),vn(2))),y);
            t(:,:,:,j) = t(:,:,:,j) + tmp{1};
        end
        t(:,:,:,end) = t(:,:,:,end) + tmp{2};
        clear y tmp
    else
        create(NU(i));
        NU(i).dat(:,:,:,:,:) = 0;
        for j=1:n1,
            vn         = NF(j,i).vn;
            dat        = NF(j,i).NI.dat(:,:,:,vn(1),vn(2));
            msk        = isfinite(dat);
            dat(~msk)  = 0;
            t(:,:,:,j) = t(:,:,:,j) + dat;
        end;
        t(:,:,:,end) = t(:,:,:,end) + msk;
        clear tmp msk
    end;
    spm_progress_bar('Set',i);
    out.files{i} = NU(i).dat.fname;
end;
spm_progress_bar('Clear');

M  = NF(1,1).NI.mat;
vx = sqrt(sum(M(1:3,1:3).^2));
if st.param(1).slam,
    for j=1:n1,
        t(:,:,:,end) = t(:,:,:,end) - t(:,:,:,j);
    end
    t = max(t,0);
    g = spm_dartel_smooth(t,st.param(1).slam*2,8,vx);
else
    for j=1:n1,
        g(:,:,:,j) = t(:,:,:,j)./(t(:,:,:,end)+eps);
    end
end

if ~isempty(tname),
    NG = NF(1,1).NI;
    NG.descrip       = sprintf('Avg of %d', n2);
    [tdir,nam,ext]   = fileparts(job.images{1}{1});
    NG.dat.fname     = fullfile(tdir,[tname, '_0.nii']);
    NG.dat.dim       = [dm n1];
    NG.dat.dtype     = 'float32-le';
    NG.dat.scl_slope = 1;
    NG.dat.scl_inter = 0;
    NG.mat0          = NG.mat;
    create(NG);
    NG.dat(:,:,:,:)  = g(:,:,:,1:n1);
    out.template     = cell(numel(job.settings.param)+1,1);
    out.template{1}  = NG.dat.fname;
end

it0 = 0;
for it=1:numel(st.param),
    param = st.param(it);
    prm   = [st.rform, param.rparam, st.optim.lmreg, ...
             st.optim.cyc, st.optim.its, param.K, code];
    drawnow

    for j=1:param.its,
        it0 = it0 + 1;
        t   = zeros([dm n1+1],'single');

        for i=1:n2,
            f = zeros([dm n1],'single');
            for j=1:n1,
                vn         = NF(j,i).vn;
                f(:,:,:,j) = single(NF(j,i).NI.dat(:,:,:,vn(1),vn(2)));
                drawnow
            end;
            [NU(i).dat,ll,f1] = spm_dartel_variable(NU(i).dat,f,g,prm);
            fprintf('%d %d\t%g\t%g\t%g\t%g\n',it0,i,ll(1),ll(2),ll(1)+ll(2),ll(3));
            t      = t + f1;
        end;
        if param.slam,
            for j=1:n1,
                t(:,:,:,end) = t(:,:,:,end) - t(:,:,:,j);
            end
            t(:,:,:,end) = max(t(:,:,:,end),0);
            g = spm_dartel_smooth(t,param.slam,8,vx);
        else
            for j=1:n1,
                g(:,:,:,j) = t(:,:,:,j)./(t(:,:,:,end)+eps);
            end
        end
        clear t
        if ~isempty(tname),
            NG.dat.fname       = fullfile(tdir,[tname '_' num2str(it) '.nii']);
            create(NG);
            NG.dat(:,:,:,:)    = g(:,:,:,1:n1);
            out.template{it+1} = NG.dat.fname;
        end
        drawnow
    end;
end;

