function spm_dartel_template(job)
% Iteratively compute a template with mean shape and intensities
% format spm_dartel_template(job)
%
% The outputs are flow fields, and a series of Template images.
%_______________________________________________________________________
% Copyright (C) 2007 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$

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
g  = zeros([dm n1],'single');
nd = zeros(dm(1:3),'single');
for i=1:n2,
    [pth,nam,ext]   = fileparts(NU(i).dat.fname);
    NU(i).dat.fname = fullfile(pth,['u_' nam '.nii']);
    NU(i).dat.dim   = [dm 1 3];
    NU(i).dat.dtype = 'float32-le';
    NU(i).dat.scl_slope = 1;
    NU(i).dat.scl_inter = 0;
    NU(i).descrip = 'Flow Field';

    vn = NF(1,i).vn;
    tmp = find(~finite(NF(1,i).NI.dat(:,:,:,vn(1),vn(2))));
    if ~isempty(tmp),
        for j=1:n2,
            vn  = NF(j,i).vn;
            dat = NF(j,i).NI.dat(:,:,:,vn(1),vn(2));
            dat(tmp)    = 0;
            NF(j,i).NI.dat = dat;
            clear dat
        end;
    end;
    if exist(NU(i).dat.fname,'file'),
        u = NU(i).dat(:,:,:,1,:);
        u = single(squeeze(u));
        [y,dt] = dartel3('Exp',u,[K -1 1]);
        dt = max(dt,0);
        clear u
        for j=1:n1,
            vn         = NF(j,i).vn;
            g(:,:,:,j) = g(:,:,:,j) + dt.*dartel3('samp',single(NF(j,i).NI.dat(:,:,:,vn(1),vn(2))),y);
        end;
        nd = nd + dt;
        clear y dt
    else
        create(NU(i));
        NU(i).dat(:,:,:,:,:) = 0;
        for j=1:n1,
            vn         = NF(j,i).vn;
            g(:,:,:,j) = g(:,:,:,j) + NF(j,i).NI.dat(:,:,:,vn(1),vn(2));
        end;
        nd = nd + 1;
    end;
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');

for j=1:n1,
    g(:,:,:,j) = g(:,:,:,j)./nd;
end;

NG = NF(1,1).NI;
NG.descrip       = sprintf('Avg of %d', n2);
[tdir,nam,ext]   = fileparts(job.images{1}{1});
NG.dat.fname     = fullfile(tdir,'Template0.nii');
NG.dat.dim       = [dm n1];
NG.dat.dtype     = 'float32-le';
NG.dat.scl_slope = 1;
NG.dat.scl_inter = 0;
NG.mat0          = NG.mat;
create(NG);
NG.dat(:,:,:,:)  = g;

it0 = 0;
for it=1:numel(st.param),
    param = st.param(it);
    prm   = [st.rform, param.rparam, st.optim.lmreg, ...
             st.optim.cyc, st.optim.its, param.K, st.code];
    drawnow

    for j=1:param.its,
        it0 = it0 + 1;
        ng = zeros(size(g),'single');
        nd = zeros([size(g,1),size(g,2),size(g,3)],'single');

        for i=1:n2,
            f = zeros([dm n1],'single');
            for j=1:n1,
                vn         = NF(j,i).vn;
                f(:,:,:,j) = single(NF(j,i).NI.dat(:,:,:,vn(1),vn(2)));
                drawnow
            end;
            u = squeeze(single(NU(i).dat(:,:,:,:,:)));
            drawnow
            [u,ll] = dartel3(u,f,g,prm);
            fprintf('%d %d\t%g\t%g\t%g\t%g\n',it0,i,ll(1),ll(2),ll(1)+ll(2),ll(3));
            drawnow

            NU(i).dat(:,:,:,:,:) = reshape(u,[dm 1 3]);
            [y,dt] = dartel3('Exp',u,[param.K -1 1]);

            dt = max(dt,0);
            clear u
            drawnow;
            for j=1:n1,
                ng(:,:,:,j) = ng(:,:,:,j) + dartel3('samp',f(:,:,:,j),y).*dt;
                drawnow
            end;
            nd = nd + dt;
            clear y dt
        end;
        for j=1:n1,
            g(:,:,:,j)  = single(ng(:,:,:,j)./nd);
        end;
        clear ng nd
        NG.dat.fname    = fullfile(tdir,['Template' num2str(it) '.nii']);
        create(NG);
        NG.dat(:,:,:,:) = g;
        drawnow
    end;
end;

