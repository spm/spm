function spm_dartel_logodds(job)
% Generate LogOdds files
% FORMAT spm_dartel_logodds(job)
% job.flowfields
% job.images
% job.K
% job.rparam
%
% The aim is to obtain better pattern recognition through using
% LogOdds.  See Killian Pohl's work for more info.
%_______________________________________________________________________
% Copyright (C) 2007 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_dartel_logodds.m 964 2007-10-19 16:35:34Z john $


PI = job.images;
PU = job.flowfields;
K  = job.K;
lam= job.rparam;
n  = numel(PI);


spm_progress_bar('Init',numel(PI{1}),'Creating LogOdds');

for i=1:numel(PI{1}),

    % Derive deformations and Jacobian determinants from the
    % flow fields.
    NU = nifti(PU{i});
    [pth,nam,ext,num] = spm_fileparts(NU.dat.fname);
    fprintf('%s: def ',nam); drawnow;
    u      = single(squeeze(NU.dat(:,:,:,1,:)));
    [y,dt] = dartel3('Exp',u,[K -1 1]);
    y1 = double(y(:,:,:,1));
    y2 = double(y(:,:,:,2));
    y3 = double(y(:,:,:,3));
    clear y
    dt = max(dt,single(0));

    f  = zeros([NU.dat.dim(1:3),n+1],'single');
 
    [pth,nam,ext,num] = spm_fileparts(PI{1}{i});
    NI = nifti(fullfile(pth,[nam ext]));
    NO = NI;
    NO.dat.fname     = fullfile('.',['logodds_' nam '.nii']);
    NO.dat.scl_slope = 1.0;
    NO.dat.scl_inter = 0.0;
    NO.dat.dtype     = 'float32-le';
    NO.dat.dim       = [NU.dat.dim(1:3) 1 n+1];
    NO.mat           = NU.mat;
    NO.mat0          = NU.mat;
    NO.mat_intent    = 'Aligned';
    NO.mat0_intent   = 'Aligned';
    NO.descrip = 'LogOdds';
    create(NO);
    vx = sqrt(sum(NU.mat(1:3,1:3).^2));

    % Compute the warped tissue probabilities
    f(:,:,:,end) = 1; % Background
    for j=1:n,
        fprintf('%d ', j); drawnow;

        [pth,nam,ext,num] = spm_fileparts(PI{j}{i});
        NI = nifti(fullfile(pth,[nam ext]));
        if sum(sum((NI.mat  - NU.mat ).^2)) < 0.0001 && ...
           sum(sum((NI.mat0 - NU.mat0).^2)) < 0.0001,
            ty1 = y1;
            ty2 = y2;
            ty3 = y3;
        else
            mat = NI.mat;
            if isfield(NI,'extras') && isfield(NI.extras,'mat'),
                mat1 = NI.extras.mat;
                if size(mat1,3) >= j && sum(sum(mat1(:,:,j).^2)) ~=0,
                    mat = mat1;
                end;
            end;
            M   = mat\NU.mat0;
            ty1 = M(1,1)*y1 + M(1,2)*y2 + M(1,3)*y3 + M(1,4);
            ty2 = M(2,1)*y1 + M(2,2)*y2 + M(2,3)*y3 + M(2,4);
            ty3 = M(3,1)*y1 + M(3,2)*y2 + M(3,3)*y3 + M(3,4);
        end;
        spl_param  = [3 3 3  1 1 1];
        cf         = spm_bsplinc(NI.dat(:,:,:,1,1),spl_param);
        f(:,:,:,j) = spm_bsplins(cf,ty1,ty2,ty3,spl_param);
        f(:,:,:,end) = f(:,:,:,end) - f(:,:,:,j);
        clear cf ty1 ty2 ty3
    end;
    clear y1 y2 y3

    a = logodds(f,dt,lam,12,vx);
    a = softmax(a);
    fprintf('\n');

    clear f dt
    for j=1:n+1,
        NO.dat(:,:,:,1,j) = log(a(:,:,:,j));
    end
    clear a
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


function a = logodds(t,s,lam,its,vx)
if nargin<5, vx  = [1 1 1]; end;
if nargin<4, its = 12;      end;
if nargin<3, lam = 1;       end;

a   = zeros(size(t),'single');
d   = size(t);
W   = zeros([d(1:3) round((d(4)*(d(4)+1))/2)],'single');
for i=1:its,
    sig = softmax(a);
    gr  = sig - t;
    for k=1:d(4), gr(:,:,:,k) = gr(:,:,:,k).*s; end
    gr  = gr + optimN('vel2mom',a,[2 vx lam 0 0]);
    jj  = d(4)+1;
    for j1=1:d(4),
        W(:,:,:,j1) = sig(:,:,:,j1).*(1-sig(:,:,:,j1)).*s;
        for j2=1:j1-1,
            W(:,:,:,jj) = -sig(:,:,:,j1).*sig(:,:,:,j2).*s;
            jj = jj+1;
        end
    end
    a   = a - optimN(W,gr,[2 vx lam 0 lam*1e-3 1 1]);
    %a  = a - mean(a(:)); % unstable
    a   = a - sum(sum(sum(sum(a))))/numel(a);
    fprintf('.');
    drawnow 
end;
return;

function sig = softmax(a)
sig = zeros(size(a),'single');
for j=1:size(a,3),
    aj = double(squeeze(a(:,:,j,:)));
    aj = min(max(aj-mean(aj(:)),-700),700);
    sj = exp(aj);
    s  = sum(sj,3);
    for i=1:size(a,4),
        sig(:,:,j,i) = single(sj(:,:,i)./s);
    end
end;

