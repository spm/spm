function fil_push_train_data(dw, Mw, Niiy, Nii1)
% Generate ``modulated categorical data'' for fil training
% FORMAT fil_push_train_data(dw, Mw, Niiy, Nii1)
% dw   - image dimensions of output
% Mw   - voxel-to-world mapping of output
% Niiy - NIfTI data structure of deformations
% Nii1 - NIfTI data structure of categorical image data to push
%        Note that the first dimension encodes the number of subjects
%        and the behaviour of the code depends on the second dimension.
%        * If the second dimension is 1, then the images are assumed
%          to be categorical labels. The output is a pcat_blah.mat file
%          containing a sparse matrix that encodes the pushed labels.
%          warped labels.
%        * If the second dimension is greater than 1, then the images
%          are assumed to encode segmentation probabilities. The output
%          in this case is a 4D image file. Note that the total
%          Number of categories is the number of dimensions + 1, accounting
%          for an implicit background class.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


assert(isnumeric(dw) && isreal(dw) && numel(dw)==3);
assert(ismatrix(Mw)  && isreal(Mw) && all(size(Mw)==[4 4]));
assert(isa(Niiy,'nifti'));
N = size(Niiy,1);
assert(isa(Nii1,'nifti') && size(Nii1,1)==N);
M = size(Nii1,2);

for n=1:N
    Ny   = Niiy(n);
    Nc   = Nii1(n,:);
    [~,nam,~] = fileparts(Nc(1).dat.fname);
    if M==1
        fname = fullfile('.',['pcat_' nam '.mat']);
        write_pushed_label(dw,Mw,...
                           single(squeeze(Ny.dat(:,:,:,:,:))),Ny.mat,...
                           Nc.dat(:,:,:),Nc(1).mat,...
                           fname);
    else
        fname = fullfile('.',['pcat_' nam '.nii']);
        fa = cat(4,Nc.dat);
        write_pushed_seg(dw,Mw,...
                         single(squeeze(Ny.dat(:,:,:,:,:))),Ny.mat,...
                         fa(:,:,:,:),Nc(1).mat,...
                         fname);
    end
end


function write_pushed_label(dw,Mw,psi1,My,X0,Mf,fname)
sh   = spm_mb_shape;
df   = [size(X0,1) size(X0,2) size(X0,3)];
psi1 = sh.matdefmul(psi1,inv(Mw));
psi  = sh.compose(psi1,sh.affine(df,My\Mf));
S    = push_label(X0,psi,dw);
save(fname,'S','dw','Mw');


function S = push_label(X0,psi,dw)
u  = unique(X0(:));
S  = sparse([],[],[],prod(dw),256,8*prod(dw))';
for i=1:numel(u)
    x0  = X0==u(i);
    r3  = squeeze(find(any(any(x0,2),1)));
    p3  = any(x0(:,:,r3),3);
    r2  = squeeze(find(any(p3,1)));
    r1  = squeeze(find(any(p3,2)));
    x   = spm_diffeo('push',single(x0(r1,r2,r3)),psi(r1,r2,r3,:),dw);
    S(u(i)+1,:) = x(:)';
    fprintf('.');
end
fprintf('\n');


function write_pushed_seg(dw,Mw,psi1,My,X0,Mf,fname)
sh   = spm_mb_shape;
df   = [size(X0,1) size(X0,2) size(X0,3)];
psi1 = sh.matdefmul(psi1,inv(Mw));
psi  = sh.compose(psi1,sh.affine(df,My\Mf));
X1   = push_seg(X0,psi,dw);
Nii  = nifti;
Nii.dat = file_array(fname,size(X1),'float32',0,1,0);
Nii.mat  = Mw;
Nii.mat0 = Mw;
Nii.descrip = 'pushed tissues';
create(Nii);
Nii.dat(:,:,:,:) = X1;


function X1 = push_seg(X0,psi,dw)
s    = sum(X0,4);
Xend = max(1 - s,0);
s    = s+Xend;
X1   = zeros([dw(1:3) size(X0,4)+1],'single');
for i=1:(size(X0,4)+1)
    if i<= size(X0,4)
        x0 = single(X0(:,:,:,i)./s);
    else
        x0 = Xend./s;
    end
    X1(:,:,:,i) = spm_diffeo('push',single(x0),single(psi),dw);
    fprintf('.');
end
fprintf('.\n');
