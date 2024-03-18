function out = spm_run_denoise(opt,cfg)
% FORMAT out = spm_run_denoise(opt,cfg)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


switch lower(opt)
    case 'run'
        N = numel(cfg.data{1});
        if any(diff(cellfun(@(x)numel(x),cfg.data)))
            error('Incompatible numbers of images.');
        end
        out.files = cell(N,1);
        for n=1:N
            P    = strvcat(cellfun(@(c)c{n},cfg.data,'UniformOutput',false));
            lam0 = cfg.lambda;
            nit  = cfg.nit;
            dev  = cfg.device;
            out.files{n} = run_denoise(P,lam0,nit,dev);
        end
    case 'vout'
            out(1)            = cfg_dep;
            out(1).sname      = 'Denoised images';
            out(1).src_output = substruct('.','files');
            out(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    otherwise
        error('Incorrect usage.');
end


function out = run_denoise(P,lam0,nit,dev)
% Run total variation denoising
if nargin<4, dev  = 'cpu'; end
if nargin<3, nit  = 200;   end
if nargin<2, lam0 = 30;    end
if nargin<1, P = spm_select(Inf,'nifti'); end

Nii = nifti(P);
vox = sqrt(sum(Nii(1).mat(1:3,1:3).^2));
x   = cellfun(@(f)single(f(:,:,:,:,1,1)),{Nii.dat},'UniformOutput',false);
x   = cat(4,x{:});

if strcmpi(dev,'gpu')
    try
        x   = gpuArray(x);
    catch
        warning('Can''t use GPU: Running on CPU.');
    end
end

% Denoising
sd  = get_sd(x);
y   = spm_TVdenoise(x, vox, lam0./sd, 1./sd.^2, nit);

mx  = gather(max(y(:)));

Nio      = nifti;
[pth,nam,ext] = fileparts(P(1,:));
oname    = fullfile(pth,['denoised_' nam '.nii']);
Nio.dat  = file_array(oname,size(y),'INT16',0,mx./32767,0);
Nio.mat  = Nii(1).mat;
Nio.mat0 = Nii(1).mat0;
Nio.mat_intent  = Nii(1).mat_intent;
Nio.mat0_intent = Nii(1).mat0_intent;
Nio.descrip = sprintf('TV denoised (%g)', lam0);
create(Nio)
Nio.dat(:,:,:,:) = gather(y);
clear Nio
out = oname;



function sd = get_sd(x)
sd = zeros(size(x,4),1);
for n=1:size(x,4)
    mu    = mean(mean(mean(x(:,:,:,n))));
    sd(n) = sqrt(mean(mean(mean((x(:,:,:,n)-mu).^2))));
end

