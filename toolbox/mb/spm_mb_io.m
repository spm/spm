function varargout = spm_mb_io(action,varargin)
% File I/O Multi-Brain functionalities
%
% FORMAT fn      = spm_mb_io('GetImage',datn)
% FORMAT [out,M] = spm_mb_io('GetData',in)
% FORMAT [d,M]   = spm_mb_io('GetSize',fin)
% FORMAT           spm_mb_io('SaveTemplate',mu,sett)
% FORMAT fout    = spm_mb_io('SetData',fin,f)
% FORMAT dat     = spm_mb_io('SavePsi',dat,sett);
%
%__________________________________________________________________________
% Copyright (C) 2019-2020 Wellcome Centre for Human Neuroimaging

% $Id: spm_mb_io.m 7852 2020-05-19 14:00:48Z spm $

switch action
    case 'GetImage'
        [varargout{1:nargout}] = GetImage(varargin{:});
    case 'GetData'
        [varargout{1:nargout}] = GetData(varargin{:});
    case 'GetSize'
        [varargout{1:nargout}] = GetSize(varargin{:});
    case 'SaveTemplate'
        [varargout{1:nargout}] = SaveTemplate(varargin{:});
    case 'SetData'
        [varargout{1:nargout}] = SetData(varargin{:});
    case 'SavePsi'
        [varargout{1:nargout}] = SavePsi(varargin{:});
    otherwise
        error('Unknown function %s.', action)
end
%==========================================================================

%==========================================================================
function out = GetScale(in)
% Return a scale for adding random numbers


if isnumeric(in)
    if isa(in,'integer')
        out = ones([1,1,1,size(in,4)]);
    else
        out = zeros([1,1,1,size(in,4)]);
    end
    return;
end
if isa(in,'char')
    in = nifti(in);
end
if isa(in,'nifti')
    C   = numel(in);
    d   = size(in(1).dat,[1 2 3 4 5]);
    if d(4)>1 && C>1, error('Don''t know what to do with this image data.'); end
    d(4) = max(d(4),C);
    out  = zeros([1,1,1,d(4)]);
    if C>1
        for c=1:C
            dt1 = in(c).dat.dtype(1);
            if dt1=='I' || dt1=='U'
                out(c) = in(c).dat.scl_slope(1);
            end
        end
    else
        dt1 = in(1).dat.dtype(1);
        if dt1=='I' || dt1=='U'
            out(:) = in(1).dat.scl_slope(1);
        end
    end
else
    error('Unknown datatype.');
end
%==========================================================================

%==========================================================================
function fn = GetImage(gmm)
% This is the place to do various image cleaning steps
fn = spm_mb_io('GetData',gmm.f);
C  = size(fn,4);
fn = Mask(fn,gmm.modality);
jitter = GetScale(gmm.f);
jitter = reshape(jitter,[1 1 1 C]);
if any(jitter~=0)
    % Data is an integer type, so to prevent aliasing in the histogram, small
    % random values are added.
    rng('default'); rng(1);
    fn = fn + bsxfun(@times,rand(size(fn)) - 1/2,jitter);
end
%==========================================================================

%==========================================================================
function fn = Mask(fn,modality)
C = size(fn,4);
for c=1:C
    fn(:,:,:,c) = ApplyMask(fn(:,:,:,c),modality(c));
end
%==========================================================================

%==========================================================================
function f = ApplyMask(f,modality)
if modality==2, f(~isfinite(f) | f == 0 | f < - 1020 | f > 3000) = NaN;
else            f(~isfinite(f) | f == 0)                         = NaN;
end
%==========================================================================

%==========================================================================
function [out,Mn] = GetData(in)
Mn = eye(4);
if isnumeric(in)
    out = single(in);
    return
end
if isa(in,'char')
    in = nifti(in);
end
if isa(in,'nifti')
    C  = numel(in);
    d  = size(in(1).dat,[1 2 3 4 5]);
    Mn = in(1).mat;
    if C>1
        d(4) = C;
        out = zeros(d,'single');
        for m=1:C
            out(:,:,:,m) = single(in(m).dat(:,:,:,:,:));
        end
    else
        out = single(in.dat(:,:,:,:,:));
        if numel(d)>4 && d(4)==1
            out = reshape(out,[d(1:3) d(5)]);
        end
    end
else
    error('Unknown datatype.');
end
%==========================================================================

%==========================================================================
function [d,M] = GetSize(fin)
d = [GetDimensions(fin) 1 1 1];
M = d(4);
d = d(1:3);
%==========================================================================

%==========================================================================
function dat  = SavePsi(dat,sett)
for n=1:numel(dat)
    dat(n) = SavePsiSub(dat(n),sett);
end
%==========================================================================

%==========================================================================
function datn = SavePsiSub(datn,sett)

% Parse function settings
B    = sett.B;
Mmu  = sett.ms.Mmu;
d    = datn.dm;
q    = double(datn.q);
Mn   = datn.Mat;
psi1 = GetData(datn.psi);
psi  = spm_mb_shape('Compose',psi1,spm_mb_shape('Affine',d,Mmu\spm_dexpm(q,B)*Mn));
psi  = reshape(bsxfun(@plus, reshape(psi,[prod(d) 3])*Mmu(1:3,1:3)', Mmu(1:3,4)'),[d 3]);
if isnumeric(datn.psi)
    datn.psi = psi;
elseif isa(datn.psi(1),'nifti')
    to_delete   = datn.psi(1).dat.fname;
    [pth,nam,~] = fileparts(datn.psi(1).dat.fname);
    nam         = datn.onam;
    fpth        = fullfile(pth,['y_' nam '.nii']);
    datn.psi(1).dat.fname = fpth;
    datn.psi(1).dat.dim   = [d 1 3];
    datn.psi(1).mat       = datn.Mat;
    datn.psi(1).descrip   = 'Deformation';
    create(datn.psi(1));
    datn.psi(1).dat(:,:,:,:,:) = reshape(psi,[d 1 3]);
    delete(to_delete);
end
%==========================================================================

%==========================================================================
function SaveTemplate(mu,sett)

if ~isfield(sett.mu,'create'), return; end

% Parse function settings
Mmu      = sett.ms.Mmu;

% Log
fa       = file_array(sett.mu.create.mu,size(mu),'float32',0);
Nmu      = nifti;
Nmu.dat  = fa;
Nmu.mat  = Mmu;
Nmu.mat0 = Mmu;
Nmu.descrip = 'Template';
create(Nmu);
Nmu.dat(:,:,:,:) = mu;

if true
    % Softmax
    mu = spm_mb_shape('TemplateK1',mu);
    mu = exp(mu);
    [pth,nam,ext] = fileparts(sett.mu.create.mu);
    nam      = ['softmax' nam(3:end)];
    f        = fullfile(pth,[nam ext]);
    fa       = file_array(f,size(mu),'float32',0);
    Nmu      = nifti;
    Nmu.dat  = fa;
    Nmu.mat  = Mmu;
    Nmu.mat0 = Mmu;
    Nmu.descrip = 'Template (softmax)';
    create(Nmu);
    Nmu.dat(:,:,:,:) = mu;
end
%==========================================================================

%==========================================================================
function fout = SetData(fin,f)
fout = fin;
if isnumeric(fin)
    fout = f;
    return
end
if isa(fin,'char')
    fin = nifti(fin);
end
if isa(fin,'nifti')
    M    = numel(fin);
    if M>1
        for m=1:M
            fout(m).dat(:,:,:,1,:) = f(:,:,:,m,:);
        end
    else
        fout(1).dat(:,:,:,:,:) = reshape(f,size(fout(1).dat));
    end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
function d = GetDimensions(fin)
if isnumeric(fin)
    d = size(fin);
    d = [d 1 1];
    return
end
if isa(fin,'char')
    fin = nifti(fin);
end
if isa(fin,'nifti')
    M    = numel(fin);
    d    = size(fin(1).dat,[1 2 3 4 5]);
    if M>1
        d(4) = M;
    else
        if numel(d)>4 && d(4)==1
            d = [d(1:3) d(5)];
        end
    end
end
%==========================================================================

%==========================================================================
