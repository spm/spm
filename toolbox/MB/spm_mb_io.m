function varargout = spm_mb_io(varargin)
% File I/O Multi-Brain functionalities
%
% FORMAT fn      = spm_mb_io('get_image',datn)
% FORMAT [out,M] = spm_mb_io('get_data',in)
% FORMAT [d,M]   = spm_mb_io('get_size',fin)
% FORMAT           spm_mb_io('save_template',mu,sett)
% FORMAT fout    = spm_mb_io('set_data',fin,f)
% FORMAT dat     = spm_mb_io('save_mat',dat,mat);
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


[varargout{1:nargout}] = spm_subfun(localfunctions,varargin{:});
%==========================================================================

%==========================================================================
function delete_files(dat)
for n=1:numel(dat)
    if isa(dat(n).v,'char') && exist(dat(n).v,'file')
        delete(dat(n).v);
    end
    if isa(dat(n).psi,'char') && exist(dat(n).psi,'file')
        delete(dat(n).psi);
    end
end
%==========================================================================

%==========================================================================
function out = save_mat(in,mat)
% Write mat to header
if isa(in,'char')
    if exist(in,'file')
        in = nifti(in);
    else
        out = in;
        return
    end
end
out = in;
if isa(in,'nifti') && numel(in)==1
    out.mat = mat;
    create(out);
end
%==========================================================================

%==========================================================================
function out = get_scale(in)
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
    d   = [in(1).dat.dim 1 1 1 1 1];
    d   = d(1:5);
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
function fn = get_image(gmm)
% This is the place to do various image cleaning steps
fn = get_data(gmm.f);
C  = size(fn,4);
fn = mask(fn,gmm.modality);
jitter = get_scale(gmm.f);
jitter = reshape(jitter,[1 1 1 C]);
if any(jitter~=0)
    % Data is an integer type, so to prevent aliasing in the histogram, small
    % random values are added.
    rng('default'); rng(1);
    fn = fn + bsxfun(@times,rand(size(fn)) - 1/2,jitter);
end
for c=1:numel(gmm.modality)
    if gmm.modality(c)==2
        fn(:,:,:,c) = fn(:,:,:,c) + 1000;
    end
end
%==========================================================================

%==========================================================================
function fn = mask(fn,modality)
C = size(fn,4);
for c=1:C
    fn(:,:,:,c) = apply_mask(fn(:,:,:,c),modality(c));
end
%==========================================================================

%==========================================================================
function f = apply_mask(f,modality)
if modality==2
    f(~isfinite(f) | f == 0 | f < - 1020 | f > 3000) = NaN;
else
    f(~isfinite(f) | f == 0)                         = NaN;
end
%==========================================================================

%==========================================================================
function [out,Mn] = get_data(in)
Mn = eye(4);
if isnumeric(in)
    out = single(in);
    return
end
if isa(in,'char')
    if ~exist(in,'file')
        out = [];
        Mn  =[]; % eye(4);
        return
    end
    in = nifti(in);
end
if isa(in,'nifti')
    C  = numel(in);
    d  = [in(1).dat.dim 1 1 1 1 1];
    d  = d(1:5);
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
function [d,M] = get_size(fin)
d = [get_dimensions(fin) 1 1 1];
M = d(4);
d = d(1:3);
%==========================================================================

%==========================================================================
function save_template(mu,sett)

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
    mu  = spm_mb_shape('softmax0',mu,4);
    d   = [size(mu) 1 1];
    [pth,nam,ext] = fileparts(sett.mu.create.mu);
    nam      = ['softmax' nam(3:end)];
    f        = fullfile(pth,[nam ext]);
    fa       = file_array(f,[d(1:3) d(4)+1],'float32',0);
    Nmu      = nifti;
    Nmu.dat  = fa;
    Nmu.mat  = Mmu;
    Nmu.mat0 = Mmu;
    Nmu.descrip = 'Template (softmax)';
    create(Nmu);
    Nmu.dat(:,:,:,1:d(4)) = mu;
    Nmu.dat(:,:,:,d(4)+1) = max(1-sum(mu,4),0);
end
%==========================================================================

%==========================================================================
function fout = set_data(fin,f)
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
function d = get_dimensions(fin)
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
    d    = [fin(1).dat.dim 1 1 1 1 1];
    d    = d(1:5);
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
