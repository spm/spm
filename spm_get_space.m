function M = spm_get_space(P,M)
% Get/set the voxel-to-world mapping of an image
% FORMAT M = spm_get_space(P)
%            spm_get_space(P,M)
% P - image filename
% M - voxel-to-world mapping
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 1996-2022 Wellcome Centre for Human Neuroimaging


[pth,nam,ext,num] = spm_fileparts(P);
if ~isempty(num), n = str2num(num(2:end)); else n = [1 1]; end
P = fullfile(pth,[nam ext]);

N = nifti(P);

if nargin==2
    if isempty(N.mat_intent), N.mat = N.mat; end
    N.mat_intent = 'Aligned';
    if n(1)==1

        % Ensure volumes 2..N have the original matrix
        if size(N.dat,4)>1 && sum(sum((N.mat-M).^2))>1e-8
            M0 = N.mat;
            if ~isfield(N.extras,'mat')
                N.extras.mat = zeros([4 4 size(N.dat,4)]);
            else
                if size(N.extras.mat,3)<size(N.dat,4)
                    N.extras.mat(:,:,size(N.dat,4)) = zeros(4);
                end
            end
            for i=2:size(N.dat,4)
                if sum(sum(N.extras.mat(:,:,i).^2))==0
                    N.extras.mat(:,:,i) = M0;
                end
            end
        end

        N.mat = M;
        if strcmp(N.mat0_intent,'Aligned'), N.mat0 = M; end
        if ~isempty(N.extras) && isstruct(N.extras) && isfield(N.extras,'mat') &&...
            size(N.extras.mat,3)>=1
            N.extras.mat(:,:,n(1)) = M;
        end
    else
        N.extras.mat(:,:,n(1)) = M;
    end
    create(N);
else
    if ~isempty(N.extras) && isstruct(N.extras) && isfield(N.extras,'mat') &&...
        size(N.extras.mat,3)>=n(1) && sum(sum(N.extras.mat(:,:,n(1)).^2))
        M = N.extras.mat(:,:,n(1));
    else
        M = N.mat;
    end
end


% Alternative code that would do everything using an additional .mat file.
% Some SPM users don't like realign, coreg etc changing the headers, so
% this code would avoid that - at the expense of an additional file.
function M = alternative_spm_get_space(P,M)
% Get/set the voxel-to-world mapping of an image
% FORMAT M = spm_get_space(P)
%            spm_get_space(P,M)
% P - image filename
% M - voxel-to-world mapping
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 1996-2022 Wellcome Centre for Human Neuroimaging

disp('Warning: Untested spm_get_space.')

[pth,nam,ext,num] = spm_fileparts(P);
if ~isempty(num), n = str2num(num(2:end)); else n = [1 1]; end
P = fullfile(pth,[nam ext]);

N = nifti(P);

if nargin==2
    if isempty(N.mat_intent), N.mat = N.mat; end
    N.mat_intent = 'Aligned';
    M0 = N.mat;
    if ~isfield(N.extras,'mat')
        N.extras.mat = repmat(M0,[1 1 size(N.dat,4)]);
    else
        if size(N.extras.mat,3)<size(N.dat,4)
            r = (size(N.extras.mat,3)+1):size(N.dat,4);
            N.extras.mat(:,:,r) = repmat(M0,[1 1 length(r)]);
        end
    end

    N.mat = M;
    if strcmp(N.mat0_intent,'Aligned'), N.mat0 = M; end
    N.extras.mat(:,:,n(1)) = M;
    N.extras = write_extras(N.dat.fname,N.extras);
else
    if ~isempty(N.extras) && isstruct(N.extras) && isfield(N.extras,'mat') &&...
        size(N.extras.mat,3)>=n(1) && sum(sum(N.extras.mat(:,:,n(1)).^2))
        M = N.extras.mat(:,:,n(1));
    else
        M = N.mat;
    end
end

function extras = write_extras(fname,extras)
% Write extra bits of information
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


if ~isstruct(extras) || isempty(fieldnames(extras))
    return;
end

[pth,nam,ext] = fileparts(fname);
switch ext
    case {'.hdr','.img','.nii'}
        mname = fullfile(pth,[nam '.mat']);
    case {'.HDR','.IMG','.NII'}
        mname = fullfile(pth,[nam '.MAT']);
    otherwise
        mname = fullfile(pth,[nam '.mat']);
end
try
    opt = spm_get_defaults('mat.format');
catch
    opt = '-v6';
end
save(mname,'-struct','extras', opt);

