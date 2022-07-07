function spm_sixel(img,col,filename)
% Display or export images in sixel format
% FORMAT spm_sixel(img,col,[filename])
% img       - m x n indexed image or m x n x 3 RGB image
% col       - colormap (three-column matrix of RGB triplets)
% filename  - output filename [default: stdout]
%
% See https://en.wikipedia.org/wiki/Sixel
%__________________________________________________________________________
%
% r = spm_read_vols(spm_vol(fullfile(spm('Dir'),'tpm','TPM.nii,1')));
% g = spm_read_vols(spm_vol(fullfile(spm('Dir'),'tpm','TPM.nii,2')));
% b = spm_read_vols(spm_vol(fullfile(spm('Dir'),'tpm','TPM.nii,3')));
% [img,col] = rgb2ind(cat(3,r(:,:,50),g(:,:,50),b(:,:,50)),64);
% spm_sixel(img,col);
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


%-Input arguments
%--------------------------------------------------------------------------
if ischar(img)
    img = spm_read_vols(spm_vol(img));
    % ANSI/VT100 Terminal Control Escape Sequences
    % http://www.termsys.demon.co.uk/vtansi.htm
    % See also clc, home and fprintf('\033[%dA',N) to move cursor up
    fprintf('\033[s');
    for z=1:size(img,3)
        spm_sixel(repmat(img(:,:,z),1,1,3));
        if z~=size(img,3), fprintf('\033[u'); end
    end
    return;
end
if size(img,3) == 3 % RGB
    if nargin < 2 || isempty(col)
        col = 128;
    end
    [img, col] = rgb2ind(img,col);
end
if nargin < 3
    fid = 1; % stdout
else
    fid = fopen(filename,'wt');
    if fid == -1
        error(sprintf('Cannot open "%s" for writing.',filename));
    end
    onclnp = onCleanup(@() fclose(fid));
end

%-Pad image to be a multiple of 6 (using extra colour beyond colourmap)
%--------------------------------------------------------------------------
extra = mod(size(img,1),6);
if extra, img(end+1:end+(6-extra),:) = size(col,1)+1; end

%-Header
%--------------------------------------------------------------------------
fprintf(fid,'\033Pq\n');

%-Colourmap
%--------------------------------------------------------------------------
col = [(0:size(col,1)-1)', round(col*100)];
n   = floor(log10(size(col,1))) + 1;
fprintf(fid,'#%d;2;%d;%d;%d',col');
fprintf(fid,'\n');

%-Data (RLE not implemented)
%--------------------------------------------------------------------------
for i=1:6:size(img,1)
    l = img(i:i+5,:);  % line of sixels
    c = unique(l(:))'; % unique colours
    if c(end) == size(col,1)+1, c(end) = []; end % padding -> background
    x = repmat(l,[1 1 numel(c)]) == repmat(reshape(c,1,1,[]),[size(l),1]);
    x = 63 + sum(x .* repmat([1 2 4 8 16 32]',[1 size(l,2) numel(c)]),1);
    x = [repmat('#',1,size(x,3)); ...
         reshape(sprintf(['%0' num2str(n) 'd'],c),n,[]); ...
         reshape(char(x),[size(x,2),size(x,3)]); ...
         repmat('$',1,size(x,3))];
    x(end) = '-';
    fprintf(fid,'%s',x);
    fprintf(fid,'\n');
end

%-Footer
%--------------------------------------------------------------------------
fprintf(fid,'\033\\');
