function Vo = spm_file_split(V, odir)
% Convert a 4D volume file into a series of 3D volume files
% FUNCTION Vo = spm_file_split(V, odir)
% V           - filename or spm_vol struct
% odir        - output directory [default: same as input]
%
% Vo          - spm_vol struct array of output files
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_file_split.m 3614 2009-12-07 18:04:24Z guillaume $

if ~nargin
    [V, sts] = spm_select(1,'image','Select a 4D volume file to split');
    if ~sts, return; end
end
if ischar(V)
    [p,n,e] = spm_fileparts(V);
    V = spm_vol(fullfile(p,[n e]));
end

[p,n,e] = spm_fileparts(V(1).fname);
if nargin < 2
    if isempty(p), p = pwd; end
    odir = p;
end

Voo = cell(numel(V),1);
spm_progress_bar('Init',numel(V),'Splitting 4D Volume','Volumes Complete');
for i=1:numel(V)
    Voo{i}        = fullfile(odir,sprintf('%s_%05d%s',n,i,e));
    ni            = nifti;
    ni.dat        = file_array(Voo{i},V(i).dim(1:3),V(i).dt,0,V(i).pinfo(1),V(i).pinfo(2));
    ni.mat        = V(i).mat;
    ni.mat0       = V(i).mat;
    ni.descrip    = [V(i).descrip sprintf(' - %d',i)];
    create(ni);
    ni.dat(:,:,:) = V(i).private.dat(:,:,:,i);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

if nargout
    Vo = spm_vol(char(Voo));
end
