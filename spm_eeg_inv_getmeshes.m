function varargout = spm_eeg_inv_getmeshes(varargin);

%=======================================================================
% Generate the tesselated surfaces of the inner-skull and scalp from binary volumes.
%
% FORMAT D = spm_eeg_inv_getmeshes(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- data struct including the new files and parameters
%
% FORMAT [head,Centre,Nvert] = spm_eeg_inv_getmeshes(Pvol,Center,Nvert);
% Input :
% Pvol      - filenames of volumes 'iskull' and 'scalp'.
% Nvert     - [1x3] provides the number of vertices on each surface
%             Nvert(i) = Npts(i)^2*5/4+2
% Output :
% head      - head structures with 1 to 3 tesselated surfaces
% Centre    - position of the head centre of mass (in mm)
% Nvert     - number of vertices on each surface
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id$

spm_defaults

def_Nvert = [2002 2002];
   
if nargout == 1
    
    try
        D = varargin{1};
    catch
        D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
        D = spm_eeg_ldata(D);
    end
 
    if ~isfield(D,'inv')
        error(sprintf('no inverse structure has been created for this data set\n'));
    end
    
    val = length(D.inv);
    
    if isempty(D.inv{val}.mesh.msk_iskull)
        Iisk = spm_select(1,'.img','Select inner-skull bin mask');
        D.inv{val}.mesh.msk_iskull = Iisk;
    else
        Iisk = D.inv{val}.mesh.msk_iskull;
    end
    
    if isempty(D.inv{val}.mesh.msk_scalp)
        Iscl = spm_select(1,'.img','Select scalp bin mask');
        D.inv{val}.mesh.msk_scalp = Iscl;
    else
        Iscl = D.inv{val}.mesh.msk_scalp;
    end
    
elseif nargout == 3
    
    Pvol = varargin{1};
    
    if length(Pvol) == 2
        Iisk = Pvol{1};
        Iscl = Pvol{2};
    else
        disp('Error: wrong entry Pvol');
        return
    end
    
else
    
    error(sprintf('Error: wrong number of arguments\n'));
    
end

if nargin < 2
    Nvert = def_Nvert;
    [Centre_vx,Centre_mm] = spm_eeg_inv_CtrBin(Iisk);
else
    Centre_mm = varargin{2};
    Nvert = varargin{3};
end  

fprintf('\n\n'), fprintf('%c','='*ones(1,80)), fprintf('\n')
fprintf(['\tGenerate surface meshes from binary volumes\n']);
    
mesh_labels = strvcat('Tess. inner skull','Tess. scalp');
head(1) = spm_eeg_inv_TesBin(Nvert(1),Centre_mm,Iisk,mesh_labels(1,:));    
head(2) = spm_eeg_inv_TesBin(Nvert(2),Centre_mm,Iscl,mesh_labels(2,:));    

fprintf(['\t\tCorrect position of vertices\n']);
head(1) = spm_eeg_inv_ElastM(head(1));
head(2) = spm_eeg_inv_ElastM(head(2));

if nargout == 1
    [pth,nam,ext] = spm_fileparts(D.inv{val}.mesh.msk_iskull);
    meshisk_name  = [nam '_Mesh_' num2str(Nvert(1)) '.mat'];
    D.inv{val}.mesh.tess_iskull = fullfile(pth,meshisk_name);
    vert = head(1).XYZmm';
    face = head(1).tri';
    norm = spm_eeg_inv_normals(vert,face);
    save(D.inv{val}.mesh.tess_iskull,'face','norm','vert');
    D.inv{val}.mesh.Iskull_Nv = Nvert(1);
    D.inv{val}.mesh.Iskull_Nf = length(face);

    [pth,nam,ext] = spm_fileparts(D.inv{val}.mesh.msk_scalp);
    meshscl_name  = [nam '_Mesh_' num2str(Nvert(2)) '.mat'];
    D.inv{val}.mesh.tess_scalp = fullfile(pth,meshscl_name);
    vert = head(2).XYZmm';
    face = head(2).tri';
    norm = spm_eeg_inv_normals(vert,face);
    save(D.inv{val}.mesh.tess_scalp,'face','norm','vert');
    D.inv{val}.mesh.Scalp_Nv = Nvert(2);
    D.inv{val}.mesh.Scalp_Nf = length(face);    

    save(fullfile(D.path, D.fname), 'D');
    varargout{1} = D;
else
    varargout{1} = head;
    varargout{2} = Centre_mm;
    varargout{3} = Nvert;
end
    
fprintf('%c','='*ones(1,80)), fprintf('\n')



%=======================================================================
function norm = spm_eeg_inv_normals(vert,face)

m = struct('Vertices',vert,'Faces',face);

h = figure('Visible','off');
n = get(patch(m),'VertexNormals');
close(h);

f = sqrt(sum(n.^2,2));

I = find(f == 0);
for i = 1:length(I)
    n(I(i)) = n(I(i) - 1);
end

f = sqrt(sum(n.^2,2));

norm(:,1) = n(:,1)./f;
norm(:,2) = n(:,2)./f;
norm(:,3) = n(:,3)./f;

clear m f

return
%=======================================================================