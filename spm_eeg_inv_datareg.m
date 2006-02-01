function [varargout] = spm_eeg_inv_datareg(typ,varargin)

%=======================================================================
% Rigid registration of the EEG/MEG data and sMRI spaces
%
% FORMAT D = spm_eeg_inv_datareg(typ,S)
% Input:
% typ       - type of rigid co-registration
%           1: fiducials based (3 landmarks: nazion, left ear, right ear)
%           2: surface matching between sensor mesh and headshap (EEG only)
%           (starts off with a typ 1 registration)
% S		    - input data struct (optional)
% Output:
% D			- data struct including the new files and parameters
%
% FORMAT [RT,sensors_reg,fid_reg,headshape_reg]
%        = spm_eeg_inv_datareg(typ,sensors_file,fid_eeg,fid_mri,headshape,scalpvert)
% Input:
% typ
% sensors_file  - mat file containing a matrix coordinate of the sensor
%               locations ([Sx1 Sy1 Sz1 ; Sx2 Sy2 Sz2 ; ...])
% fid_eeg       - mat file containing the fiducial coordinates in sensor
%               space ([Nx Ny Nz ; LEx LEy LEz ; REx REy REz])
% fid_mri       - mat file containing the fiducial coordinates in sMRI
%               space ([nx ny nz ; lex ley lez ; rex rey rez])
% headshape     - mat file containing the headshape point coordinates in
%               sensor space ([hx1 hy1 hz1 ; hx2 hy2 hz2 ; ...])
%               (only if option typ = 2)
% scalpvert     - mat file containing the vertices coordinates of the scalp
%               tesselation in mri space ([vx1 vy1 vz1 ; vx2 vy2 vz2 ; ...])
%               (only if option typ = 2)
%
% IMPORTANT: all the coordinates must be in the same unit (usually mm).
%
% Output:
% RT            - mat file containing the applied rigid transformation
%                 (Rotation + Translation)
% sensors_reg   - mat file containing the registrated sensor coordinates
%                 in sMRI space
% fid_reg       - mat file containing the registrated file coordinates
%                 in sMRI space
% headshape_reg - mat file conaining the registrated headshap point coordinates
%                 in sMRI space (only if option typ = 2)
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_datareg.m 418 2006-02-01 14:37:08Z jeremie $

spm_defaults

if nargout == 1 
    try
        D = varargin{1};
    catch
        D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
        D = spm_eeg_ldata(D);
    end
    
    if ~isfield(D,'inv')
        disp('Error: no inverse structure has been created for this data set');
        return
    end

    val = length(D.inv);

    if isempty(D.inv{val}.datareg.sens) & D.modality == 'EEG'
        Fpol = spm_input('Read Polhemus files?','+1','Yes|No',[1 2]);
        if Fpol == 1
            D = spm_eeg_inv_ReadPolhemus(D);
        end
    end    
    
    if nargin <= 2
        if isempty(D.inv{val}.datareg.sens)
            sensors_file = spm_select(1,'.mat','Select EEG/MEG sensor file');
            D.inv{val}.datareg.sens = sensors_file;
        else
            sensors_file = D.inv{val}.datareg.sens;
        end
        if isempty(D.inv{val}.datareg.fid)
            fid_eeg = spm_select(1,'.mat','Select fiducials in EEG/MEG space');
            D.inv{val}.datareg.fid = fid_eeg;
        else
            fid_eeg = D.inv{val}.datareg.fid;
        end
        if isempty(D.inv{val}.datareg.fidmri)
            fid_mri = spm_select(1,'.mat','Select fiducials in sMRI space');
            if isempty(fid_mri)
                spm_image('init',D.inv{val}.mesh.sMRI);
                disp(sprintf('\nYou Should Build a .mat file containing the fiducial mm coordinates in MRI space\n'));
                return
            end
            D.inv{val}.datareg.fidmri = fid_mri;
        else
            fid_mri = D.inv{val}.datareg.fidmri;
        end
        if typ == 2
            if isempty(D.inv{val}.datareg.hsp)
                headshape = spm_select(1,'.mat','Select headshape file');
                D.inv{val}.datareg.hsp = headshape;
            else
                headshape = D.inv{val}.datareg.hsp;
            end
            if isempty(D.inv{val}.datareg.scalpvert)
                scalpvert = spm_select(1,'.mat','Select scalp vertices');
                D.inv{val}.datareg.scalpvert = scalpvert;
            else
                scalpvert = D.inv{val}.datareg.scalpvert;
            end
        end
    end           
else
    sensors_file = varargin{1};
    fid_eeg      = varargin{2};
    fid_mri      = varargin{3};    
    if typ == 2
        headshape   = varargin{4};
        scalpvert   = varargin{5};
    end
end

fidmrivar   = load(fid_mri);
name        = fieldnames(fidmrivar);
fidmriloc   = getfield(fidmrivar,name{1});
clear fidmrivar name

fideegvar   = load(fid_eeg);
name        = fieldnames(fideegvar);
fideegloc   = getfield(fideegvar,name{1});
clear fideegvar name

% Landmark based registration (fiducials only)
% The fiducial coordinates must be in the same order in both files
% (usually: NZ & LE, RE)
% If option typ = 1, the registration is limited to this step
% If option typ = 2, the registration is initialised with this step

nfid = size(fideegloc,1);
if nfid ~= size(fidmriloc,1)
    error(sprintf('You need to specify as much MRI as EEG/MEG fiducials'));
end
    
corr = [1 2 3]'*ones(1,2);
[R1, t1] = spm_eeg_inv_rigidreg(fidmriloc', fideegloc', corr);
fideegreg = R1*fideegloc';
fideegreg = (fideegreg + t1*ones(1,nfid))';

% Surface matching between the scalp vertices in MRI space and
% the headshape or sensor positions (EEG only) in data space
% (option typ = 2)

if typ == 2    
    hspvar      = load(headshape);
    name        = fieldnames(hspvar);
    hspvarloc   = getfield(hspvar,name{1});
    if size(hspvarloc,2) > size(hspvarloc,1)
        hspvarloc = hspvarloc';
    end
    clear hspvar name
    
    scalpvar    = load(scalpvert);
    name        = fieldnames(scalpvar);
    scalpmesh   = getfield(scalpvar,'vert');
    if size(scalpmesh,2) > size(scalpmesh,1)
        scalpmesh = scalpmesh';
    end    
    clear scalpvar name  
    
    h = spm_figure;
    set(h,'DoubleBuffer','on','BackingStore','on');
%     cameratoolbar;cameramenu;
    plot3(scalpmesh(:,1),scalpmesh(:,2),scalpmesh(:,3),'ro','MarkerFaceColor','r');
    hold on;
    g = plot3(hspvarloc(:,1),hspvarloc(:,2),hspvarloc(:,3),'bs','MarkerFaceColor','b');
    axis off
    pause(1)

    hspvarloc = R1*hspvarloc';
    hspvarloc = (hspvarloc + t1*ones(1,length(hspvarloc)))';
    
    Tol = max(max(fidmriloc))/1000;
    [Rot,Trans,coor,data2reg] = spm_eeg_inv_icp(scalpmesh',hspvarloc',Tol,g);
    Rot = R1*Rot;
    Trans = R1*Trans + t1;
else
    Rot = R1;
    Trans = t1;
end

sensvar = load(sensors_file);
name = fieldnames(sensvar);
sensloc = getfield(sensvar,name{1});
clear sensvar name
sensreg = Rot*sensloc';
sensreg = (sensreg + Trans*ones(1,length(sensreg)))';

[fpath1,fname,fext] = fileparts(fid_eeg);
fname_fid = [fname '_coreg.mat'];
if nargout == 1
    if str2num(version('-release'))>=14
        save(fullfile(D.path, fname_fid), '-V6', 'fideegreg');
    else
     	save(fullfile(D.path, fname_fid), 'fideegreg');
    end
else
    if str2num(version('-release'))>=14
        save(fullfile(fpath1, fname_fid), '-V6', 'fideegreg');
    else
     	save(fullfile(fpath1, fname_fid), 'fideegreg');
    end
end

[fpath2,fname,fext] = fileparts(sensors_file);
fname_sens = [fname '_coreg.mat'];
if nargout == 1
    if str2num(version('-release'))>=14
        save(fullfile(D.path, fname_sens), '-V6', 'sensreg');
    else
     	save(fullfile(D.path, fname_sens), 'sensreg');
    end
else
    if str2num(version('-release'))>=14
        save(fullfile(fpath2, fname_sens), '-V6', 'sensreg');
    else
     	save(fullfile(fpath2, fname_sens), 'sensreg');
    end
end

fname_transf = ['RegMat.mat'];
if nargout == 1
    if str2num(version('-release'))>=14
        save(fullfile(D.path, fname_transf), '-V6', 'Rot', 'Trans');
    else
     	save(fullfile(D.path, fname_transf), 'Rot','Trans');
    end
else
    if str2num(version('-release'))>=14
        save(fullfile(fpath1, fname_transf), '-V6', 'Rot', 'Trans');
    else
     	save(fullfile(fpath1, fname_transf), 'Rot','Trans');
    end
end

if typ == 2
   [fpath3,fname,fext] = fileparts(headshape);
   fname_hsp = [fname '_coreg.mat'];
   if nargout == 1
       if str2num(version('-release'))>=14
           save(fullfile(D.path, fname_hsp), '-V6', 'data2reg');
       else
       	   save(fullfile(D.path, fname_transf), 'data2reg');
       end
   else
       if str2num(version('-release'))>=14
           save(fullfile(fpath3, fname_hsp), '-V6', 'data2reg');
       else
       	   save(fullfile(fpath3, fname_transf), 'data2reg');
       end
   end   
end

if nargout == 1     
    D.inv{val}.datareg.fid_coreg     = fullfile(fpath1,fname_fid);
    D.inv{val}.datareg.eeg2mri       = fullfile(fpath1,fname_transf);
    D.inv{val}.datareg.sens_coreg    = fullfile(fpath2,fname_sens);   
    if typ == 2
        D.inv{val}.datareg.hsp_coreg     = fullfile(fpath3,fname_hsp);    
    end
    save(fullfile(D.path, D.fname), 'D');
    varargout{1} = D;
elseif nargout == 3
    varargout{1} = fname_transf;
    varargout{2} = fname_sens;
    varargout{3} = fname_fid;
elseif nargout == 4
    varargout{1} = fname_transf;
    varargout{2} = fname_sens;
    varargout{3} = fname_fid;
    varargout{4} = fname_hsp;
end

return

%=======================================================================
function [R, t, corr, data2reg] = spm_eeg_inv_icp(data1, data2, tol, varargin)

% Iterative Closest Point (ICP) registration algorithm.
% Surface matching computation: registration from one 3D surface (set data2 = [Dx1 Dy1 Dz1 ; Dx2 Dy2 Dz2 ; ...])
% onto another 3D surface (set data1 = [dx1 dy1 dz1 ; dx2 dy2 dz2 ; ...])
%
% FORMAT [R,t,corr,data2reg] = spm_eeg_inv_icp(data1,data2,tol)
% Input:
% data1      - locations of the first set of points corresponding to the
%            3D surface to register onto (p points)
% data2      - locations of the second set of points corresponding to the
%            second 3D surface to be registered (m points)
% tol        - convergence criterion: tolerated distance between the
%            coregistrated datasets at iteration k and k+1
% fig_h      - figure handle (optional)
%
% Output:
% R          - 3 x 3 accumulative rotation matrix used to register points in the
%            eeg_file
% t          - 3 x 1 accumulative translation vector used to register
%            eeg_file
% corr       - p x 2 matrix of the index no.s of the corresponding points of
%            scalptess and eeg_file
% data2reg    - 3 x m matrix of the registered points in data2
%=======================================================================
% Adapted from (http://www.csse.uwa.edu.au/~ajmal/icp.m) written by Ajmal Saeed Mian {ajmal@csse.uwa.edu.au}
% Computer Science, The University of Western Australia.

% Jeremie Mattout & Guillaume Flandin

% Landmarks (fiduciales) based registration
% Fiducial coordinates must be given in the same order in both files

if nargin == 4
    fig_h = varargin{1};
end

R = eye(3);
t = zeros(3,1);
tri = delaunayn(data1');
dist_iter = 1;
while dist_iter > tol
    data2reg = data2;
    [corr, D] = dsearchn(data1', tri, data2reg');
    corr(:,2) = [1 : length(corr)]';
    ii = find(D < tol);
    corr(ii,:) = [];
    [R1, t1] = spm_eeg_inv_rigidreg(data1, data2, corr);
    data2 = R1*data2;
    data2 = data2 + t1*ones(1,length(data2));
    R = R1*R;
    t = R1*t + t1;    
	
  	set(fig_h,'XData', data2reg(1,:), 'YData', data2reg(2,:), 'ZData', data2reg(3,:));
    pause(1)
    
    dist_iter = max(sqrt(sum((data2reg - data2).^2)));
end
pause(3)
return
%=======================================================================


%=======================================================================
function [R1, t1] = spm_eeg_inv_rigidreg(data1, data2, corr)
n = length(corr); 
M = data1(:,corr(:,1)); 
mm = mean(M,2); mmstd = std(mm,1,2);
S = data2(:,corr(:,2));
ms = mean(S,2); mstd = std(ms,1,2);
Sshifted = [S(1,:)-ms(1); S(2,:)-ms(2); S(3,:)-ms(3)];
Mshifted = [M(1,:)-mm(1); M(2,:)-mm(2); M(3,:)-mm(3)];
K = Sshifted*Mshifted';
K = K/n;
[U A V] = svd(K);
R1 = V*U';
if det(R1)<0
    B = eye(3);
    B(3,3) = det(V*U');
    R1 = V*B*U';
end
t1 = mm - R1*ms;
return
%=======================================================================
