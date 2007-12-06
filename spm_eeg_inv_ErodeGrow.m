function varargout = spm_eeg_inv_ErodeGrow(varargin)

%=======================================================================
% FORMAT [Iout/Vout] = spm_eeg_inv_ErodeGrow(Iin/Vin,ne,ng,thr_im)
% 
% This routine erodes then grows an image after thresholding it.
% Inputs :
% Iin       - file name of image to erode/grow
% Vin       - full volume of image (as loaded by spm_read_vols)
% ne        - # of erosion steps (default = 3)
% ng        - # of growing steps (default = 6)
% thr_im    - threshold value to apply (default .8)
%             (compared to the maximal value of image)
% Output :
% Iout      - generated image file name
% Vout      - full volume of eroded-grown image
%     
% Notes:
% - if requested, the output filename is created from Iin as
%   [Iin,'_e',num2str(ne),'g',num2str(ng),'.img']
% - If a file name is passed the output is a filename.
%   If a matrix of values is passed, the output is a matrix of values.
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips & Jeremie Mattout
% $Id: spm_eeg_inv_ErodeGrow.m 1020 2007-12-06 20:20:31Z john $

fl_rvol = 0; % Need to load (1) or not (0) the volume from a file
if nargin<2
    Iin     = spm_select(1,'image','Image to erode-grow');
    nIin    = size(Iin,1);
    Vin     = spm_vol(Iin);
    fl_rvol = 1;
    dim     = Vin.dim;
    cr_file = 1; % Create (1) or not (0) a file at the end of the process.
end

if length(varargin)>=1
    if ischar(varargin{1})
        Iin     = varargin{1};
        nIin    = size(Iin,1);
        Vin     = spm_vol(Iin(1,:));
        fl_rvol = 1;
        dim     = Vin.dim;
        cr_file = 1;
    elseif isa(varargin{1},'uint8')
        val     = varargin{1};
        dim     = size(val);
        cr_file = 0;
    elseif isnumeric(varargin{1})
        val     = uint8(varargin{1});
        dim     = size(val);
        cr_file = 0;
    else
        error('Wrong data format');
    end
end

if fl_rvol
    val = spm_loaduint8(Vin);
end

ne = 1;
ng = 1;
thr_im = .1;
if length(varargin)==4
    ne      =  varargin{2}; 
    ng      =  varargin{3}; 
    thr_im  =  varargin{4};
else
    disp('Using default erode/grow parameters: ne=1, ng=1, thr=.1 !')
end

p1 = zeros(1,dim(2),dim(3));
p2 = zeros(dim(1),1,dim(3));
p3 = zeros(dim(1),dim(2),1);

p1_3 = zeros(1,dim(2),dim(3)-1);
p2_1 = zeros(dim(1)-1,1,dim(3));
p3_2 = zeros(dim(1),dim(2)-1,1);

% Erosion 
val  = val>thr_im*double(max(val(:))); 
Nbin = length(find(val(:)));
fprintf('\tOriginal number of voxels %d .\n',Nbin)
for i=1:ne
    temp = val ...
        + cat(1,p1,val(1:end-1,:,:)) + cat(1,val(2:end,:,:),p1) ...
        + cat(2,p2,val(:,1:end-1,:)) + cat(2,val(:,2:end,:),p2) ...
        + cat(3,p3,val(:,:,1:end-1)) + cat(3,val(:,:,2:end),p3) ...
        + cat(1,p1,cat(2,p2_1,val(1:end-1,1:end-1,:))) ...
        + cat(1,p1,cat(2,val(1:end-1,2:end,:),p2_1)) ...
        + cat(1,cat(2,p2_1,val(2:end,1:end-1,:)),p1) ...
        + cat(1,cat(2,val(2:end,2:end,:),p2_1),p1) ...        
        + cat(2,p2,cat(3,p3_2,val(:,1:end-1,1:end-1))) ...
        + cat(2,p2,cat(3,val(:,1:end-1,2:end),p3_2)) ...
        + cat(2,cat(3,p3_2,val(:,1:end-1,1:end-1)),p2) ...
        + cat(2,cat(3,val(:,1:end-1,2:end),p3_2),p2) ...        
        + cat(3,p3,cat(1,p1_3,val(1:end-1,:,1:end-1))) ...
        + cat(3,p3,cat(1,val(2:end,:,1:end-1),p1_3)) ...
        + cat(3,cat(1,p1_3,val(1:end-1,:,1:end-1)),p3) ...
        + cat(3,cat(1,val(2:end,:,1:end-1),p1_3),p3) ;
    val  = temp > 17;
    Nbin = length(find(val(:)));
    fprintf('\tErosion step %d , %d voxels left.\n',i,Nbin)
end

% Growing
for i=1:ng
    temp = val ...
        + cat(1,p1,val(1:end-1,:,:)) + cat(1,val(2:end,:,:),p1) ...
        + cat(2,p2,val(:,1:end-1,:)) + cat(2,val(:,2:end,:),p2) ...
        + cat(3,p3,val(:,:,1:end-1)) + cat(3,val(:,:,2:end),p3) ...
        + cat(1,p1,cat(2,p2_1,val(1:end-1,1:end-1,:))) ...
        + cat(1,p1,cat(2,val(1:end-1,2:end,:),p2_1)) ...
        + cat(1,cat(2,p2_1,val(2:end,1:end-1,:)),p1) ...
        + cat(1,cat(2,val(2:end,2:end,:),p2_1),p1) ...        
        + cat(2,p2,cat(3,p3_2,val(:,1:end-1,1:end-1))) ...
        + cat(2,p2,cat(3,val(:,1:end-1,2:end),p3_2)) ...
        + cat(2,cat(3,p3_2,val(:,1:end-1,1:end-1)),p2) ...
        + cat(2,cat(3,val(:,1:end-1,2:end),p3_2),p2) ...        
        + cat(3,p3,cat(1,p1_3,val(1:end-1,:,1:end-1))) ...
        + cat(3,p3,cat(1,val(2:end,:,1:end-1),p1_3)) ...
        + cat(3,cat(1,p1_3,val(1:end-1,:,1:end-1)),p3) ...
        + cat(3,cat(1,val(2:end,:,1:end-1),p1_3),p3) ;
    val  = temp>1;
    Nbin = length(find(val(:)));
    fprintf('\tGrowing step %d , %d voxels left.\n',i,Nbin)
end

if cr_file
    Vout        = Vin;
    %Vout.pinfo = [1 0 0]';
    Vout        = rmfield(Vout,'pinfo');
    if nIin==1
        Vout.fname = [spm_str_manip(Vin.fname,'r'),'_e',num2str(ne),'g',num2str(ng),'.img'];
    else
        Vout.fname = Iin(2,:);
    end
    spm_write_vol(Vout,val);
    varargout{1} = Vout.fname ;
else
    varargout{1} = val ;
end
