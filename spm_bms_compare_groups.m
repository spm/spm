function con_image = spm_bms_compare_groups(BMSfiles,name,contrast)
% Compare BMS maps for different groups
% FORMAT spm_bms_compare_groups()
%
% Input (interactive):
% BMS             - BMS.mat files for the two groups to compare
% contrast (name) - name of contrast image that will be save in the current
%                   directory 
% contrast (comp) - comparison between groups. options: 'A>B' (posterior 
%                   probability for group 1 > posterior group 2) 
%                   or 'A<B' (posterior probability group 1 < 
%                   posterior for group 2)
%
% Output: contrast image (path)
%__________________________________________________________________________
% Copyright (C) 2009-2011 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_bms_compare_groups.m 4489 2011-09-14 11:27:38Z guillaume $

%-Parameters
%--------------------------------------------------------------------------
if nargin < 3
    BMSfiles = spm_select(2,'^BMS\.mat$','select BMS.mat files');
    name     = spm_input('Contrast name ? ',1,'s');
    contrast = spm_input('Contrast:','+1','b','A>B|A<B',['A>B';'A<B']);    
end

%-Sort out log-evidence images dimensions
%--------------------------------------------------------------------------
load(deblank(BMSfiles(1,:)))
Vol = spm_vol(BMS.map.rfx.alpha{1});

M               = Vol.mat;
DIM             = Vol.dim(1:3)'; 
xdim            = DIM(1); 
ydim            = DIM(2); 
zdim            = DIM(3);
[xords,yords]   = ndgrid(1:xdim,1:ydim);
xords           = xords(:)';  
yords           = yords(:)';
I               = 1:xdim*ydim;
zords_init      = ones(1,xdim*ydim);

%-Setup images
%--------------------------------------------------------------------------
con_image = struct(...
    'fname',    fullfile(pwd,[name spm_file_ext]),...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'n',        [1 1], ...
    'descrip',  sprintf('Contrat image: %s',name));

% Create files
con_image = spm_create_vol(con_image);

% Contrast
switch contrast
    case 'A>B'
        ind = [1 2];
    case 'A<B'
        ind = [2 1];
end
ncon = length(ind);
ngrp = size(BMSfiles,1);

%-Progress bar
%--------------------------------------------------------------------------
spm_progress_bar('Init',zdim,'BMS Maps (Inference)','Slices complete');

%-Loop through image slices
%--------------------------------------------------------------------------
for z = 1:zdim
    
    spm_progress_bar('Set',z);                  % Update progress bar
    j = NaN(xdim,ydim);                         % Init. image values
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'Computing contrast...')
    str   = sprintf('Slice %d out of %d',z,zdim); % Display slice nb.
    fprintf('\r%-40s: %30s',str,' ')
 
    zords   = z*zords_init;                     % Slice z
    xyz     = [xords(I); yords(I); zords(I)];   % Slice coordinates
    nVox    = size(xyz,2);                      % Nb. of voxels per slice
    
    alpha   = NaN(ngrp,ncon,nVox);              % Data 
    for jj=1:ngrp
        load(deblank(BMSfiles(jj,:)))
        for kk=1:ncon
            Vol            = spm_vol(BMS.map.rfx.alpha{ind(kk)});
            alpha(jj,kk,:) = spm_get_data(Vol,xyz); % Data: all subs/mods
        end
    end
    non_nan         = find(~isnan(alpha(1,1,:)));   % Voxels ~NaN
    Nvoxels         = length(non_nan);

    if Nvoxels > 0
        % Initialise results
        con_total = NaN(1,Nvoxels);
        
        % Do BMS in all voxels of slice z
        for n  = 1:Nvoxels,
            alpha1 = alpha(1,:,non_nan(n));
            alpha2 = alpha(2,:,non_nan(n));
            xp = spm_beta_compare(alpha1,alpha2);
            con_total(1,n) = xp;          % Cond. Expecta.
        end
        % Write images
        j(non_nan) = con_total(1,:);
        con_image  = spm_write_plane(con_image,j,z);

    else
        % Write images when Nvoxels = 0
        con_image  = spm_write_plane(con_image,j,z);
    end
        

end % Loop over slices

%-Clear progress bar
%--------------------------------------------------------------------------
spm_progress_bar('Clear');

% Output: path to new file
con_image = con_image.fname;
