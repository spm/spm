function spm_bms_partition(BMS)
% compute model partitioning for BMS
% FORMAT spm_bms_partition(BMS)
%
% Input:
% BMS structure
%
% Output:
% alpha and xppm images (example: xppm_subset1.img)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa

% Find graphics window
% -------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');

% Contrast vector
% -------------------------------------------------------------------------
spm_input('Specify contrast vector. Example: [1 1 2 2 3 3]',1,'d');
contrast = spm_input('Contrast vector',2,'e',[]);

% Inference method to plot
% -------------------------------------------------------------------------
method = spm_input('Inference method',3,'b','FFX|RFX',['FFX';'RFX']);

nb_subsets = max(contrast);
nb_models  = length(contrast);

switch method
    
    % Fixed effects
    % ---------------------------------------------------------------------
    case 'FFX'
        
         % Check if ffx exists
         % ================================================================
         if isfield(BMS.map,'ffx')
                     
           nmodels = size(BMS.map.ffx.ppm,2);
           
           % Check number of subsets and nb of models
           % ==============================================================
           if nb_models ~= nmodels || nb_subsets == 1
              msgbox('Error: contrast vector incorrect!')
              return
           end             
        
           data = cell(1,nb_subsets);
           
           % Get data for each subset
           % ==============================================================
           for i = 1:nmodels,
               num = contrast(i);
               data{num} = [data{num};BMS.map.ffx.ppm{i}];
           end
          
           % Create new images by summing old the ppms
           % ==============================================================
           tmp = data{1};
           vec = find(tmp(1,:) == '\');
           dir = tmp(1,1:vec(end));
           
           data_vol = cell(nb_subsets,1);
           ftmp     = cell(nb_subsets,1);
           
           for j = 1:nb_subsets,
               data_vol{j}  = spm_vol(char(data{j}));
               n_models_sub = size(data{j},1);

               ftmp{j}    = 'i1';
               for jj  = 1:n_models_sub-1
                  ftmp{j} = [ftmp{j},sprintf(' + i%d',jj+1)];
               end
               fname      = sprintf('%ssubset%d_ppm.img',dir,j);
               save_fn{j} = fname;
               Vo = calc_im(j,data_vol,fname,ftmp);

           end
           
           BMS.map.ffx.subsets = save_fn;
           file_name           = BMS.fname;
           BMS.xSPM            = [];
           save(file_name,'BMS')
   
           % Return to results
           % ==============================================================
           spm_input('Done',1,'d');
           return 
            
         else
        
           % No ffx found in BMS.mat
           % ==============================================================
           msgbox('Error: no FFX analysis in current BMS.mat!')
           return
        
        end
     
    % Random effects
    % ---------------------------------------------------------------------
    case 'RFX'
        
                % Check if ffx exists
         % ================================================================
         if isfield(BMS.map,'rfx')
                     
           nmodels = size(BMS.map.rfx.ppm,2);
           
           % Check number of subsets and nb of models
           % ==============================================================
           if nb_models ~= nmodels || nb_subsets == 1
              msgbox('Error: contrast vector incorrect!')
              return
           end             
        
           data = cell(1,nb_subsets);
           
           % Get data for each subset
           % ==============================================================
           for i = 1:nmodels,
               num = contrast(i);
               data{num} = [data{num};BMS.map.rfx.ppm{i}];
           end
          
           % Create new images by summing old the ppms
           % ==============================================================
           tmp = data{1};
           vec = find(tmp(1,:) == '\');
           dir = tmp(1,1:vec(end));
           
           data_vol = cell(nb_subsets,1);
           ftmp     = cell(nb_subsets,1);
           
           for j = 1:nb_subsets,
               data_vol{j}  = spm_vol(char(data{j}));
               n_models_sub = size(data{j},1);

               ftmp{j}    = 'i1';
               for jj  = 1:n_models_sub-1
                  ftmp{j} = [ftmp{j},sprintf(' + i%d',jj+1)];
               end
               fname      = sprintf('%ssubset%d_xppm.img',dir,j);
               save_fn{j} = fname;

               Vo = calc_im(j,data_vol,fname,ftmp);

           end
           
           BMS.map.rfx.subsets = save_fn;
           file_name           = BMS.fname;
           BMS.xSPM            = [];
           save(file_name,'BMS')
   
           % Return to results
           % ==============================================================
           spm_input('Done',1,'d');
           return 
            
         else
        
           % No rfx found in BMS.mat
           % ==============================================================
           msgbox('Error: no RFX analysis in current BMS.mat!')
           return
        
        end
end

end

% Function to sum the data (taken from spm_imcalc)
% -------------------------------------------------------------------------
function out = calc_im(j,data_vol,fname,ftmp)

Vi_tmp    = data_vol{j};
Vi        = Vi_tmp(1);

Vo(j) = struct(...
        'fname',    fname,...
        'dim',      Vi.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      Vi.mat,...
        'descrip',  'spm - algebra');

hold = 1; mask = 0; dmtx = 0;
Vi   = data_vol{j};
n    = numel(Vi);
Y    = zeros(Vo(j).dim(1:3));
f    = ftmp{j};

for p = 1:Vo(j).dim(3),
    B = spm_matrix([0 0 -p 0 0 0 1 1 1]);

    if dmtx, X=zeros(n,prod(Vo(j).dim(1:2))); end
    for i = 1:n
        M = inv(B*inv(Vo(j).mat)*Vi(i).mat);
        d = spm_slice_vol(Vi(i),M,Vo(j).dim(1:2),[hold,NaN]);
        if (mask<0), d(isnan(d))=0; end;
        if (mask>0) && ~spm_type(Vi(i).dt(1),'nanrep'), d(d==0)=NaN; end
        if dmtx, X(i,:) = d(:)'; else eval(['i',num2str(i),'=d;']); end
    end

    eval(['Yp = ' f ';'],['error([''Can''''t evaluate "'',f,''".'']);']);
    if prod(Vo(j).dim(1:2)) ~= numel(Yp),
       error(['"',f,'" produced incompatible image.']); end
    if (mask<0), Yp(isnan(Yp))=0; end
    Y(:,:,p) = reshape(Yp,Vo(j).dim(1:2));

end

temp   = [];
temp   = Vo(j);
temp   = spm_write_vol(temp,Y);
out(j) = temp;
                
end