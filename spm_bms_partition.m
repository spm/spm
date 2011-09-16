function spm_bms_partition(BMS)
% Compute model partitioning for BMS
% FORMAT spm_bms_partition(BMS)
%
% Input:
% BMS structure (BMS.mat)
%
% Output:
% PPM (images) for each of the subsets defined
% xppm_subsetn.img (RFX) and ppm_subsetn.img (FFX)
%__________________________________________________________________________
% Copyright (C) 2009-2011 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_bms_partition.m 4492 2011-09-16 12:11:09Z guillaume $

% Contrast vector
% -------------------------------------------------------------------------
spm_input('Specify contrast vector. Example: [1 1 2 2 3 3]',1,'d');
contrast = spm_input('Contrast vector',2,'e',[]);

% Inference method to plot
% -------------------------------------------------------------------------
method = spm_input('Inference method',3,'b','FFX|RFX',['FFX';'RFX']);

nb_subsets = length(unique(contrast));
max_cont   = max(contrast);
nb_models  = length(contrast);
switch method
    
   case 'FFX'
       
       str_method = 'ffx';
       str_output = 'ppm';     
       
   case 'RFX'
       
       str_method = 'rfx';
       str_output = 'xppm';
       
   otherwise
        
       error('Unknown inference method.');
       
end

% Check if ffx exists
% -------------------------------------------------------------------------
if ~isfield(BMS.map,str_method)
    msgbox(sprintf('No %s analysis in current BMS.mat.',method));
    return
end

% Check number of subsets and nb of models
% -------------------------------------------------------------------------
bms_fields = eval(sprintf('BMS.map.%s.ppm',str_method));
nmodels    = size(bms_fields,2);
        
if nb_models ~= nmodels || nb_subsets == 1 || max_cont ~= nb_subsets
   msgbox('Invalid contrast vector!')
   return
end

% Get data for each subset
% -------------------------------------------------------------------------
data = cell(1,nb_subsets);
        
for i = 1:nmodels,
    num = contrast(i);
    data{num} = [data{num};bms_fields{i}];
end

% Create new images by summing old the ppms
% -------------------------------------------------------------------------
pth      = fileparts(BMS.fname);

data_vol = cell(nb_subsets,1);
ftmp     = cell(nb_subsets,1);

for j = 1:nb_subsets,
    data_vol{j}  = spm_vol(char(data{j}));
    n_models_sub = size(data{j},1);

    ftmp{j}     = 'i1';
    for jj = 1:n_models_sub-1
        ftmp{j} = [ftmp{j},sprintf(' + i%d',jj+1)];
    end
    fname       = fullfile(pth,[sprintf('subset%d_%s',j,str_output) spm_file_ext]);
    save_fn{j}  = fname;
    Vo = calc_im(j,data_vol,fname,ftmp);
end


% Save new BMS structure
% -------------------------------------------------------------------------
bms_struct          = eval(sprintf('BMS.map.%s',str_method));
bms_struct.subsets  = save_fn;
switch method
    case 'FFX'
        BMS.map.ffx = bms_struct;
    case 'RFX'
        BMS.map.rfx = bms_struct;
end
file_name           = BMS.fname;
BMS.xSPM            = [];
save(file_name,'BMS', spm_get_defaults('mat.format'))

% Return to results
%==========================================================================
spm_input('Done',1,'d');

return;

%==========================================================================
% out = calc_im(j,data_vol,fname,ftmp)
%==========================================================================
% Function to sum the data (taken from spm_imcalc)
%--------------------------------------------------------------------------
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
    
    try
        eval(['Yp = ' f ';']);
    catch
        error(['Can''t evaluate "',f,'".']);
    end
    if prod(Vo(j).dim(1:2)) ~= numel(Yp),
       error(['"',f,'" produced incompatible image.']); end
    if (mask<0), Yp(isnan(Yp))=0; end
    Y(:,:,p) = reshape(Yp,Vo(j).dim(1:2));

end

temp   = Vo(j);
temp   = spm_write_vol(temp,Y);
out(j) = temp;
