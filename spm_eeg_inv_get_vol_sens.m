function data = spm_eeg_inv_get_vol_sens(D, val, space, gradsource, modality)
% Retrieves data for leadfield computation from D.inv structure
% Format:
%    D   -  @meeg object
%    val - inversion index (overrides D.val)
%    space - one of 'MNI-aligned', 'Head', 'Native' (default 'MNI-aligned')
%    gradsource - 'inv' (default) to get MEG grad from D.inv
%                 otherwise from D.sensors (useful for reusing head-models
%                 for different runs in the same session)
%    modality - 'EEG' or 'MEG' to force only one modality for multimodal
%                datasets
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_get_vol_sens.m 5803 2013-12-11 16:18:12Z vladimir $

data   = [];

D      = spm_eeg_load(D);

if isempty(D)
    return;
end

if nargin <5 
    modality = [];
end

if nargin < 4 || isempty(gradsource)
    gradsource = 'inv';
end  

if nargin < 3
    space = [];
end

if nargin < 2 || isempty(val)
    val = D.val;
end
   
if ~isfield(D, 'inv')
    error('Please run head model specification.');
end

if numel(D.inv) < val
    error('Invalid inversion index');
end

eegind = 0;
megind = 0;
for m = 1:numel(D.inv{val}.forward)
    if strncmp('EEG', D.inv{val}.forward(m).modality, 3)
        eegind = m;
    elseif strncmp('MEG', D.inv{val}.forward(m).modality, 3)
        megind = m;
    end
end

istemplate = D.inv{val}.mesh.template;

if megind > 0 && ~isequal(modality, 'EEG')
    
    siunits  = isfield(D.inv{val}.forward(megind), 'siunits') &&...
        D.inv{val}.forward(megind).siunits;
    
    datareg  = D.inv{val}.datareg(megind);
    forward  = D.inv{val}.forward(megind);
    
    vol      = forward.vol;
    
    
    if siunits      
        toMNI    =  forward.toMNI;
        to_mm    = diag([1e3 1e3 1e3 1]);
    else
        toMNI    = datareg.toMNI;
        to_mm    = eye(4);
    end
    
    if isequal(gradsource, 'inv')
        if siunits
            sens     =  forward.sensors;           
        else
            sens     = datareg.sensors; 
        end
    else
        sens     = D.sensors('MEG');
    end      
    
    if isfield(forward, 'mesh_correction')
        data.MEG.mesh_correction = forward.mesh_correction;
    else
        data.MEG.mesh_correction = [];
    end
      
    if siunits
        sens  = ft_convert_units(sens, 'm');
    end
    
    M          = to_mm\toMNI;
    [U, L, V]  = svd(M(1:3, 1:3));
    M(1:3,1:3) = U*V';    
    
    if isempty(space)
        space = 'Head';
    end
    
    switch space
        case 'MNI-aligned'            
            data.MEG.vol  = ft_transform_vol(M, vol);
            data.MEG.sens = ft_transform_sens(M, sens);            
            
            data.transforms.toMNI         = toMNI/M;
            data.transforms.toMNI_aligned = to_mm;
            data.transforms.toHead        = inv(M);
            data.transforms.toNative      = D.inv{val}.mesh.Affine\data.transforms.toMNI;
        case 'Head'
            data.MEG.vol  = vol;
            data.MEG.sens = sens;     
            
            data.transforms.toMNI         = toMNI;
            data.transforms.toMNI_aligned = to_mm*M;
            data.transforms.toHead        = eye(4);
            data.transforms.toNative      = D.inv{val}.mesh.Affine\data.transforms.toMNI;
        case 'Native'
           error('Native coordinates option is deprecated for MEG.');
    end
end
            

if eegind > 0 && ~strncmp(modality, 'MEG', 3)
    siunits  = isfield(D.inv{val}.forward(eegind), 'siunits') &&...
        D.inv{val}.forward(eegind).siunits;
        
    forward  = D.inv{val}.forward(eegind);    
    datareg  = D.inv{val}.datareg(eegind);
    
    vol      = forward.vol;
    
    if siunits
        sens     = forward.sensors;
        toMNI    = forward.toMNI;
        to_mm    = diag([1e3 1e3 1e3 1]);
    else
        sens     = datareg.sensors;
        toMNI    = datareg.toMNI;
        to_mm    = eye(4);
    end   
    
    data.EEG.vol  = vol;        
    data.EEG.sens = sens;             
                
    if isfield(forward, 'mesh_correction')
        data.EEG.mesh_correction = forward.mesh_correction;
    else
        data.EEG.mesh_correction = [];
    end
    
    if isfield(data, 'transforms')  % With MEG
        if istemplate
            error('Combining EEG and MEG cannot be done with template head for now.');
        else
            if isa(vol, 'char')
                vol = ft_read_vol(vol);
            end
            
            fromNative = data.transforms.toNative\to_mm;
            
            data.EEG.vol  = ft_transform_vol(fromNative, vol);
            data.EEG.sens = ft_transform_sens(fromNative, sens);
        end
    else                             % EEG only        
        M          = to_mm\toMNI;
        [U, L, V]  = svd(M(1:3, 1:3));
        M(1:3,1:3) = U*V';
        
        if isempty(space)
            space = 'Native';
        end
        
        switch space
            case 'Native'
                data.EEG.vol  = vol;
                data.EEG.sens = sens;
                
                data.transforms.toMNI         = toMNI;
                data.transforms.toMNI_aligned = to_mm*M;
                data.transforms.toHead        = eye(4); 
                data.transforms.toNative      = to_mm; 
            case {'MNI-aligned'}
                if isa(vol, 'char')
                    vol = ft_read_vol(vol);
                end
                
                data.EEG.vol  = ft_transform_vol(M, vol);
                data.EEG.sens = ft_transform_sens(M, sens);
                
                data.transforms.toMNI         = toMNI/M;
                data.transforms.toMNI_aligned = to_mm;
                data.transforms.toHead        = inv(M);
                data.transforms.toNative      = to_mm/M;
            case {'Head'}              
               error('Head space is not defined for EEG data');
        end
    end
end

data.space   = space;
data.siunits = siunits;
