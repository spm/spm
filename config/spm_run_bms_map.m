function out = spm_run_bms_map (job)
% Run Bayesian Model Selection Maps
% SPM job execution function
% takes a harvested job data structure and calls SPM functions to perform
% Bayesian Inference for Model Selection of Log. Evidence Maps  
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%
%
% Bayesian Inference on Model Space:
%
% The Random-effects 'RFX' method is described in Stephan et al. [1] 
% 'Bayesian Model Selection for Group Studies'.
% Output files (for each model): 
%       BMS.mat 
%       Exceedance Probability Maps (*epm.<ext>),
%       Posterior Probability Maps (*ppm.<ext>),
%       Dirichlet Paramters (alpha) Maps (*alpha.<ext>).
%
% The Fixed-effects 'FFX' method adds together the log-evidences over 
% subjects/sessions for each group, then compares the group log-ev's. 
% This is also known as the Group Bayes Factor (GBF) approach [2]. 
% Output files (for each model):
%       BMS.mat 
%       Posterior Probability Maps (*ppm.<ext>).
%
% BMS contains:
%     BMS.fname
%     BMS.map.ffx(rfx).data
%     BMS.map.ffx(rfx).ppm 
%     BMS.map.ffx(rfx).xppm     - only for RFX (this is the expected posterior
%                                 probability map ie. posterior mean)
%     BMS.map.ffx(rfx).epm      - only for RFX (optional) - this is the 
%                                 exceedance probability map 
%     BMS.map.ffx(rfx).alpha    - only for RFX
%
% [1] Rosa et al., 2009, Bayesian Model Selection Maps for Group Studies,
% NeuroImage.
% [2] Stephan et al., 2009, Bayesian Model Selection for Group Studies,
% NeuroImage.
% [3] Penny et al., 2004, Comparing Dynamic Causal Models, NeuroImage.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_run_bms_map.m 5740 2013-11-13 12:00:04Z guillaume $

% Input
% -------------------------------------------------------------------------
direct  = job.dir{1};
fname   = fullfile(direct,'BMS.mat'); % Output filename (with full path)
mask    = length(job.mask{1});        % Mask image
if mask
   mask_image = spm_vol(job.mask);    % Mask image Vol
end
nsamps    = str2num(job.nsamp);       % Number of samples (nmodels > 3)
do_maps   = job.out_file;           
do_ecp    = do_maps > 0;              % Compute Exceedance Probability
do_alpha  = do_maps > 1;              % Compute Alpha Parameters

% Nb. of subjects and models
% -------------------------------------------------------------------------
nsubjs    = size(job.sess_map,2);
nmodels   = size(job.sess_map{1}(1).mod_map,1);
nsess     = size(job.sess_map{1},2);
nnames    = size(job.mod_name,2);
names     = job.mod_name;

% method
% -------------------------------------------------------------------------
if strcmp(job.method_maps,'FFX');
    method = 'FFX';
else
    method = 'RFX';
end

% Name models
% -------------------------------------------------------------------------
if nnames < nmodels
    for nn=1:nmodels-nnames
        names = [names, sprintf('m%d',nn)];
    end
end

if size(unique(names(1:nmodels)),2) < nmodels,
    id = 'Indentical names for different models!';  % Same name!
    error(id);
end

if nsubjs == 1
   method = 'FFX';                              % If only 1 subject do FFX
end

% Sort out log-evidence images dimensions
% -------------------------------------------------------------------------
Vol_models(1,1) = spm_vol(job.sess_map{1}(1).mod_map(1));

first_vol       = Vol_models(1,1);
M               = first_vol{1}.mat;
DIM             = first_vol{1}.dim(1:3)'; 

xdim            = DIM(1); 
ydim            = DIM(2); 
zdim            = DIM(3);
[xords,yords]   = ndgrid(1:xdim,1:ydim);
xords           = xords(:)';  
yords           = yords(:)';
I               = 1:xdim*ydim;
zords_init      = ones(1,xdim*ydim);

% Setup images
% -------------------------------------------------------------------------
switch method

    case 'FFX',   % Fixed Effects
        
        % Check if BMS.mat exists
        if exist(fullfile(job.dir{1},'BMS.mat'),'file')
           load(fname);
           if  isfield(BMS,'map') && isfield(BMS.map,'ffx')
               str = { 'Warning: existing BMS.mat file has been over-written!'};
               msgbox(str)
           end
        end
        
        % Save BMS data
        out.files{1} = fname;
            
        % Create PPM image files for each model
        model_ppm(1:nmodels) = struct(...
        'fname',    '',...
        'dim',      DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      M,...
        'pinfo',    [1 0 0]',...
        'n', [1 1], ...
        'descrip',  '');

        % Load Vols for all subjects/models 
        for i = 1:nmodels
            model_ppm(i).fname   = fullfile(direct,[sprintf('%s_model_ppm',names{i}) spm_file_ext]);
            model_ppm(i).descrip = sprintf('PPM: %s model',names{i});
            BMS.map.ffx.ppm{i}   = model_ppm(i).fname;
            
            for s = 1:nsubjs,
                for se = 1:nsess,
                    nsessi      = size(job.sess_map{s},2);
                    nmodelsi    = size(job.sess_map{s}(se).mod_map,1);
                    if (nsess == nsessi && nmodels == nmodelsi)
                        Vol_models(s,i,se) = spm_vol(job.sess_map{s}(se).mod_map(i));
                        tmp = Vol_models(s,i,se);
                    else
                        msgbox('The number of sessions/models should be the same for all subjects!')
                        return
                    end
                    % Stop if log-ev images have different dimensions
                    if tmp{1}.dim(1)~=xdim || tmp{1}.dim(2)~=ydim || tmp{1}.dim(3)~=zdim
                       error('Log-evidence images must have the same dimensions!')
                    end
                end
            end      
        end

        % Create files
        model_ppm = spm_create_vol(model_ppm);
        BMS.fname = fname;
        
        % Save data and BMS
        BMS.fname = fname;
        BMS.map.ffx.data = job.sess_map;
        save(out.files{1},'BMS', spm_get_defaults('mat.format'));

    case 'RFX',  % Random Effects
        
        % Check if BMS.mat exists
        if exist(fullfile(job.dir{1},'BMS.mat'),'file')
           load(fname);
           if  isfield(BMS,'map') && isfield(BMS.map,'rfx')
               str = { 'Warning: existing BMS.mat file has been over-written!'};
               msgbox(str)
           end
        end
                
        % BMS structure
        out.files{1}   = fname; 
    
        % Create PPM image files for each model
        model_exp_r(1:nmodels) = struct(...
        'fname',    '',...
        'dim',      DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      M,...
        'pinfo',    [1 0 0]',...
        'n', [1 1], ...
        'descrip',  '');
   
        if do_ecp
            % Create EPM image files for each model
            model_xp(1:nmodels) = struct(...
            'fname',    '',...
            'dim',      DIM',...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      M,...
            'pinfo',    [1 0 0]',...
            'n', [1 1], ...
            'descrip',  '');   
        end
        
        if do_alpha
            % Create alpha image files for each model
            model_alpha(1:nmodels) = struct(...
            'fname',    '',...
            'dim',      DIM',...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      M,...
            'pinfo',    [1 0 0]',...
            'n', [1 1], ...
            'descrip',  '');
        end
        
        % Load Vols for all subjects/models
        for i = 1:nmodels
            model_exp_r(i).fname   = fullfile(direct,[sprintf('%s_model_xppm',names{i}) spm_file_ext]);
            model_exp_r(i).descrip = sprintf('Exp_r: %s model',names{i});
            BMS.map.rfx.ppm{i}     = model_exp_r(i).fname;
            if do_ecp
            model_xp(i).fname      = fullfile(direct,[sprintf('%s_model_epm',names{i}) spm_file_ext]);
            model_xp(i).descrip    = sprintf('XP: %s model',names{i});
            BMS.map.rfx.epm{i}     = model_xp(i).fname;
            end
            if do_alpha
            model_alpha(i).fname   = fullfile(direct,[sprintf('%s_model_alpha',names{i}) spm_file_ext]);
            model_alpha(i).descrip = sprintf('Alpha: %s model',names{i});
            BMS.map.rfx.alpha{i}   = model_alpha(i).fname;
            end
            for s = 1:nsubjs,
                for se = 1:nsess,
                    nsessi      = size(job.sess_map{s},2);
                    nmodelsi    = size(job.sess_map{s}(se).mod_map,1);
                    if (nsess == nsessi && nmodels == nmodelsi)
                        Vol_models(s,i,se) = spm_vol(job.sess_map{s}(se).mod_map(i));
                        tmp = Vol_models(s,i,se);
                    else
                        msgbox('The number of sessions/models should be the same for all subjects!')
                        return
                    end
                    % Stop if log-ev images have different dimensions
                    if tmp{1}.dim(1)~=xdim || tmp{1}.dim(2)~=ydim || tmp{1}.dim(3)~=zdim
                       error('Log-evidence images must have the same dimensions!')
                    end
                end
            end 
        end
        
        % Create files     
        model_exp_r              = spm_create_vol(model_exp_r);
        if do_ecp, model_xp      = spm_create_vol(model_xp); end
        if do_alpha, model_alpha = spm_create_vol(model_alpha); end
        
        % Save data and BMS
        BMS.fname = fname;
        BMS.map.rfx.data = job.sess_map;
        save(out.files{1},'BMS', spm_get_defaults('mat.format')); 
    
end


% Progress bar
% -------------------------------------------------------------------------
spm_progress_bar('Init',zdim,'BMS Maps (Inference)','Slices complete');


% Loop through image slices
% -------------------------------------------------------------------------
for z = 1:zdim,
    
    spm_progress_bar('Set',z);                  % Update progress bar
    j = repmat(NaN,xdim,ydim);                  % Init. image values
    
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'Computing maps...')
    str   = sprintf('Slice %d out of %d',z,zdim); % Display slice nb.
    fprintf('\r%-40s: %30s',str,' ')
 
    zords   = z*zords_init;                     % Slice z
    xyz     = [xords(I); yords(I); zords(I)];   % Slice coordinates
    nVox    = size(xyz,2);                      % Nb. of voxels per slice
    
    if mask
        % Voxels inside mask
        mask_xyz  = mask_image{1}.mat\M*[xyz(:,1:nVox);ones(1,nVox)];
        gamma     = spm_get_data(mask_image{1},mask_xyz);
        b         = find(gamma>0.5);            % Voxels in the mask
    else
        b         = 1:nVox;                     % All voxels
    end
    
    z_models        = NaN(nsubjs,nmodels,nVox);       % Data 
    z_models(1,1,:) = spm_get_data(first_vol{1},xyz); % Data: all subs/mods  
    non_nan         = find(~isnan(z_models(1,1,:)));  % Voxels ~NaN


    % Find voxels ~NaN and sum sessions
    % ---------------------------------------------------------------------
    for s = 1:nsubjs,
        for k = 1:nmodels,
                sum_tmp_data    = [];
            for ns = 1:nsess,
                tmp_data        = Vol_models(s,k,ns);
                sum_tmp_data    = [sum_tmp_data; spm_get_data(tmp_data{1},xyz)];
            end
                z_models(s,k,:) = sum(sum_tmp_data,1);
                non_nani        = find(~isnan(z_models(s,k,:)));
                non_nan         = intersect(non_nan,non_nani);
        end
    end

    % Voxels to be analysed
    non_nan = intersect(non_nan,b);    
    Nvoxels = length(non_nan);

    % Method
    % ---------------------------------------------------------------------
    switch method
            
          % Fixed Effects
          % ---------------------------------------------------------------
          case 'FFX',            
                
            
              
              if Nvoxels > 0                % Slice with ~NaN voxels
                  
                zz     = sum(z_models,1);   % Sum all subjects/sessions
                mz     = mean(zz,2);        % Get mean of all models
                zzmean = zeros(1,nmodels,length(mz));
                for jj = 1:nmodels
                    zzmean(1,jj,:) = mz; 
                end
                
                if nmodels==1
                    % Process log Bayes factor image
                    zz  = exp(zz);              % Exponentiate log-bf values
                    % Calculate posterior probabiliy
                    pz  = zeros(1,1,length(zz));
                    pz  = zz./(1+zz);
                    j(non_nan)   = pz(1,1,non_nan);
                    model_ppm(1) = spm_write_plane(model_ppm(k),j,z);
                else
                    zz  = zz-zzmean;            % Subtract mean
                    zz  = exp(zz);              % Exponentiate log-ev values
                    tzz = sum(zz,2);            % Sum exp(log-ev.)
                    
                    % Calculate posterior probabiliy
                    pz  = zeros(1,nmodels,length(zz));
                    for k = 1:nmodels,
                        pz(1,k,:)    = zz(1,k,:)./tzz;
                        j(non_nan)   = pz(1,k,non_nan);
                        model_ppm(k) = spm_write_plane(model_ppm(k),j,z);
                    end
                end            
                
              else
                % Nvoxels = 0
                for k = 1:nmodels,
                    % Write NaN for slice z
                    model_ppm(k) = spm_write_plane(model_ppm(k),j,z);
                end
              end
              
          % Fixed Effects
          % ---------------------------------------------------------------
          case 'RFX',
                
                if Nvoxels > 0
                    % Initialise results
                    exp_r_total              = zeros(Nvoxels,nmodels);
                    if do_ecp, xp_total      = zeros(Nvoxels,nmodels); end
                    if do_alpha, alpha_total = zeros(Nvoxels,nmodels); end

                    % Do BMS in all voxels of slice z
                    for n = 1:Nvoxels,
                        lme = z_models(:,:,non_nan(n));
                        
                        if nmodels==1
                            % Provide evidence for dummy null model
                            lme=0.5*[lme, -lme];
                            [alpha,exp_r,xp] = spm_BMS(lme,nsamps,0,0,do_ecp);
                            
                            exp_r_total(n,:)              = exp_r(1);  % Cond. Expecta.
                            if do_ecp, xp_total(n,:)      = xp(1); end % Exceeda. Prob.
                            if do_alpha, alpha_total(n,:) = alpha(1); end % Dirichlet par.
                        else
                            
                            % Group BMS
                            [alpha,exp_r,xp] = spm_BMS(lme,nsamps,0,0,do_ecp);
                            
                            exp_r_total(n,:)              = exp_r;  % Cond. Expecta.
                            if do_ecp, xp_total(n,:)      = xp; end % Exceeda. Prob.
                            if do_alpha, alpha_total(n,:) = alpha; end % Dirichlet par.
                        end
                    end

                    % Write images
                    for i = 1:nmodels,
                        j(non_nan)     = exp_r_total(:,i);
                        model_exp_r(i) = spm_write_plane(model_exp_r(i),j,z);
                        if do_ecp
                        j(non_nan)     = xp_total(:,i);
                        model_xp(i)    = spm_write_plane(model_xp(i),j,z);
                        end
                        if do_alpha
                        j(non_nan)     = alpha_total(:,i);
                        model_alpha(i) = spm_write_plane(model_alpha(i),j,z);
                        end
                    end
                else
                    % Write images when Nvoxels = 0
                    for i = 1:nmodels,
                        model_exp_r(i) = spm_write_plane(model_exp_r(i),j,z);
                        if do_ecp, model_xp(i) = spm_write_plane(model_xp(i),j,z); end
                        if do_alpha, model_alpha(i) = spm_write_plane(model_alpha(i),j,z); end
                    end
                end
      
    end

end % Loop over slices

% Clear progress bar
% -------------------------------------------------------------------------
spm_progress_bar('Clear');
disp('Done.');
