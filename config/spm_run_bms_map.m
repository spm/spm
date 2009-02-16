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
%       Exceedance Probability Maps (*epm.img),
%       Posterior Probability Maps (*ppm.img),
%       Dirichlet Paramters (alpha) Maps (*alpha.img).
%
% The Fixed-effects 'FFX' method adds together the log-evidences over 
% subjects/sessions for each group, then compares the group log-ev's. 
% This is also known as the Group Bayes Factor (GBF) approach [2]. 
% Output files (for each model):
%       BMS.mat 
%       Posterior Probability Maps (*ppm.img).
%
% BMS contains:
%     BMS.fname
%     BMS.map.ffx(rfx).data
%     BMS.map.ffx(rfx).ppm 
%     BMS.map.ffx(rfx).xppm     - only for RFX
%     BMS.map.ffx(rfx).epm      - only for RFX
%     BMS.map.ffx(rfx).alpha    - only for RFX
%
% [1] Stephan et al., (under review), Bayesian Model Selection for Group 
% Studies, NeuroImage.
% [2] Penny et al., 2004, Comparing Dynamic Causal Models, NeuroImage.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_run_bms_map.m 2751 2009-02-16 15:50:26Z maria $

% Input
% -------------------------------------------------------------------------
method = job.method;                % Inference method
direct = job.dir{1};
fname  = [direct,'BMS.mat'];        % Output filename (including directory)
mask   = length(job.mask{1});       % Mask image
if mask
   mask_image = spm_vol(job.mask);  % Mask image Vol
end

% Nb. of subjects and models
% -------------------------------------------------------------------------
nsubjs  = size(job.sess_map,2);
nmodels = size(job.sess_map{1}(1).mod_map,1);
nsess   = size(job.sess_map{1},2);

if nmodels < 2
    msgbox('Please select more than one file')  % Models must be > 1!
    return
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
zords_init      = repmat(1,1,xdim*ydim);

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
            
        % Create PPM .img files for each model
        model_ppm(1:nmodels) = struct(...
        'fname',    '',...
        'dim',      DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      M,...
        'pinfo',    [1 0 0]',...
        'n', [1 1], ...
        'descrip',  '');

        % Load Vols for all subjects/models 
        for i = 1:nmodels,
            model_ppm(i).fname   = sprintf('%smodel%d_ppm.img',direct,i);
            model_ppm(i).descrip = sprintf('PPM: model %d',i);
            BMS.map.ffx.ppm{i} = model_ppm(i).fname;
            
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
        save(out.files{1},'BMS');

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
        
        % Create alpha .img files for each model
        model_alpha(1:nmodels) = struct(...
        'fname',    '',...
        'dim',      DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      M,...
        'pinfo',    [1 0 0]',...
        'n', [1 1], ...
        'descrip',  '');
    
        % Create PPM .img files for each model
        model_exp_r(1:nmodels) = struct(...
        'fname',    '',...
        'dim',      DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      M,...
        'pinfo',    [1 0 0]',...
        'n', [1 1], ...
        'descrip',  '');
   
        % Create EPM .img files for each model
        model_xp(1:nmodels) = struct(...
        'fname',    '',...
        'dim',      DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      M,...
        'pinfo',    [1 0 0]',...
        'n', [1 1], ...
        'descrip',  '');   

        % Load Vols for all subjects/models
        for i = 1:nmodels,
            model_alpha(i).fname   = sprintf('%smodel%d_alpha.img',direct,i);
            model_alpha(i).descrip = sprintf('Alpha: model %d',i);
            BMS.map.rfx.alpha{i} = model_alpha(i).fname;
            model_exp_r(i).fname   = sprintf('%smodel%d_xppm.img',direct,i);
            model_exp_r(i).descrip = sprintf('Exp_r: model %d',i);
            BMS.map.rfx.ppm{i}   = model_exp_r(i).fname;
            model_xp(i).fname      = sprintf('%smodel%d_epm.img',direct,i);
            model_xp(i).descrip    = sprintf('XP: model %d',i);
            BMS.map.rfx.epm{i}   = model_xp(i).fname;
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
        model_alpha   = spm_create_vol(model_alpha);
        model_exp_r   = spm_create_vol(model_exp_r);
        model_xp      = spm_create_vol(model_xp);
        
        % Save data and BMS
        BMS.fname = fname;
        BMS.map.rfx.data = job.sess_map;
        save(out.files{1},'BMS'); 
    
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
                    alpha_total = zeros(Nvoxels,nmodels);
                    exp_r_total = zeros(Nvoxels,nmodels);
                    xp_total    = zeros(Nvoxels,nmodels);

                    % Do BMS in all voxels of slice z
                    for n = 1:Nvoxels,
                        lme = z_models(:,:,non_nan(n));
                        [alpha,exp_r,xp] = spm_BMS(lme); % Group BMS
                        alpha_total(n,:) = alpha;        % Dirichlet par.
                        exp_r_total(n,:) = exp_r;        % Cond. Expecta.
                        xp_total(n,:)    = xp;           % Exceedance Prob.
                    end
                    
                    % Write images
                    for i = 1:nmodels,
                        j(non_nan)     = alpha_total(:,i);
                        model_alpha(i) = spm_write_plane(model_alpha(i),j,z);
                        j(non_nan)     = exp_r_total(:,i);
                        model_exp_r(i) = spm_write_plane(model_exp_r(i),j,z);
                        j(non_nan)     = xp_total(:,i);
                        model_xp(i)    = spm_write_plane(model_xp(i),j,z);
                    end
                else
                    % Write images when Nvoxels = 0
                    for i = 1:nmodels,
                        model_alpha(i) = spm_write_plane(model_alpha(i),j,z);
                        model_exp_r(i) = spm_write_plane(model_exp_r(i),j,z);
                        model_xp(i)    = spm_write_plane(model_xp(i),j,z);
                    end
                end
   
             
    end

end % Loop over slices

% Clear progress bar
% -------------------------------------------------------------------------
spm_progress_bar('Clear');
disp('Done.');