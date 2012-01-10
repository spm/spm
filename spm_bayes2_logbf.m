function [xCon,SPM]= spm_bayes2_logbf(SPM,XYZ,xCon,ic)
% Compute and write log Bayes factor image
% FORMAT [xCon,SPM]= spm_bayes2_logbf(SPM,XYZ,xCon,ic)
%
% SPM  - SPM data structure
% XYZ  - voxel list
% xCon - contrast info
% ic   - contrast number
%__________________________________________________________________________
% Copyright (C) 2010-2011 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_bayes2_logbf.m 4615 2012-01-10 16:56:25Z will $
        

%-Compound Contrast
%--------------------------------------------------------------------------
c  = xCon(ic).c;
kc = size(c,2);

%-Get posterior beta's
%--------------------------------------------------------------------------
Nk = size(SPM.xX.X,2);

for k=1:Nk
    beta(k,:) = spm_get_data(SPM.VCbeta(k),XYZ);
end

%-Get noise hyperparameters
%--------------------------------------------------------------------------
Nh=length(SPM.PPM.l);
for jj = 1:Nh
    hyper(jj).l = spm_get_data(SPM.VHp(jj),XYZ);
end
 
%-Get posterior SD beta's
%--------------------------------------------------------------------------
Nk = size(SPM.xX.X,2);

%-Loop over voxels
%==========================================================================
Nvoxels = size(XYZ,2);
D       = NaN(reshape(SPM.xVol.DIM(1:3),1,[]));

prior_cov=c'*SPM.PPM.Cb*c;
a=zeros(kc,1);

spm_progress_bar('Init',100,'Estimating Bayes factor','');

for v=1:Nvoxels
    
    % Reconstruct approximation to voxel wise covariance matrix
    Sigma_post   = SPM.PPM.Cby;
    for jj = 1:Nh,
        % Taylor approximation to posterior covariance
        Sigma_post = Sigma_post + SPM.PPM.dC{jj}*(hyper(jj).l(v) - SPM.PPM.l(jj));
    end
     
    % Get posterior for this voxel
    post_mean = c'*beta(:,v);
    post_cov = c' * Sigma_post * c;
    
    en=a-post_mean;
    iC=inv(post_cov);
    logbf = - 0.5*spm_logdet(prior_cov)+0.5*spm_logdet(post_cov)+0.5*en'*iC*en;
    
    D(XYZ(1,v),XYZ(2,v),XYZ(3,v)) = logbf;
    if rem(v,100)==0
        % update progress bar every 100th voxel
        spm_progress_bar('Set',100*v/Nvoxels);
    end
    
end

xCon(ic).eidf=rank(post_cov);

spm_progress_bar('Clear');   

%-Create handle
%--------------------------------------------------------------------------
Vhandle = struct(...
    'fname',  [sprintf('logbf_%04d',ic) spm_file_ext],...
    'dim',    SPM.xVol.DIM',...
    'dt',     [spm_type('float32') spm_platform('bigend')],... 
    'mat',    SPM.xVol.M,...
    'pinfo',  [1,0,0]',...
    'descrip',sprintf('Log Bayes factor con %d: %s',ic,xCon(ic).name));

%-Write image
%--------------------------------------------------------------------------
Vhandle = spm_create_vol(Vhandle);
Vhandle = spm_write_vol(Vhandle,D);

xCon(ic).Vcon = Vhandle;

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),...
    sprintf('...written %s',spm_file(Vhandle.fname,'filename')));       %-#
