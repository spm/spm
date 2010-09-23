function [xCon,SPM]= spm_bayes2_x2(SPM,XYZ,xCon,ic)
% Compute and write Chi^2 image for 2nd level Bayes
% FORMAT [xCon,SPM]= spm_bayes2_x2(SPM,XYZ,xCon,ic)
%
% SPM  - SPM data structure
% XYZ  - voxel list
% xCon - contrast info
% ic   - contrast number
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $$
        

%- Compound Contrast
%-----------------------------------------------------------------------
c  = xCon(ic).c;
kc = size(c,2);

%-Get posterior beta's
%-----------------------------------------------------------------------
Nk = size(SPM.xX.X,2);

for k=1:Nk,
    beta(k,:) = spm_get_data(SPM.VCbeta(k),XYZ);
end

% Get noise hyperparameters
Nh=length(SPM.PPM.l);
for jj = 1:Nh,
    hyper(jj).l   = spm_get_data(SPM.VHp(jj),XYZ);
end
 
%-Get posterior SD beta's
%-----------------------------------------------------------------------
Nk = size(SPM.xX.X,2);

%-Loop over voxels
%=======================================================================
Nvoxels = size(XYZ,2);
D       = repmat(NaN,reshape(SPM.xVol.DIM(1:3),1,[]));

spm_progress_bar('Init',100,'Estimating posterior contrast covariance','');

for v=1:Nvoxels,
    
    % Reconstruct approximation to voxel wise covariance matrix
    Sigma_post   = SPM.PPM.Cby;
    for jj = 1:Nh,
        % Taylor approximation to posterior covariance
        Sigma_post = Sigma_post + SPM.PPM.dC{jj}*(hyper(jj).l(v) - SPM.PPM.l(jj));
    end
    
    % Get posterior covariance of contrast
    V  =  c' * Sigma_post * c;
   
    % Get posterior mean of contrast
    m = c'*beta(:,v);
    
    % Get Chi^2 value
    d = m'*inv(V)*m;
    
    D(XYZ(1,v),XYZ(2,v),XYZ(3,v)) = d;
    if rem(v,100)==0
        % update progress bar every 100th voxel
        spm_progress_bar('Set',100*v/Nvoxels);
    end
    
end

xCon(ic).eidf=rank(V);

spm_progress_bar('Clear');   

% Create handle
%-----------------------------------------------------------------------
Vhandle = struct(...
    'fname',  sprintf('x2_%04d.img',ic),...
    'dim',    SPM.xVol.DIM',...
    'dt',     [spm_type('float32') spm_platform('bigend')],... 
    'mat',    SPM.xVol.M,...
    'pinfo',  [1,0,0]',...
    'descrip',sprintf('Chi^2 stat for Bayes multivar con %d: %s',ic,xCon(ic).name));

%-Write image
%-----------------------------------------------------------------------
Vhandle = spm_create_vol(Vhandle);
Vhandle = spm_write_vol(Vhandle,D);

xCon(ic).Vcon = Vhandle;

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),...
    sprintf('...written %s',spm_str_manip(Vhandle.fname,'t')));            %-#
