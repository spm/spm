function [SPM] = spm_contrasts(SPM,Ic)
% Fills in SPM.xCon and writes con_????.img, ess_????.img and SPM?_????.img
% FORMAT [SPM] = spm_contrasts(SPM,Ic)
%
% SPM - SPM data structure
% Ic  - indices of xCon to compute
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes, Karl Friston & Jean-Baptiste Poline
% $Id: spm_contrasts.m 4403 2011-07-21 14:04:33Z guillaume $

% Temporary SPM variable to check for any changes to SPM. We want to avoid
% always having to save SPM.mat unless it has changed, because this is
% slow. A side benefit is one can look at results with just read
% privileges.
%--------------------------------------------------------------------------
tmpSPM = SPM;

%-Get and change to results directory
%--------------------------------------------------------------------------
try
    cd(SPM.swd);
end

%-Get contrast definitions (if available)
%--------------------------------------------------------------------------
try
    xCon = SPM.xCon;
catch
    xCon = [];
end

%-Set all contrasts by default
%--------------------------------------------------------------------------
if nargin < 2
    Ic   = 1:length(xCon);
end
Ic(Ic == 0) = [];

%-Map parameter and hyperarameter files
%--------------------------------------------------------------------------
if ~isempty(xCon) && xCon(Ic(1)).STAT == 'P'
    
    %-Conditional estimators and error variance hyperparameters
    %----------------------------------------------------------------------
    Vbeta = SPM.VCbeta;
else
    
    %-OLS estimators and error variance estimate
    %----------------------------------------------------------------------
    Vbeta = SPM.Vbeta;
    VHp   = SPM.VResMS;
end


%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%==========================================================================
spm('Pointer','Watch')
XYZ   = SPM.xVol.XYZ;
for i = 1:length(Ic)
    
    
    %-Canonicalise contrast structure with required fields
    %----------------------------------------------------------------------
    ic  = Ic(i);
    if isempty(xCon(ic).eidf)
        X1o           = spm_FcUtil('X1o',xCon(ic),SPM.xX.xKXs);
        [trMV,trMVMV] = spm_SpUtil('trMV',X1o,SPM.xX.V);
        xCon(ic).eidf = trMV^2/trMVMV;
    end
    
    
    %-Write contrast/ESS images?
    %======================================================================
    if isempty(xCon(ic).Vcon)
        
        switch(xCon(ic).STAT)
            
            case {'T','P'}
                
                if strcmp(xCon(ic).STAT,'P') && strcmp(SPM.PPM.xCon(ic).PSTAT,'F')
                    
                    % Chi^2 Bayesian inference for compound contrast
                    %------------------------------------------------------
                    disp('Chi^2 Bayesian inference for compound contrast');
                    fprintf('\t%-32s: %30s',sprintf('X2 image %2d',ic),...
                        '...computing');                                %-#
                    
                    if isfield(SPM.PPM,'VB')
                        % First level Bayes
                        xCon = spm_vb_x2(SPM,XYZ,xCon,ic);
                    else
                        % Second level Bayes
                        xCon = spm_bayes2_x2(SPM,XYZ,xCon,ic);
                    end
                else
                    %-Implement contrast as sum of scaled beta images
                    %------------------------------------------------------
                    fprintf('\t%-32s: %-10s%20s',sprintf('contrast image %2d',ic),...
                        '(spm_add)','...initialising');                 %-#
                    
                    Q     = find(abs(xCon(ic).c) > 0);
                    V     = Vbeta(Q);
                    for j = 1:length(Q)
                        V(j).pinfo(1:2,:) = V(j).pinfo(1:2,:)*xCon(ic).c(Q(j));
                    end
                    
                    %-Prepare handle for contrast image
                    %------------------------------------------------------
                    xCon(ic).Vcon = struct(...
                        'fname',  sprintf('con_%04d.img',ic),...
                        'dim',    SPM.xVol.DIM',...
                        'dt',     [16 spm_platform('bigend')],...
                        'mat',    SPM.xVol.M,...
                        'pinfo',  [1,0,0]',...
                        'descrip',sprintf('SPM contrast - %d: %s',ic,xCon(ic).name));
                    
                    %-Write image
                    %------------------------------------------------------
                    fprintf('%s%20s',repmat(sprintf('\b'),1,20),'...computing')%-#
                    xCon(ic).Vcon            = spm_create_vol(xCon(ic).Vcon);
                    xCon(ic).Vcon.pinfo(1,1) = spm_add(V,xCon(ic).Vcon);
                    xCon(ic).Vcon            = spm_create_vol(xCon(ic).Vcon);
                    
                    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
                        '...written %s',spm_str_manip(xCon(ic).Vcon.fname,'t')))%-#
                    
                end
                
            case 'F'  %-Implement ESS as sum of squared weighted beta images
                %----------------------------------------------------------
                fprintf('\t%-32s: %30s',sprintf('ESS image %2d',ic),...
                    '...computing');                                    %-#
                
                %-Residual (in parameter space) forming mtx
                %----------------------------------------------------------
                h       = spm_FcUtil('Hsqr',xCon(ic),SPM.xX.xKXs);
                
                %-Prepare handle for ESS image
                %----------------------------------------------------------
                xCon(ic).Vcon = struct(...
                    'fname',  sprintf('ess_%04d.img',ic),...
                    'dim',    SPM.xVol.DIM',...
                    'dt',     [16, spm_platform('bigend')],...
                    'mat',    SPM.xVol.M,...
                    'pinfo',  [1,0,0]',...
                    'descrip',sprintf('SPM ESS -contrast %d: %s',ic,xCon(ic).name));
                
                %-Write image
                %----------------------------------------------------------
                fprintf('%s',repmat(sprintf('\b'),1,30))                %-#
                xCon(ic).Vcon = spm_create_vol(xCon(ic).Vcon);
                xCon(ic).Vcon = spm_resss(Vbeta,xCon(ic).Vcon,h);
                xCon(ic).Vcon = spm_create_vol(xCon(ic).Vcon);
                
                
            otherwise
                %----------------------------------------------------------
                error(['unknown STAT "',xCon(ic).STAT,'"'])
                
        end % (switch(xCon...)
        
    end % (if isempty(xCon(ic)...)
    
    
    %-Write inference SPM/PPM
    %======================================================================
    if isempty(xCon(ic).Vspm) || xCon(ic).STAT == 'P'
        
        % Always update PPM as size threshold, gamma, may have changed
        %------------------------------------------------------------------
        fprintf('\t%-32s: %30s',sprintf('spm{%c} image %2d',xCon(ic).STAT,ic),...
            '...computing');                                            %-#
        
        switch(xCon(ic).STAT)
            
            case 'T'                                 %-Compute SPM{t} image
                %----------------------------------------------------------
                cB    = spm_get_data(xCon(ic).Vcon,XYZ);
                l     = spm_get_data(VHp,XYZ);       % get hyperparamters
                Vc    = xCon(ic).c'*SPM.xX.Bcov*xCon(ic).c;
                SE    = sqrt(l*Vc);                  % and standard error
                Z     = cB./SE;
                str   = sprintf('[%.1f]',SPM.xX.erdf);
                
                
            case 'P'                                 %-Compute PPM{P} image
                %----------------------------------------------------------
                
                if all(strcmp({SPM.PPM.xCon(ic).PSTAT},'T'))
                    % Simple contrast - Gaussian distributed
                    
                    c     = xCon(ic).c;
                    cB    = spm_get_data(xCon(ic).Vcon,XYZ);
                    if isfield(SPM.PPM,'VB');
                        % If posterior sd image for that contrast does
                        % not already exist, then compute it
                        try
                            SPM.PPM.Vcon_sd(ic);
                        catch
                            SPM = spm_vb_contrasts(SPM,XYZ,xCon,ic);
                        end
                        % Read in posterior sd image for contrast
                        Vsd = spm_get_data(SPM.PPM.Vcon_sd(ic),XYZ);
                        VcB = Vsd.^2;
                    else
                        VcB   = c'*SPM.PPM.Cby*c;
                        for j = 1:length(SPM.PPM.l)
                            
                            % hyperparameter and Taylor approximation
                            %----------------------------------------------
                            l   = spm_get_data(SPM.VHp(j),XYZ);
                            VcB = VcB + (c'*SPM.PPM.dC{j}*c)*(l - SPM.PPM.l(j));
                        end
                    end
                    
                    % posterior probability cB > g
                    %------------------------------------------------------
                    Gamma          = xCon(ic).eidf;
                    Z              = 1 - spm_Ncdf(Gamma,cB,VcB);
                    str            = sprintf('[%.2f]',Gamma);
                    %xCon(ic).name = [xCon(ic).name ' ' str];
                else
                    % Compound contrast - Chi^2 distributed
                    fprintf('\t\t%-75s\n','Chi^2 Bayesian inference for compound contrast');
                    fprintf('\t%-32s: %29s\n',' ',' ');
                    d = spm_get_data(xCon(ic).Vcon,XYZ);
                    Z = spm_Xcdf(d,xCon(ic).eidf);
                    
                    str = sprintf('[%1.2f]',xCon(ic).eidf);
                end
                
                
            case 'F'                                 %-Compute SPM{F} image
                %----------------------------------------------------------
                MVM = spm_get_data(xCon(ic).Vcon,XYZ)/trMV;
                RVR = spm_get_data(VHp,XYZ);
                Z   = MVM./RVR;
                str = sprintf('[%.1f,%.1f]',xCon(ic).eidf,SPM.xX.erdf);
                
            otherwise
                %----------------------------------------------------------
                error(['unknown STAT "',xCon(ic).STAT,'"']);
                
        end % (switch(xCon(ic)...)
        
        
        %-Write SPM - statistic image
        %------------------------------------------------------------------
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...writing');      %-#
        
        xCon(ic).Vspm = struct(...
            'fname',  sprintf('spm%c_%04d.img',xCon(ic).STAT,ic),...
            'dim',    SPM.xVol.DIM',...
            'dt',     [16, spm_platform('bigend')],...
            'mat',    SPM.xVol.M,...
            'pinfo',  [1,0,0]',...
            'descrip',sprintf('SPM{%c_%s} - contrast %d: %s',...
            xCon(ic).STAT,str,ic,xCon(ic).name));
        
        tmp           = zeros(SPM.xVol.DIM');
        Q             = cumprod([1,SPM.xVol.DIM(1:2)'])*XYZ - ...
            sum(cumprod(SPM.xVol.DIM(1:2)'));
        tmp(Q)        = Z;
        xCon(ic).Vspm = spm_write_vol(xCon(ic).Vspm,tmp);
        
        clear tmp Z
        fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(ic).Vspm.fname,'t')));   %-#
        
    end % (if isempty(xCon(ic)...)
    
end % (for i = 1:length(Ic))
spm('Pointer','Arrow')

% place xCon back in SPM
%--------------------------------------------------------------------------
SPM.xCon = xCon;

% Check if SPM has changed. Save only if it has.
%--------------------------------------------------------------------------
if ~isequal(tmpSPM,SPM)
    if spm_check_version('matlab','7') >=0
        save('SPM.mat', 'SPM', '-V6');
    else
        save('SPM.mat', 'SPM');
    end
end
