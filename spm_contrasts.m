function SPM = spm_contrasts(SPM,Ic)
% Compute and store contrast parameters and inference SPM{.}
% FORMAT SPM = spm_contrasts(SPM,Ic)
%
% SPM  - SPM data structure
% Ic   - indices of xCon to compute
%
% This function fills in SPM.xCon and writes con_????, ess_???? and
% spm?_???? images.
%__________________________________________________________________________
% Copyright (C) 2002-2017 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Will Penny & Guillaume Flandin
% $Id: spm_contrasts.m 7738 2019-12-02 12:45:37Z guillaume $


% Temporary copy of the SPM variable, to avoid saving it in SPM.mat unless
% it has changed (faster, read-only access)
%--------------------------------------------------------------------------
tmpSPM = SPM;

%-Change to results directory
%--------------------------------------------------------------------------
try, cd(SPM.swd); end

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
    
    %-Conditional estimators
    %----------------------------------------------------------------------
    Vbeta = SPM.VCbeta;
else
    
    %-OLS estimators and error variance estimate
    %----------------------------------------------------------------------
    Vbeta = SPM.Vbeta;
    VHp   = SPM.VResMS;
end

if spm_mesh_detect(Vbeta)
    file_ext = '.gii';
    g        = SPM.xY.VY(1).private;
    metadata = g.private.metadata;
    name     = {metadata.name};
    if any(ismember(name,'SurfaceID'))
        metadata = metadata(ismember(name,'SurfaceID'));
        metadata = {metadata.name, metadata.value};
    elseif isfield(g,'faces') && ~isempty(g.faces)
        metadata = {'SurfaceID', SPM.xY.VY(1).fname};
    else
        metadata = {};
    end
else
    file_ext = spm_file_ext;
    metadata = {};
end

%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%==========================================================================
spm('Pointer','Watch')
XYZ  = SPM.xVol.XYZ;
iXYZ = cumprod([1,SPM.xVol.DIM(1:2)'])*XYZ - sum(cumprod(SPM.xVol.DIM(1:2)'));
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
        
        switch xCon(ic).STAT
            
            case {'T','P'}
                
                if strcmp(xCon(ic).STAT,'P') && strcmp(SPM.PPM.xCon(ic).PSTAT,'F')
                    
                    % Bayes Factor for compound contrast
                    %------------------------------------------------------
                    disp('Bayes factor for compound contrast');
                    fprintf('\t%-32s: %30s',sprintf('LogBF image %2d',ic),...
                        '...computing');                                %-#
                    
                    if isfield(SPM.PPM,'VB')
                        % First level Bayes
                        xCon = spm_vb_logbf(SPM,XYZ,xCon,ic);
                    else
                        % Second level Bayes
                        xCon = spm_bayes2_logbf(SPM,XYZ,xCon,ic);
                    end
                else
                    %-Implement contrast as linear combination of beta images
                    %------------------------------------------------------
                    fprintf('\t%-32s: %30s',sprintf('contrast image %2d',ic),...
                        '...computing');                                %-#
                    
                    %-Prepare handle for contrast image
                    %------------------------------------------------------
                    xCon(ic).Vcon = struct(...
                        'fname',  [sprintf('con_%04d',ic) file_ext],...
                        'dim',    SPM.xVol.DIM',...
                        'dt',     [spm_type('float32'), spm_platform('bigend')],...
                        'mat',    SPM.xVol.M,...
                        'pinfo',  [1,0,0]',...
                        'descrip',sprintf('Contrast %d: %s',ic,xCon(ic).name),...
                        metadata{:});
                    
                    xCon(ic).Vcon = spm_data_hdr_write(xCon(ic).Vcon);
                    
                    %-Compute contrast
                    %------------------------------------------------------
                    Q      = find(abs(xCon(ic).c) > 0);
                    V      = Vbeta(Q);
                    
                    cB     = zeros(1,size(XYZ,2));
                    for j=1:numel(V)
                        cB = cB + xCon(ic).c(Q(j)) * spm_data_read(V(j),'xyz',XYZ);
                    end
                    
                    %-Write contrast image
                    %------------------------------------------------------
                    tmp = NaN(SPM.xVol.DIM');
                    tmp(iXYZ) = cB;
                    xCon(ic).Vcon = spm_data_write(xCon(ic).Vcon,tmp);
                    
                    clear tmp cB
                    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
                        '...written %s',spm_file(xCon(ic).Vcon.fname,'filename')))%-#
                    
                end
                
            case 'F' %-Implement ESS as sum of squared weighted beta images
                %----------------------------------------------------------
                fprintf('\t%-32s: %30s',sprintf('ESS image %2d',ic),...
                    '...computing');                                    %-#
                                
                %-Prepare handle for ESS image
                %----------------------------------------------------------
                xCon(ic).Vcon = struct(...
                    'fname',  [sprintf('ess_%04d',ic) file_ext],...
                    'dim',    SPM.xVol.DIM',...
                    'dt',     [spm_type('float32'), spm_platform('bigend')],...
                    'mat',    SPM.xVol.M,...
                    'pinfo',  [1,0,0]',...
                    'descrip',sprintf('ESS contrast %d: %s',ic,xCon(ic).name),...
                    metadata{:});
                
                xCon(ic).Vcon = spm_data_hdr_write(xCon(ic).Vcon);
                
                %-Compute ESS
                %----------------------------------------------------------
                % Residual (in parameter space) forming matrix
                h  = spm_FcUtil('Hsqr',xCon(ic),SPM.xX.xKXs);
                ss = zeros(numel(Vbeta),size(XYZ,2));
                for j=1:numel(Vbeta)
                    ss(j,:) = spm_data_read(Vbeta(j),'xyz',XYZ);
                end
                ss = sum((h*ss).^2,1);
                
                %-Write ESS image
                %----------------------------------------------------------
                tmp = NaN(SPM.xVol.DIM');
                tmp(iXYZ) = ss;
                xCon(ic).Vcon = spm_data_write(xCon(ic).Vcon,tmp);
                
                clear tmp ss
                fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
                        '...written %s',spm_file(xCon(ic).Vcon.fname,'filename')))%-#
                
            otherwise
                %----------------------------------------------------------
                error(['unknown STAT "',xCon(ic).STAT,'"'])
                
        end % (switch(xCon...)
        
    end % (if isempty(xCon(ic)...)
    
    
    %-Write inference SPM/PPM
    %======================================================================
    if isempty(xCon(ic).Vspm) || xCon(ic).STAT == 'P'
        % (always update PPM as size threshold, gamma, may have changed)

        fprintf('\t%-32s: %30s',sprintf('spm{%s} image %2d',xCon(ic).STAT,ic),...
            '...computing');                                            %-#
        
        switch(xCon(ic).STAT)
            
            case 'T'                                 %-Compute SPM{t} image
                %----------------------------------------------------------
                cB  = spm_data_read(xCon(ic).Vcon,'xyz',XYZ);
                l   = spm_data_read(VHp,'xyz',XYZ);    % get hyperparamters
                Vc  = xCon(ic).c'*SPM.xX.Bcov*xCon(ic).c;
                SE  = sqrt(l*Vc);                      % and standard error
                Z   = cB./SE;
                str = sprintf('[%.1f]',SPM.xX.erdf);
                
                
            case 'P'                                 %-Compute PPM{P} image
                %----------------------------------------------------------
                
                if all(strcmp({SPM.PPM.xCon(ic).PSTAT},'T'))
                    % Simple contrast - Gaussian distributed
                    
                    c     = xCon(ic).c;
                    cB    = spm_data_read(xCon(ic).Vcon,'xyz',XYZ);
                    if isfield(SPM.PPM,'VB');
                        % If posterior sd image for that contrast does
                        % not already exist, then compute it
                        try
                            SPM.PPM.Vcon_sd(ic);
                        catch
                            SPM = spm_vb_contrasts(SPM,XYZ,xCon,ic);
                        end
                        % Read in posterior sd image for contrast
                        Vsd = spm_data_read(SPM.PPM.Vcon_sd(ic),'xyz',XYZ);
                        VcB = Vsd.^2;
                    else
                        VcB   = c'*SPM.PPM.Cby*c;
                        for j = 1:length(SPM.PPM.l)
                            
                            % hyperparameter and Taylor approximation
                            %----------------------------------------------
                            l   = spm_data_read(SPM.VHp(j),'xyz',XYZ);
                            VcB = VcB + (c'*SPM.PPM.dC{j}*c)*(l - SPM.PPM.l(j));
                        end
                    end
                    
                    % posterior probability cB > g
                    %------------------------------------------------------
                    Gamma          = xCon(ic).eidf;
                    Z              = 1 - spm_Ncdf(Gamma,cB,VcB);
                    
                    % Convert probability to Log Odds Ratio
                    Z              = log( Z ./ (1 - Z+eps) ); 
                    str            = sprintf('[%.2f]',Gamma);
                    %xCon(ic).name = [xCon(ic).name ' ' str];
                else
                    % Compound contrast - Log Bayes Factor
                    fprintf('\t\t%-75s\n','Log Bayes Factor for compound contrast');
                    fprintf('\t%-32s: %29s\n',' ',' ');
                    Z = spm_data_read(xCon(ic).Vcon,'xyz',XYZ);
                    
                    str = sprintf('[%1.2f]',xCon(ic).eidf);
                end
                
                
            case 'F'                                 %-Compute SPM{F} image
                %----------------------------------------------------------
                MVM = spm_data_read(xCon(ic).Vcon,'xyz',XYZ)/trMV;
                RVR = spm_data_read(VHp,'xyz',XYZ);
                Z   = MVM./RVR;
                str = sprintf('[%.1f,%.1f]',xCon(ic).eidf,SPM.xX.erdf);
                
            otherwise
                %----------------------------------------------------------
                error(['unknown STAT "',xCon(ic).STAT,'"']);
                
        end % (switch(xCon(ic)...)
        
        
        %-Write SPM - statistic image
        %------------------------------------------------------------------
        xCon(ic).Vspm = struct(...
            'fname',  [sprintf('spm%s_%04d',xCon(ic).STAT,ic) file_ext],...
            'dim',    SPM.xVol.DIM',...
            'dt',     [spm_type('float32'), spm_platform('bigend')],...
            'mat',    SPM.xVol.M,...
            'pinfo',  [1,0,0]',...
            'descrip',sprintf('SPM{%s_%s} - contrast %d: %s',...
                xCon(ic).STAT,str,ic,xCon(ic).name),...
            metadata{:});
        
        xCon(ic).Vspm = spm_data_hdr_write(xCon(ic).Vspm);
        
        tmp           = zeros(SPM.xVol.DIM');
        tmp(iXYZ)     = Z;
        xCon(ic).Vspm = spm_data_write(xCon(ic).Vspm,tmp);
        
        clear tmp Z
        cmd = sprintf(['[hReg,xSPM,SPM] = spm_results_ui(''Setup'',',...
            'struct(''swd'',''%s'',''Ic'',%d));',...
            'TabDat = spm_list(''List'',xSPM,hReg);'],pwd,ic);
        img = spm_file(spm_file(xCon(ic).Vspm.fname,'filename'),'link',cmd);
        n = 30; if length(img)>n, n = length(img)+n-13; end
        fprintf('%s%*s\n',repmat(sprintf('\b'),1,30),n,sprintf(...
            '...written %s',img)); %-#
        
    end % (if isempty(xCon(ic)...)
    
end % (for i = 1:length(Ic))
spm('Pointer','Arrow')

% place xCon back in SPM
%--------------------------------------------------------------------------
SPM.xCon = xCon;

% Check if SPM has changed. Save only if it has.
%--------------------------------------------------------------------------
if spm_check_version('matlab','8.0') >= 0, my_isequaln = @isequaln;
else my_isequaln = @isequalwithequalnans; end
if ~my_isequaln(tmpSPM,SPM)
    fprintf('\t%-32s: %30s','Saving SPM.mat','...writing');             %-#
    fmt = spm_get_defaults('mat.format');
    s = whos('SPM');
    if s.bytes > 2147483647, fmt = '-v7.3'; end
    save('SPM.mat', 'SPM', fmt);
    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...SPM.mat saved')   %-#
end
