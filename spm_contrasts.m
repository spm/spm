function [SPM] = spm_contrasts(SPM,Ic)
% Fills in SPM.xCon and writes con_????.img and SPM?_????.img
% FORMAT [SPM] = spm_contrasts(SPM,Ic);
% Ic  - indices of xCon to compute
%_______________________________________________________________________
% @(#)spm_contrasts.m	2.3 Andrew Holmes, Karl Friston & Jean-Baptiste Poline 02/12/30


%-Get and change to results directory
%-----------------------------------------------------------------------
try
	swd   = SPM.swd;
	cd(swd)
catch
	swd   = pwd;
end

%-Get contrast definitions (if available)
%-----------------------------------------------------------------------
try
	xCon  = SPM.xCon;
catch
	xCon  = [];
end

%-set all contrasts by default
%-----------------------------------------------------------------------
if nargin < 2
	Ic    = 1:length(xCon);
end


% map parameter and hyperarameter files
%-----------------------------------------------------------------------
if xCon(Ic(1)).STAT == 'P' 

	% Conditional estimators and error variance hyperparameters
	%----------------------------------------------------------------
    if isfield(SPM.PPM,'VB')
        Vbeta = SPM.VCbeta;
        VHp   = SPM.VHp;
        VPsd  = SPM.VPsd;
        VAR   = SPM.VAR;
    else
        Vbeta = SPM.VCbeta;
        VHp   = SPM.VHp;
    end

else
	% OLS estimators and error variance estimate
	%----------------------------------------------------------------
	Vbeta = SPM.Vbeta;
	VHp   = SPM.VResMS;
end



%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%=======================================================================
spm('Pointer','Watch')
XYZ   = SPM.xVol.XYZ;
for i = 1:length(Ic)


    %-Canonicalise contrast structure with required fields
    %-------------------------------------------------------------------
    ic  = Ic(i);
    if isempty(xCon(ic).eidf)
	X1o           = spm_FcUtil('X1o',xCon(ic),SPM.xX.xKXs);
	[trMV,trMVMV] = spm_SpUtil('trMV',X1o,SPM.xX.V);
        xCon(ic).eidf = trMV^2/trMVMV;
    end


    %-Write contrast/ESS images?
    %-------------------------------------------------------------------
    if isempty(xCon(ic).Vcon)
 
	switch(xCon(ic).STAT)

	    case {'T','P'} %-Implement contrast as sum of scaled beta images
            %-----------------------------------------------------------
            fprintf('\t%-32s: %-10s%20s',sprintf('contrast image %2d',i),...
                                      '(spm_add)','...initialising') %-#

	    Q     = find(abs(xCon(ic).c) > 0);
	    V     = Vbeta(Q);
	    for j = 1:length(Q)
	        V(j).pinfo(1:2,:) = V(j).pinfo(1:2,:)*xCon(ic).c(Q(j));
	    end
	    
	    %-Prepare handle for contrast image
	    %-----------------------------------------------------------
	    xCon(ic).Vcon = struct(...
	        'fname',  sprintf('con_%04d.img',ic),...
                'dim',    [SPM.xVol.DIM',16],...
                'mat',    SPM.xVol.M,...
                'pinfo',  [1,0,0]',...
                'descrip',sprintf('SPM contrast - %d: %s',ic,xCon(ic).name));

        %-Write image
	    %-----------------------------------------------------------
            fprintf('%s%20s',sprintf('\b')*ones(1,20),'...computing')%-#
            xCon(ic).Vcon            = spm_create_vol(xCon(ic).Vcon);
            xCon(ic).Vcon.pinfo(1,1) = spm_add(V,xCon(ic).Vcon);
	    xCon(ic).Vcon            = spm_close_vol(xCon(ic).Vcon);
            xCon(ic).Vcon            = spm_create_vol(xCon(ic).Vcon,'noopen');
	    xCon(ic).Vcon            = spm_close_vol(xCon(ic).Vcon);
            
            fprintf('%s%30s\n',sprintf('\b')*ones(1,30),sprintf(...
                '...written %s',spm_str_manip(xCon(ic).Vcon.fname,'t')))%-#



	    case 'F'  %-Implement ESS as sum of squared weighted beta images
            %-----------------------------------------------------------
            fprintf('\t%-32s: %30s',sprintf('ESS image %2d',i),...
                                                     '...computing') %-#

            %-Residual (in parameter space) forming mtx
	    %-----------------------------------------------------------
            h       = spm_FcUtil('Hsqr',xCon(ic),SPM.xX.xKXs);

	    %-Prepare handle for ESS image
	    %-----------------------------------------------------------
	    xCon(ic).Vcon = struct(...
	        'fname',  sprintf('ess_%04d.img',ic),...
                'dim',    [SPM.xVol.DIM',16],...
                'mat',    SPM.xVol.M,...
                'pinfo',  [1,0,0]',...
                'descrip',sprintf('SPM ESS -contrast %d: %s',ic,xCon(ic).name));

            %-Write image
	    %-----------------------------------------------------------
            fprintf('%s',sprintf('\b')*ones(1,30))                   %-#
            xCon(ic).Vcon = spm_create_vol(xCon(ic).Vcon);
            xCon(ic).Vcon = spm_resss(Vbeta,xCon(ic).Vcon,h);
	    xCon(ic).Vcon = spm_close_vol(xCon(ic).Vcon);
            xCon(ic).Vcon = spm_create_vol(xCon(ic).Vcon,'noopen');
	    xCon(ic).Vcon = spm_close_vol(xCon(ic).Vcon);


	otherwise
        %---------------------------------------------------------------
	    error(['unknown STAT "',xCon(ic).STAT,'"'])

	end % (switch(xCon...)

    end % (if ~isfield...)

    
    
    
    %-Write inference SPM/PPM
    %===================================================================
    if isempty(xCon(ic).Vspm)
	
        fprintf('\t%-32s: %30s',sprintf('spm{%c} image %2d',xCon(ic).STAT,i),...
                                                    '...computing')  %-#

	switch(xCon(ic).STAT)

	case 'T'                                  %-Compute SPM{t} image
        %---------------------------------------------------------------
	cB    = spm_get_data(xCon(ic).Vcon,XYZ);
	l     = spm_get_data(VHp,XYZ);
	VcB   = xCon(ic).c'*SPM.xX.Bcov*xCon(ic).c; 
	Z     = cB./sqrt(l*VcB);
    str   = sprintf('[%.1f]',SPM.xX.erdf);
    

	case 'P'                                  %-Compute PPM{P} image
        %---------------------------------------------------------------
	c     = xCon(ic).c;
	cB    = spm_get_data(xCon(ic).Vcon,XYZ);
    if isfield(SPM.PPM,'VB');
        % Get approximate posterior covariance at ic
        % using Taylor-series approximation
        
        % Get posterior SD beta's
        Nk=size(SPM.xX.X,2);
        for k=1:Nk,
            sd_beta(k,:) = spm_get_data(VPsd(k),XYZ);
        end
        
        % Get AR coefficients
        for p=1:SPM.PPM.AR_P
            a(p,:) = spm_get_data(VAR(p),XYZ);
        end
        
        % Get noise SD
        lambda = spm_get_data(VHp,XYZ);
        
        % Loop over voxels
        Nvoxels=size(XYZ,2);
        for v=1:Nvoxels,
            % Which slice are we in ?
            slice_index=XYZ(3,v);
            
            % Reconstuct approximation to voxel wise correlation matrix
            R=SPM.PPM.slice(slice_index).mean.R;
            dh=a(:,v)'-SPM.PPM.slice(slice_index).mean.a;
            dh=[dh lambda(v)-SPM.PPM.slice(slice_index).mean.lambda];
            for i=1:length(dh),
                R=R+SPM.PPM.slice(slice_index).mean.dR(:,:,i)*dh(i);
            end 
            % Reconstuct approximation to voxel wise covariance matrix
            Sigma_post = (sd_beta(:,v)*sd_beta(:,v)').*R;
            
            VcB (v) = c'*Sigma_post*c;
        end
            
  
    else
        VcB   = c'*SPM.PPM.Cby*c;
        for j = 1:length(SPM.PPM.l)
            
            % hyperparameter and Taylor approximation
            %-------------------------------------------------------
            l   = spm_get_data(VHp(j),XYZ);
            VcB = VcB + (c'*SPM.PPM.dC{j}*c)*(l - SPM.PPM.l(j));
        end
    end
    
	% posterior probability cB > g
        %---------------------------------------------------------------
 	Gamma         = xCon(ic).eidf;
	Z             = 1 - spm_Ncdf(Gamma,cB,VcB);
    str           = sprintf('[%.2f]',Gamma);
	xCon(ic).name = [xCon(ic).name ' ' str];


	case 'F'                                  %-Compute SPM{F} image
        %---------------------------------------------------------------
	MVM   = spm_get_data(xCon(ic).Vcon,XYZ)/trMV;
	RVR   = spm_get_data(VHp,XYZ);
	Z     = MVM./RVR;
        str   = sprintf('[%.1f,%.1f]',xCon(ic).eidf,SPM.xX.erdf);

	otherwise
        %---------------------------------------------------------------
	    error(['unknown STAT "',xCon(ic).STAT,'"'])
	end


        %-Write SPM - statistic image
        %---------------------------------------------------------------
        fprintf('%s%30s',sprintf('\b')*ones(1,30),'...writing')      %-#

        xCon(ic).Vspm = struct(...
	    'fname',  sprintf('spm%c_%04d.img',xCon(ic).STAT,ic),...
	    'dim',    [SPM.xVol.DIM',16],...
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
        fprintf('%s%30s\n',sprintf('\b')*ones(1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(ic).Vspm.fname,'t')))  %-#

    end % (if ~isfield...)

end % (for i = 1:length(Ic))


% place xCon back in SPM
%-----------------------------------------------------------------------
SPM.xCon = xCon;
save SPM SPM
