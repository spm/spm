function [SPM] = spm_eeg_contrasts_conv(SPM, Ic)
% function equivalent to spm_contrasts for EEG/MEG designs which use a full
% rank design matrix.
% FORMAT [SPM] = spm_eeg_contrasts_conv(SPM, Ic)
%
% SPM - SPM struct
% Ic  - index of contrast weights in xCon struct
%_______________________________________________________________________
%
% spm_eeg_contrasts_conv is the equivalent to spm_contrasts for designs
% that have a full rank design matrix. spm_contrasts wouldn't work on full
% rank designs, because there are zero degrees of freedom.
% spm_eeg_contrasts_conv just computes contrast images.
%_______________________________________________________________________
% Stefan Kiebel $Id$


% Get and change to results directory
%-----------------------------------------------------------------------
try
    swd   = SPM.swd;
    cd(swd)
catch
    swd   = pwd;
end

% Get contrast definitions (if available)
%-----------------------------------------------------------------------
try
    xCon  = SPM.xCon;
catch
    xCon  = [];
end

xCon(Ic).eidf = 0;


switch(xCon(Ic).STAT)
    
    case 'T'
        fprintf('\t%-32s: %-10s%20s',sprintf('contrast image %2d', Ic),...
            '(spm_add)','...initialising') %-#
        
        Q = find(abs(xCon(Ic).c) > 0);
        V = SPM.Vbeta(Q);
        
        for j = 1:length(Q)
            V(j).pinfo(1:2,:) = V(j).pinfo(1:2,:)*xCon(Ic).c(Q(j));
        end
        
        %-Prepare handle for contrast image
        %-----------------------------------------------------------
        xCon(Ic).Vcon = struct(...
            'fname', sprintf('con_%04d.img', Ic),...
            'dim', SPM.xVol.DIM',...
	    'dt',  [16, spm_platform('bigend')],...
            'mat', SPM.xVol.M,...
            'pinfo', [1, 0, 0]',...
            'descrip', sprintf('SPM contrast - %d: %s', Ic, xCon(Ic).name));
        
        %-Write image
        %-----------------------------------------------------------
        fprintf('%s%20s',repmat(sprintf('\b'),1,20),'...computing')%-#
        xCon(Ic).Vcon = spm_create_vol(xCon(Ic).Vcon);
        xCon(Ic).Vcon.pinfo(1,1) = spm_add(V, xCon(Ic).Vcon);
        xCon(Ic).Vcon = spm_close_vol(xCon(Ic).Vcon);
        xCon(Ic).Vcon = spm_create_vol(xCon(Ic).Vcon, 'noopen');
        xCon(Ic).Vcon = spm_close_vol(xCon(Ic).Vcon);
        
        fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(Ic).Vcon.fname, 't')))%-#
        
        % multi-dimensional contrast, ESS
    case 'F'  %-Implement ESS as sum of squared weighted beta images
        %-----------------------------------------------------------
        fprintf('\t%-32s: %30s',sprintf('ESS image %2d',i),...
            '...computing') %-#
        
        %-Residual (in parameter space) forming mtx
        %-----------------------------------------------------------
        h       = spm_FcUtil('Hsqr',xCon(Ic),SPM.xX.xKXs);
        
        %-Prepare handle for ESS image
        %-----------------------------------------------------------
        xCon(Ic).Vcon = struct(...
            'fname',  sprintf('ess_%04d.img',Ic),...
            'dim',    SPM.xVol.DIM',...
	    'dt',     [16, spm_platform('bigend')],...
            'mat',    SPM.xVol.M,...
            'pinfo',  [1,0,0]',...
            'descrip',sprintf('SPM ESS -contrast %d: %s',Ic,xCon(Ic).name));
        
        %-Write image
        %-----------------------------------------------------------
        fprintf('%s',repmat(sprintf('\b'),1,30))                   %-#
        xCon(Ic).Vcon = spm_create_vol(xCon(Ic).Vcon);
        xCon(Ic).Vcon = spm_resss(SPM.Vbeta,xCon(Ic).Vcon,h);
        xCon(Ic).Vcon = spm_close_vol(xCon(Ic).Vcon);
        xCon(Ic).Vcon = spm_create_vol(xCon(Ic).Vcon,'noopen');
        xCon(Ic).Vcon = spm_close_vol(xCon(Ic).Vcon);

end
SPM.xCon = xCon;
save SPM SPM
