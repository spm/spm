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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_contrasts_conv.m 539 2006-05-19 17:59:30Z Darren $

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


for i = 1:length(Ic)
        ic  = Ic(i);

        xCon(ic).eidf = 0;

        switch(xCon(ic).STAT)

            case 'T'
                fprintf('\t%-32s: %-10s%20s',sprintf('contrast image %2d', ic),...
                    '(spm_add)','...initialising') %-#

                Q = find(abs(xCon(ic).c) > 0);
                V = SPM.Vbeta(Q);

                for j = 1:length(Q)
                    V(j).pinfo(1:2,:) = V(j).pinfo(1:2,:)*xCon(ic).c(Q(j));
                end

                %-Prepare handle for contrast image
                %-----------------------------------------------------------
                xCon(ic).Vcon = struct(...
                    'fname', sprintf('con_%04d.img', ic),...
                    'dim', SPM.xVol.DIM',...
		    'dt',  [spm_type('float32') spm_platform('bigend')],...
                    'mat', SPM.xVol.M,...
                    'pinfo', [1, 0, 0]',...
                    'descrip', sprintf('SPM contrast - %d: %s', ic, xCon(ic).name));

                %-Write image
                %-----------------------------------------------------------
                fprintf('%s%20s',repmat(sprintf('\b'),1,20),'...computing')%-#
                xCon(ic).Vcon = spm_create_vol(xCon(ic).Vcon);
                xCon(ic).Vcon.pinfo(1,1) = spm_add(V, xCon(ic).Vcon);
                xCon(ic).Vcon = spm_create_vol(xCon(ic).Vcon, 'noopen');

                fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
                    '...written %s',spm_str_manip(xCon(ic).Vcon.fname, 't')))%-#

                % multi-dimensional contrast, ESS
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
                    'dim',    SPM.xVol.DIM',...
		    'dt',     [spm_type('float32') spm_platform('bigend')],...
                    'mat',    SPM.xVol.M,...
                    'pinfo',  [1,0,0]',...
                    'descrip',sprintf('SPM ESS -contrast %d: %s',ic,xCon(ic).name));

                %-Write image
                %-----------------------------------------------------------
                fprintf('%s',repmat(sprintf('\b'),1,30))                   %-#
                xCon(ic).Vcon = spm_create_vol(xCon(ic).Vcon);
                xCon(ic).Vcon = spm_resss(SPM.Vbeta,xCon(ic).Vcon,h);
                xCon(ic).Vcon = spm_create_vol(xCon(ic).Vcon,'noopen');
        end
end
SPM.xCon = xCon;

if spm_matlab_version_chk('7') >= 0
	save('SPM', 'SPM', '-V6');
else
	save('SPM', 'SPM');
end;
