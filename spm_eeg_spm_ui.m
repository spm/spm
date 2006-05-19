function [SPM] = spm_eeg_spm_ui(SPM)
% user interface for calling general linear model specification for EEG
% data
% FORMAT [SPM] = spm_eeg_spm_ui(SPM)
%
%_______________________________________________________________________
%
% Specification of M/EEG designs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel, Karl Friston
% $Id: spm_eeg_spm_ui.m 539 2006-05-19 17:59:30Z Darren $

%-GUI setup
%-----------------------------------------------------------------------
[Finter, Fgraph, CmdLine] = spm('FnUIsetup', 'EEG stats model setup', 0);
spm_help('!ContextHelp', mfilename)


% get design matrix and/or data
%=======================================================================
if ~nargin
	str = 'specify design or data';
	if spm_input(str, 1, 'b', {'design','data'}, [1 0]);
		% specify a design
		%-------------------------------------------------------
		if sf_abort, spm_clf(Finter), return, end
        
        % choose either normal design specification or shortcut
        Oanalysis = spm_input('Choose design options', 1,'m',{'all options', 'ERP/ERF'}, [0 1]);
        SPM.eeg.Oanalysis = Oanalysis;

		SPM = spm_eeg_design(SPM);

		return
	else
		% get design
		%-------------------------------------------------------
		load(spm_select(1,'^SPM\.mat$','Select SPM.mat'));
	end
else
	% get design matrix
	%---------------------------------------------------------------
	SPM = spm_eeg_design(SPM);
end

% check data are specified
%-----------------------------------------------------------------------
try 
    if ~nargin
        % was called by GUI, i.e. user _wants_ to input data
        SPM.xY = rmfield(SPM.xY, 'P');
    end
    
	SPM.xY.P;
    
catch
    if SPM.eeg.Oanalysis == 1
        % ERP analysis, only 2 factors (condition and time)
        Nobs = SPM.eeg.Nlevels{1}; % Nerps
        tmp = spm_select(Nobs, 'image', 'Select ERPs (1st frame only)');
        clear q
        for i = 1:Nobs
            [p1, p2, p3] = spm_fileparts(deblank(tmp(i, :)));
            q{i} = fullfile(p1, [p2 p3]); % removes ,1
        end
        P = strvcat(q);
    else
        % defaults to normal design specification
        % get filenames
        %---------------------------------------------------------------
        
        % Number of observations of factor 1
        Nobs = size(SPM.eeg.Xind{end-1}, 1);
        
        P     = [];
        oldpwd = pwd;
        
        for i = 1:Nobs
            str = sprintf('Select images for ');
            for j = 1:SPM.eeg.Nfactors-1
                str = [str sprintf('%s(%d)', SPM.eeg.factor{j}, SPM.eeg.Xind{end-1}(i, j))];
                if j < SPM.eeg.Nfactors-1, str = [str ', ']; end
            end
            Nimages = sum(all(kron(ones(size(SPM.eeg.Xind{end}, 1), 1), SPM.eeg.Xind{end-1}(i, :)) == SPM.eeg.Xind{end}(:, 1:end-1), 2));
            
            q = spm_select(Nimages, 'image', str, '', oldpwd);
            P = strvcat(P, q);
        end
    end
	% place in data field
	%---------------------------------------------------------------
	SPM.xY.P = P;

end

% Assemble remaining design parameters
%=======================================================================
spm_help('!ContextHelp',mfilename)

%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spm_clf(Finter);
spm('FigName','Configuring, please wait...',Finter,CmdLine);
spm('Pointer','Watch');


% get file identifiers
%=======================================================================

%-Map files
%-----------------------------------------------------------------------
fprintf('%-40s: ','Mapping files')                          	     %-#
VY    = spm_vol(SPM.xY.P);
fprintf('%30s\n','...done')                                 	     %-#

%-check internal consistency of images
%-----------------------------------------------------------------------
spm_check_orientations(VY);

%-place in xY
%-----------------------------------------------------------------------
SPM.xY.VY = VY;

%-Only implicit mask
%=======================================================================
SPM.xM        = struct('T',	ones(length(VY), 1),...
			'TH',	-inf*ones(length(VY), 1),...
			'I',	0,...
			'VM',	{[]},...
			'xs',	struct('Masking','analysis threshold'));
	
%-Design description - for saving and display
%=======================================================================
% SPM.xsDes = struct(...
% 	'Basis_functions',	SPM.xBF.name,...
% 	'Number_of_ERPs',	sprintf('%d', sum(SPM.eeg.Nsub)*SPM.eeg.Ntypes),...
% 	'Sampling_frequency',	sprintf('%0.2f {s}',SPM.xY.RT)...
% 	);
% 

%-Save SPM.mat
%-----------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                        %-#
if spm_matlab_version_chk('7') >= 0
    save('SPM', 'SPM', '-V6');
else
    save('SPM', 'SPM');
end
fprintf('%30s\n','...SPM.mat saved')                                 %-#

 
%-Display Design report
%=======================================================================
fprintf('%-40s: ','Design reporting')                                %-#
fname     = cat(1,{SPM.xY.VY.fname}');
% spm_DesRep('DesMtx',SPM.xX, SPM.xY.P);
fprintf('%30s\n','...done')                                          %-#


%-End: Cleanup GUI
%=======================================================================
spm_clf(Finter)
spm('FigName','Stats: configured',Finter,CmdLine);
spm('Pointer','Arrow')
fprintf('\n\n')


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function abort = sf_abort
%=======================================================================
if exist(fullfile('.','SPM.mat'))
	str = {	'Current directory contains existing SPM file:',...
		'Continuing will overwrite existing file!'};

	abort = spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
	if abort, fprintf('%-40s: %30s\n\n',...
		'Abort...   (existing SPM files)',spm('time')), end
else
	abort = 0;
end
