function [SPM] = spm_eeg_spm_ui(SPM)
% user interface for calling general linear model specification for EEG
% data
% FORMAT [SPM] = spm_eeg_spm_ui(SPM)
%
%_______________________________________________________________________
%
% Most of the code is taken from spm_fmri_spm_ui and the functionality is
% the same.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id$

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
		SPM = spm_eeg_design;
% 		spm_DesRep('DesMtx',SPM.xX);
		return
	else
		% get design
		%-------------------------------------------------------
		load(spm_get(1,'SPM.mat','Select SPM.mat'));
	end
else
	% get design matrix
	%---------------------------------------------------------------
	SPM = spm_eeg_design(SPM);
end

% check data are specified
%-----------------------------------------------------------------------
try 
	SPM.xY.P;
catch
	% get filenames
	%---------------------------------------------------------------
    
    % Number of observations of factor 1
    Nobs = size(SPM.eeg.Xind{end-1}, 1);

	P     = [];
	for i = 1:Nobs
		str = sprintf('Select images for ');
        for j = 1:SPM.eeg.Nfactors-1
            str = [str sprintf('%s (%d)', SPM.eeg.factor{j}, SPM.eeg.Xind{end-1}(i, j))];
            if j < 1:SPM.eeg.Nfactors-1, str = [str ', ']; end
        end
        Nimages = sum(all(kron(ones(size(SPM.eeg.Xind{end}, 1), 1), SPM.eeg.Xind{end-1}(i, :)) == SPM.eeg.Xind{end}(:, 1:end-1), 2));
        
        % Problem, spm_get doesn't know about 4D-images
        % q = spm_get(Nimages, '.img', str);
        q = spm_get(inf, '.img', str);
		P = strvcat(P, q);
	end

	% place in data field
	%---------------------------------------------------------------
	SPM.xY.P = P;

end

% Assemble remaining design parameters
%=======================================================================
spm_help('!ContextHelp',mfilename)
% SPM.SPMid = spm('FnBanner',mfilename,SCCSid);


% Global normalization
% definitely skip proportional scaling
% maybe do ERP specific scaling, but still have to look at data
%-----------------------------------------------------------------------
SPM.xGX.sGXcalc = 'mean voxel value';
SPM.xGX.sGMsca  = 'session specific';

% 
%=======================================================================

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% that's the place where the lowpass filter should go.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


% for now, we don't filter the data
% highpass filter has no effect
% HParam = Inf*ones(1,nsess);

% % create and set filter struct
% %---------------------------------------------------------------
% for  i = 1:nsess
% 	K(i) = struct(	'HParam',	HParam(i),...
% 			'row',		SPM.Sess(i).row,...
% 			'RT',		SPM.xY.RT);
% end
% SPM.xX.K = spm_filter(K);

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
if any(any(diff(cat(1,VY.dim),1,1),1) & [1,1,1,0])
error('images do not all have the same dimensions'),           end
if any(any(any(diff(cat(3,VY.mat),1,3),3)))
error('images do not all have same orientation & voxel size'), end

%-place in xY
%-----------------------------------------------------------------------
SPM.xY.VY = VY;


%-Compute Global variate
%=======================================================================
GM    = 100;
q     = length(VY);
% g     = zeros(q, 1);
% fprintf('%-40s: %30s','Calculating globals',' ')                     %-#
% for i = 1:q
% 	fprintf('%s%30s',sprintf('\b')*ones(1,30),sprintf('%4d/%-4d',i,q)) %-#
% 	g(i) = spm_global(VY(i));
% end
% fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#
% 
% gSF = ones(sum(nscan), 1);
% 
% ERP specific scaling
%-----------------------------------------------------------------------
% for i = 1:nsess
%	 gSF(SPM.Sess(i).row) = GM./mean(g(SPM.Sess(i).row));
% end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%-----------------------------------------------------------------------
% for i = 1:q
% 	 SPM.xY.VY(i).pinfo(1:2,:) = SPM.xY.VY(i).pinfo(1:2,:)*gSF(i);
% end

%-place global variates in global structure
%-----------------------------------------------------------------------
% SPM.xGX.rg    = g;
% SPM.xGX.GM    = GM;
% SPM.xGX.gSF   = gSF;
% 

% Don't use mask for 2D interpolated images
% (maybe make it an option when having the choice between 2D and 3D)

%-Masking structure automatically set to 80% of mean
%=======================================================================
if 0
SPM.xM        = struct(	'T',	ones(q,1),...
			'TH',	g.*gSF*0.8,...
			'I',	0,...
			'VM',	{[]},...
			'xs',	struct('Masking','analysis threshold'));
else
SPM.xM        = struct('T',	ones(q,1),...
			'TH',	-inf*ones(q,1),...
			'I',	0,...
			'VM',	{[]},...
			'xs',	struct('Masking','analysis threshold'));
	
end

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
if str2num(version('-release'))>=14
    save('SPM, '-V6', 'SPM');
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
