function varargout = spm_batch(varargin)
% Function for batch estimation of pre-configures SPM stats designs
% FORMAT I = spm_batch('StatsEst',SPMcfg1,SPMcfg2,SPMcfg3,...
% SPMcfg1... - full file pathnames to SPM.cfg files containing design configs
% I          - logical vector indicating designs estimated OK
%_______________________________________________________________________
%
% The function allows batch-mode estimation of stats models: Given file
% pathnames for a set of SPMcfg.mat files, spm_spm_run loops through
% them in turn, passing on the design to spm_spm for estimation.
%
% Model estimates are written in the same directory as the SPMcfg.mat
% file. If the model has already been estimated (as evidenced by a
% SPM.mat file next to the SPMcfg.mat file), then the design is skipped
% so that existing results are not overwritten.
%
%                           ----------------
%
% The SPM stats estimation does not need any graphics windows, and
% therefore can be run on a text-only connection, or in batch mode with
% the text output redirected into a file. (Note that some of the text
% output may look funny, as the dynamic status updating that appears in
% the Matlab command window under normal usage uses special characters
% which will look a bit of a mess in a plain text file.)
%
%                           ----------------
%
% To run as a batch job on UNIX from outside MatLab, you can use the
% redirection operators and/or piping to specify input commands and
% output files. The simplest usage is (for example):
%  echo "spm_batch StatsEst /home/joesmith/analysis/SPMcfg.mat" | matlab \
%	> /home/joesmith/analysis/SPMbatch.log
%
% Alternatively, you can create a small file containing the calls to
% spm_batch . and any other Matlab commands you want to run. Say this
% is called SPMbatch.m, then your batch command would be:
%       matlab < SPMbatch.m > SPMbatch.log
% SPMbatch.log will contain the output that would have appeared in the
% MatLab command window.
%
% This can be combined with the UNIX job scheduling and priority
% facilities (batch, at, nice), to run SPM batch jobs in queues or at
% pre-determined times. Consult your system administrator for further
% details.
%
%-----------------------------------------------------------------------
%
% PLEASE NOTE: SPM is academic software, the authors are not software
% engineers, and we don't have the resources to support SPM as a formal
% software package.
%
% In particular, interface niceties such as these batch facilities are
% not a key part of the academic effort, and are hence provided as is
% without promise of support, continuation, or further development.
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%
SCCSid = '%I%';


if nargin==0, error('No action specified!'), end


%-SPMcfg.mat to run
%=======================================================================
if nargin==0, error('Specify filenames of some SPMcfg.mat files!'), end
SPMcfg = varargin(2:end);

%-SPM setup preliminaries...
%=======================================================================
%-Read SPM defaults if not already set
if isempty(spm('getglobal','CMDLINE')), spm_defaults, end
cwd  = pwd;

%-Function start banner:
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);
fprintf('Batch SPM stats estimation (spm_spm) for designs:\n')
fprintf('\t%s\n',SPMcfg{:}), fprintf('\n')


%-Loop through designs specified, estimating them...
%=======================================================================
str = cell(1,length(SPMcfg));
for i=1:length(SPMcfg)
	swd = spm_str_manip(SPMcfg{i},'h');
	if exist(fullfile(swd,'SPM.mat'),'file')==2
		str{i} = 'Model already estimated!';
	else
		try
			cfg=load(SPMcfg{i});
		catch,
			str{i}='Problem loading SPMcfg';
		end
	end

	if ~isempty(str{i})
		spm('alert!',{['Skipping design: ',SPMcfg{i}],' ',str{i}},...
			'Batch SPM stats estimation',1)
	else
		spm('alert"',{['Design: ',SPMcfg{i}]},...
			'Batch SPM stats estimation',1)
		%-Change to SPMcfg's directory for estimation
		cd(swd)
		
		%-Run spm_spm.m with appropriate parameters
		if isfield(cfg,'Sess') & ~isempty(cfg.Sess)	%-fMRI
			spm_spm(cfg.VY,cfg.xX,cfg.xM,cfg.F_iX0,...
				cfg.Sess,cfg.xsDes)
		elseif isfield(cfg,'xC')			%-PET/SPECT/Basic
			spm_spm(cfg.VY,cfg.xX,cfg.xM,cfg.F_iX0,...
				cfg.xC,cfg.xsDes)
		end
		str{i}='OK';
	end
end


%-Post-estimation cleanup & report
%=======================================================================
cd(cwd)
for i=1:length(SPMcfg), tmp(i) = {sprintf('%s: %s',SPMcfg{i},str{i})}; end
SPMid = spm('FnBanner',mfilename,SCCSid);
spm('alert"',tmp,'Batch SPM stats estimation',1)
varargout = {strcmp(str,'OK')};
