function spm_bch(varargin)
% SPM batch system: Batch running program.
% FORMAT spm_bch(bch_mfile);
%
% bch_mfile - m-file containing the description of the batch analyses
%
%__________________________________________________________________________
%
% bch_mfile has to contain a variable 'analyses' that describes
% the kind of analyses to perform in batch, and variables:
%             'type', 'work_dir', 'file'
% See spm_bch.man for indication on how to fill these and what they are. 
%
% spm_bch will basically loop over different analyses (eg: 'smooth',
% ..).  Each analyses inputs will be taken from an m-file, (file
% variable) and performed in the working directory describe in
% work_dir.  m-file names can be relative to the working directory.  If
% several analyses of the same kind are to be computed, this is
% resolved with the analyses.index
% 
% The top m-file may also contain the structures for the input of one
% or several analyses
%
%_______________________________________________________________________
% %W% Jean-Baptiste Poline & Stephanie Rouquette %E%

%=======================================================================
% Programmer's notes
%=======================================================================
% To add a new kind of analysis, eg NEW_Analysis
%    1- include it in the switch, 
%       (with BCH.index0  = {'NEW_Analysis',iA(cA)}) 
%    2- modify spm_bch_bchmat.m (see its documentation)
%    3- include 'NEW_Analysis' in the variable 'type' in the 
%       description of the sequences of analyses (top m-file).
%    4- the m-file describing the entries for this new anlysis
%       should contain a top level variable 'NEW_Analysis' 



%- Global batch variables
%-----------------------------------------------------------------------

global BCH

BCH = struct(...
    'bch_mat', '', ...         % mat file or structure in which inputs are read
    'index0'  ,{{'', 0}}, ...  % index for the root level variable
    'flag'   ,1 ...            % flag : 1 -> batch mode, 0 -> GUI mode
);

%- Inputs and Defaults
%-----------------------------------------------------------------------
global MODALITY; 

if nargin == 0, 
   fprintf('FORMAT : spm_bch(bch_mfile [,''FMRI''|''PET'']\n');
   return;
end
if nargin == 1, 
   MODALITY = 'FMRI'; 
   bch_mfile = varargin{1};
end
if nargin == 2, 
   bch_mfile = varargin{1};
   MODALITY  = varargin{2};
   if ~ismember(MODALITY,{'FMRI','PET'}),
      fprintf('FORMAT : spm_bch(bch_mfile [,''FMRI''|''PET'']\n');
      return;
   end
end

%- launch spm defaults for that modality
%-----------------------------------------------------------------------
spm('defaults',MODALITY);

%- m->mat file for the upper level analysis  
%-----------------------------------------------------------------------
try 
   %--------- put it in mat file 
   ana_mat = spm_bch_bchmat(bch_mfile,'analyses');
   %--------- put it in memory. spm_input supports mat file or struct.
   ana_mat = load(ana_mat);
catch
   error(sprintf('\nCan''t evaluate %s of ''analyses'' ', bch_mfile));
   fprintf('\n Try to execute %s in matlab', bch_mfile)
end

BCH.bch_mat = ana_mat;
BCH.index0  = {'analyses',1};

%- save current directory 
%-----------------------------------------------------------------------
cwd = pwd;

%- get indexes for all analyses and the m-files, 
%- working dir, type of analysis 
%-----------------------------------------------------------------------
iA = spm_input('batch',{},'index');
try
   wk_dir = spm_input('batch',{'work_dir',':'});
%- same as: for cA = 1:length(iA)
%-              wk_dir{cA} = spm_input('batch',{'work_dir',cA},1); ...
   typeA  = spm_input('batch',{'type',':'});
   mfileA = spm_input('batch',{'mfile',':'});
catch 
   error(['check the top mfile .... ' bch_mfile]); 
end


for cA = 1:length(iA) % length(iA) == number of analyses to be done

	% go back to current dir (useful if path are relative) 
	%---------------------------------------------------------------
	eval(['cd '  cwd]); 
	%- go into working directory 
	%---------------------------------------------------------------
	try 
	   eval(['cd ' wk_dir{cA}]);
	catch 
	   pwd
	   error(['can''t go to work dir ' ...
		 wk_dir{cA} ' or can''t get in mfile']); 
	end

	%- m->mat file for sub analyses  
   	BCH.bch_mat = spm_bch_bchmat(mfileA{cA},typeA{cA});

	switch typeA{cA}

	    %-----------------------------------------------------------
	    case 'defaults_edit'  
		BCH.index0  = {'defaults_edit',iA(cA)};

		%- first, get the indexes of the default areas
		%-------------------------------------------------------
		index = spm_input('batch',{},'index')
                %- get the list of defaults areas to work on 
		%-------------------------------------------------------
		areas = {};
		areas = spm_input('batch',{'type_area',':'});
                % for i_area = 1:length(index)
                %    areas{i_area} = spm_input('batch',{'type_area',i_area},1);
                % end
		%-------------------------------------------------------
		for i_area = 1:length(areas)
                    %--- little trick : the top variable changes here ... 
                    %--- reading info from {areas{i_area},index(i_area)}
                    BCH.index0 = {areas{i_area},index(i_area)};
                    spm_defaults_edit(areas{i_area});
                end

	    %-----------------------------------------------------------
	    case 'model'
		BCH.index0  = {'model',iA(cA)};
		spm_fmri_spm_ui;
    
	    %-----------------------------------------------------------
	    case 'contrasts'
		BCH.index0  = {'contrasts',iA(cA)};
		s = spm_bch_GetCont; 	
		s = spm_bch_DoCont;		
    
	    %-----------------------------------------------------------
	    case 'headers'
		BCH.index0  = {'headers',iA(cA)};
		s = spm_bch_headers;
    
	    %-----------------------------------------------------------
	    case 'means'
		BCH.index0  = {'means',iA(cA)};
		s = spm_means;
	    
	    %-----------------------------------------------------------
	    case 'realign'
		BCH.index0  = {'realign',iA(cA)};
		spm_realign_ui;
    
	    %-----------------------------------------------------------
	    case 'coreg'
		BCH.index0  = {'coreg',iA(cA)};
		spm_coreg_ui;
    
	    %-----------------------------------------------------------
	    case 'normalize'
		BCH.index0  = {'normalize',iA(cA)};
		spm_sn3d;
    
	    %-----------------------------------------------------------
	    case 'smooth'
		BCH.index0  = {'smooth',iA(cA)};
		spm_smooth_ui;    
    
	    %-----------------------------------------------------------
	    otherwise
		warning(sprintf('unknown type of analyse %s',typeA{cA}))
    
	end %-  switch typeA{cA}

end %- for cA = 1:length(iA)

cd(cwd);
