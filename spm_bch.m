function  spm_bch(varargin);
% FORMAT spm_bch(bch_mfile);
%---------------------------------------------------------------
% bch_mfile : 	m file containing the description of the batch
% 		analyses
%
% bch_mfile will contain a variable 'analyses' with fields 
% 
%  'type',  	 INDICES
%  'contrastes', INDICES
%  'model',      INDICES
%  'work_dir',   INDICES               
%
% and variables : 'type', 'contrastes', 'model', 'work_dir'
% See batch_doc.txt for indication on how to fill these
%
% 
% Jean-Baptiste Poline & Stephanie Rouquette 
%---------------------------------------------------------------

%---------------------------------------------------------------
% Programmer's guide :
%
%
%


%- Global batch variables
%---------------------------------------------------------------

global batch_mat;
global iA;

%- Inputs and Defaults
%---------------------------------------------------------------
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
   MODALITY = varargin{2};
   if ~ismember(MODALITY,{'FMRI','PET'}),
      fprintf('FORMAT : spm_bch(bch_mfile [,''FMRI''|''PET'']\n');
      return;
   end
end

spm('defaults',MODALITY);

%- m->mat file for the upper level analysis  
%---------------------------------------------------------------
try 
   ana_mat = spm_bch_bchmat(bch_mfile,'analyses');
   %--------- put it in memory 
   ana_mat = load(ana_mat);
catch
   error(sprintf('\nCan''t evaluate %s of ''analyses'' ', bch_mfile));
end

%- save current directory 
%---------------------------------------------------------------
cwd = pwd;

for cA = 1:length(spm_input_b('batch',ana_mat,{'analyses',1},'index'))

	% go back to current dir (useful if path are relative) 
	%---------------------------------------------------------------
	eval(['cd '  cwd]); 
	%- go into working directory 
	%---------------------------------------------------------------
	try 
	  wk_dir = spm_input_b('batch',ana_mat,{'analyses',1,'work_dir',cA},1);
	  eval(['cd ' wk_dir]);
	catch 
	   pwd
	   error(['can''t go to work dir '  wk_dir]); 
	end
	typeA 	= spm_input_b('batch',ana_mat,{'analyses',1,'type',cA},1);
   	mfileA 	= spm_input_b('batch',ana_mat,{'analyses',1,'mfile',cA},1);
	iA  	= spm_input_b('batch',ana_mat,{'analyses',1},'index',cA);

	%- m->mat file for sub analyses  
   	batch_mat = spm_bch_bchmat(mfileA,typeA);

	switch typeA

	    %---------------------------------------------------------------
	    case 'defaults_edit'  

		%------------ save iA
		iA_save = iA;

		%- get the default areas (number) to work on
		%-------------------------------------------
		areas = spm_input_b('batch',batch_mat, ...
				{'defaults_edit',iA_save},'type_area');

		%- get the indexes of these default areas
		%-------------------------------------------
		index = spm_input_b('batch',batch_mat, ...
				{'defaults_edit',iA_save},'index');

		for i_area = 1:length(areas)
		    str_area = spm_input_b('batch',batch_mat,...
				{'defaults_edit',iA_save,'type_area',i_area},1);
		    iA 	= index(i_area)
		    str_area
		    spm_defaults_edit(str_area);
	  	end

		%------------ restore iA 
		% iA = iA_save;


	    %---------------------------------------------------------------
	    case 'model'
		spm_fmri_spm_ui;
    
	    %---------------------------------------------------------------
	    case 'contrastes'
		s = spm_bch_GetCont(batch_mat,iA); 	
		s = spm_bch_DoCont(batch_mat); 			
    
	    %---------------------------------------------------------------
	    case 'headers'
		s = spm_bch_headers(batch_mat,iA)
    
	    %---------------------------------------------------------------
	    case 'means'
		s = spm_means(batch_mat,iA)
	    
	    %---------------------------------------------------------------
	    case 'realignment'
		spm_realign;
    
	    %---------------------------------------------------------------
	    case 'normalisation'
		spm_sn3d;
    
	    %---------------------------------------------------------------
	    case 'smooth'
		spm_smooth_ui;    
    
	    %---------------------------------------------------------------
	    otherwise
		warning(sprintf('unknown type of analyse %s',typeA))
    
	end %-  switch typeA

end %- for cA = 1:length

cd(cwd);

