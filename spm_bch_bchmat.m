function bch_mat = spm_bch_bchmat(bch_mfile,typeA)
% SPM batch system: Extract batch definitions from Mat file
% FORMAT bch_mat = spm_bch_bchmat(bch_mfile,typeA)
%
% bch_mfile  - 
% typeA      - Type of analysis
% bch_mat    - 
%
%_______________________________________________________________________
%
% This function gets a m-file that describes the bach process and
% returns a .mat file containing the variables used during the batch.
%_______________________________________________________________________
% %W% Jean-Baptiste Poline & Stephanie Rouquette %E%


%=======================================================================
% Programmers guide: Adding a new type of analysis
%=======================================================================
% When a new type of analysis is added, this file should be edited to
% define the variables that will be used by spm_input in BCH mode  AND
% add them in the 'bch_names' variable.  If  you wish to check the
% consistency of the variables (errors in the m-file specifications),
% this should be done here and NOT in the function called in spm_bch to
% do the actual job.


% temporarily goes to the m file dir and comes back at the end
%-----------------------------------------------------------------------
cwd = pwd;
try 
	cd(spm_str_manip(bch_mfile,'H'));
catch
	error([spm_str_manip(bch_mfile,'H') 'doesn''t exist']);
end

%-check that the m-file exists
%-----------------------------------------------------------------------
if ~(exist(spm_str_manip(bch_mfile,'rt')) == 2) 
   error(sprintf('%s not a m file in path', ...
	 spm_str_manip(bch_mfile,'rt')));
end

%-try to evaluate m-file
%-----------------------------------------------------------------------
try 
  	feval(spm_str_manip(bch_mfile,'rt'));
catch
	str = sprintf('cannot eval %s',spm_str_manip(bch_mfile,'rt'));
  	error([str ' : check your paths or parameter file']);
end

cd(cwd); %- come back here 
clear cwd 



switch typeA, case 'analyses'
%=======================================================================
  bch_names = {typeA,'type','work_dir','mfile'};


case 'defaults_edit'
%=======================================================================
  bch_names = {typeA, 'type_area' ...
	     'Misc','Printing','Hdr','Statistics', ...
	     'Normalisation','RealignCoreg','Reset'};


case {'headers','means','normalize','smooth','coreg','slicetime'}
%=======================================================================
  bch_names = {typeA};


case 'realign' 
%=======================================================================
  bch_names = {typeA,'sessions'};



case 'model'
%=======================================================================

% missing fields added here for model       	  
%-----------------------------------------------------------------------
% Added and modified fields are 
%   nodel(k).nscans
%   model(k).remain
%   model(k).interp (if not defined before)
%   model(k).time_sampl (if not defined before)
%-----------------------------------------------------------------------
%
if exist('model') == 1, 

  %- 
  if ~isfield(model,'time_sampl')
      for k =1:length(model)
          model(k).time_sampl = cell(1,model(k).nsess);
      end
      [model.interp] = deal('');
  end 

  for k = 1:length(model)
     for l = 1:model(k).nsess
        model(k).nscans(l) = model(k).nscans(l)+length(model(k).time_sampl{l});
        model(k).remain{l} = setdiff(1:model(k).nscans(l),...
                             model(k).time_sampl{l});     
     end %- for k = 1:length(model)
  end %- if ~isfield(model(k),'time_sampl')


end % if exist('model') == 1, 


bch_names = {typeA,'conditions','stochastics',...
  'regressors','parametrics','bf_ev','bf_ep'};



case 'contrasts'
%=======================================================================
% missing fields added here for contrasts
%-----------------------------------------------------------------------
% Added fields are 
%   contrasts(k).set_action
%-----------------------------------------------------------------------

if exist('contrasts') == 1
   for k = 1:length(contrasts)
      if ~isfield(contrasts(k),'set_action')
       	ncont = length(contrasts(k).names);
       	if ncont
              tmp = cell(1,ncont);
              [tmp{:}] = deal('c');
              contrasts(k).set_action = tmp;
              clear tmp;
       	end
      end
   end
end %- if exist('contrasts') == 1

bch_names = {typeA};



otherwise
%=======================================================================
   warning(sprintf('unknown type of analyse %s',typeA))

end
%=======================================================================







%-create bch_names with existing variables only
%-----------------------------------------------------------------------
str = {};
names = {};

for i=1:length(bch_names), 
   if exist(bch_names{i}) == 1,
       str = {str{:},bch_names{i}};
       names{length(names)+1} = bch_names{i};
   end
end

bch_names = names;
str = {str{:},'bch_names'};

% .mat saved here: here bch_mat does contain the path of the 
% bch_mfile (because we are in the working dir).
%-----------------------------------------------------------------------
save(spm_str_manip(bch_mfile,'rp'),str{:});

bch_mat = spm_str_manip(bch_mfile,'rp'); %-here contains the path

