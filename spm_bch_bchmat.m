function bch_mat = spm_bch_bchmat(bch_mfile,typeA)
%---------------------------------------------------------------
%
% This function gets a m-file that describes the 
% bach process and returns a .mat file containing 
% the variables used during the batch.
%---------------------------------------------------------------
% %W% Jean-Baptiste Poline & Stephanie Rouquette %E%

%---------------------------------------------------------------
% Programmers guide  :  when a new type of analysis is added, 
% this file should be edited to define the  
% variables that will be used by spm_input in BCH mode  AND 
% add them in the 'bch_names' variable.
% If  you wish to check the consistency of the 
% variables (errors in the m-file specifications), this 
% should be done here and NOT in the function called in spm_bch
% to do the actual job. 


%---------------------------------------------------------------
% temporarily goes to the m file dir and comes back at the end

cwd = pwd;
try 
	cd(spm_str_manip(bch_mfile,'H'));
catch
	error([spm_str_manip(bch_mfile,'H') 'doesn''t exist']);
end

%-------------------- check that the m-file exists -------------
if ~(exist(spm_str_manip(bch_mfile,'rt')) == 2) 
   error(sprintf('%s not a m file in path', ...
	 spm_str_manip(bch_mfile,'rt')));
end

%-------------------- try to evaluate m-file -------------------
try 
  	feval(spm_str_manip(bch_mfile,'rt'));
catch
	str = sprintf('cannot eval %s',spm_str_manip(bch_mfile,'rt'));
  	error([str ' : check your paths or parameter file']);
end

cd(cwd); %- come back here 
clear cwd 


%===============================================================
switch typeA

%===============================================================
case 'analyses'

  bch_names = {typeA,'type','work_dir','mfile'};

%===============================================================
case 'defaults_edit'

  bch_names = {typeA, 'type_area' ...
	     'Misc','Printing','Hdr','Statistics', ...
	     'Normalisation','RealignCoreg','Reset'};

%===============================================================
case {'headers','means','normalize','smooth'}

  bch_names = {typeA};

%===============================================================
case 'realign' 
  bch_names = {typeA,'sessions'};


%===============================================================
case 'model'


% missing fields added here for model       	  
%---------------------------------------------------------------
% Added fields are 
%   model(k).nsess
%   model(k).nscans
%   model(k).remain
%   model(k).interp
%---------------------------------------------------------------
%
if exist('model') == 1, 

   for k = 1:length(model)
      model(k).nsess = length(model(k).conditions);
   end

   if ~isfield(model(k),'time_sampl')
      for k = 1:length(model)
         for l = 1:model(k).nsess

   	    try %- if it can be evaluated...
   	         files = eval(model(k).files{l});
   	         e_nscans = size(files,1);
   	    catch
         	 e_nscans = size(model(k).files{l},1);
   	    end

            model(k).nscans(l) = e_nscans;
            model(k).remain{l} = 1:model(k).nscans(l);
            model(k).time_sampl{l} = [];
            model(k).interp = '';
         
          end %- for l = 1:model(k).nsess

      end %- for k = 1:length(model)
   end %- if ~isfield(model(k),'time_sampl')

   for k = 1:length(model)
      for l = 1:model(k).nsess
	   
	   try %- if it can be evaluated...
	      files = eval(model(k).files{l});
	      e_nscans = size(files,1);
	   catch
	      e_nscans = size(model(k).files{l},1);
	   end
      	   model(k).nscans(l) = e_nscans+length(model(k).time_sampl{l});
           model(k).remain{l} = setdiff(1:model(k).nscans(l),...
                                     model(k).time_sampl{l});
      
      end % for l = 1:model(k).nsess
   end % for k = 1:length(model)

end % if exist('model') == 1, 


% create regressors, parametrics and stochastics if necessary
%------------------------------------------------------------
if exist('model') == 1, for k = 1:length(model)

  if isempty(model(k).regressors) | ~any(model(k).regressors) ...
     | ~(exist('regressors')==1) %- not a variable in WS

	if (exist('regressors')==1) & length(regressors),
	   len = length(regressors);
	   regressors(len+1) = regressors(len);
	   regressors(len+1).number = 0;
	   model(k).regressors = (len+1)*ones(1,model(k).nsess);
	else
     	   regressors(1)  	= struct('number', 0);
	   model(k).regressors 	= ones(1,model(k).nsess);
	end
  end

  if isempty(model(k).parametrics) | ~any(model(k).parametrics) ...
     | ~(exist('parametrics')==1) %- not a variable in WS

	if (exist('parametrics')==1) & length(parametrics),
	   len = length(parametrics);
	   parametrics(len+1) = parametrics(len);
	   parametrics(len+1).type = 'none';
	   model(k).parametrics = (len+1)*ones(1,model(k).nsess);
	else
     	   parametrics(1)  	= struct('type', 'none');
	   model(k).parametrics = ones(1,model(k).nsess);
	end
  end

  if isempty(model(k).stochastics) | ~any(model(k).stochastics) ...
     | ~(exist('stochastics')==1) %- not a variable in WS

	if (exist('stochastics')==1) & length(stochastics),
	   len = length(stochastics);
	   stochastics(len+1) = stochastics(len);
	   stochastics(len+1).specify = 0;
	   model(k).stochastics = (len+1)*ones(1,model(k).nsess);
	else
     	   stochastics(1)  	= struct('specify', 0);
	   model(k).stochastics = ones(1,model(k).nsess);
	end
  end

end, end % if exist('model') == 1, for k = 1:length(model)


% missing fields added here for conditions       	  
%------------------------------------------------------------
% Added fields are 
%   conditions(k,l).number
%   conditions(k,l).fix_var_SOA
%   conditions(k,l).types
%------------------------------------------------------------

if exist('model') == 1, 

   if exist('conditions') == 1, for l = 1:length(conditions)

       conditions(l).number = length(conditions(l).names);	
       conditions(l).fix_var_SOA = 'Variable';
       conditions(l).types = cell(1,conditions(l).number);

       i_ev = find(conditions(l).bf_ev);
       i_ep = find(conditions(l).bf_ep); 
	    
       if ~isempty(i_ev)
          [conditions(l).types{i_ev}] = deal('events'); 
       end
       if ~isempty(i_ep)
          [conditions(l).types{i_ep}] = deal('epochs'); 
       end       

  end, end %- if exist('conditions') == 1, for l = 1:length(conditions) 

end % if exist('model') == 1

%- clear variables ....
clear i_ev i_ep k l str


bch_names = {typeA,'conditions','stochastics',...
  'regressors','parametrics','bf_ev','bf_ep'};


%===============================================================
case 'contrasts'


% missing fields added here for contrasts
%---------------------------------------------------------------
% Added fields are 
%   contrasts(k).set_action
%---------------------------------------------------------------

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

%===============================================================
case 'realignment'

   warning(sprintf('not implemented %s',typeA))


%===============================================================
otherwise
   warning(sprintf('unknown type of analyse %s',typeA))


%===============================================================
end


%-------- create bch_names with existing variables only
str = '';
names = {};

for i=1:length(bch_names), 
   if exist(bch_names{i})
       str = [str ' ' bch_names{i}];
       names{length(names)+1} = bch_names{i};
   end
end

bch_names = names;
str = [str  ' bch_names'];

% .mat saved here: here bch_mat does contain the path of the 
% bch_mfile (because we are in the working dir).
%---------------------------------------------------------------
eval(['save ' spm_str_manip(bch_mfile,'rp') ' ' str]);

bch_mat = spm_str_manip(bch_mfile,'rp'); %- here contains the path
%---------------------------------------------------------------



%===============================================================
%===============================================================
% sub function to check the consistency of the variables
%===============================================================

function ok = sf_chk_model(model)
ok=1; return;

function ok = sf_chk_conditions(conditions)

for l = 1:length(conditions)
  i_ev = find(conditions(l).bf_ev);
  i_ep = find(conditions(l).bf_ep); 
  if (any(intersect(i_ev,i_ep)) | ...
     any(setdiff(1:conditions(l).number,union(i_ev,i_ep))))
     str = sprintf('conditions(%d): bad event/epoch specif',l);
     error([str ' : check your parameter file']);
  end
end
ok = 1; return;
