function [BF, matlabbatch, output] = bf_wizard_output(S)
% A handy command-line based batch filler with some defaults for DAiSS
% output module, pick a few options, and it will default for unpopulated
% fields
%
% Current *definitely* supported output methods include:
%   - image_dics
%   - image_mv
%   - image_power
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


if ~isfield(S,'batch'), matlabbatch = []; else; matlabbatch = S.batch;  end
if ~isfield(S,'BF'),        error('I need a BF.mat file specified!');   end
if ~isfield(S,'method'),    error('You need to specify a method!');     end
if ~isfield(S,S.method),    S.(S.method) = struct();                    end
if ~isfield(S,'run'),       S.run = 1;                                  end

% specify BF, ensure its a cell...
if ~iscell(S.BF)
    S.BF = {S.BF};
end

output = struct();
output.BF = S.BF;

% there's a few common arguments which might be specified above the method
% options structure (erroneously) so spot them and put them in the right place;
targets = {'conditions','contrast','woi','foi'};
for ii = 1:numel(targets)
    if isfield(S,targets{ii})
        S.(S.method).(targets{ii}) = S.(targets{ii});
    end
end

% whatconditions structure is common to many methods, get it in order now
% to save coding multiple times later
if isfield(S.(S.method),'conditions')
    % check the response to conditions
    if sum(strcmp(S.conditions,'all'))
        S.(S.method).whatconditions.all = 1;
    else
        % check its a cell aray
        if ~iscell(S.(S.method).conditions)
            error('custom condition labels must be in a cell array!')
        end
        S.(S.method).whatconditions.condlabel = S.conditions;
    end
end

% check method actually exists
try
    opts = feval(['bf_output_' S.method]);
catch
    error('not a valid output method!')
end

output.plugin.(S.method) = struct();

% The output module is a lot more complicated than others, so special cases
% need to be implemented for each specific routine
switch S.method
    
    case 'image_dics'
        
        % check if a reference field is present
        if ~isfield(S.(S.method),'reference')
            warning('no reference suggestion found, assuming power image');
            S.(S.method).reference.power = 1;
        else
            if isfield(S.(S.method).reference,'refchan')
                if ~isfield(S.(S.method).reference.refchan,'name')
                    error('please enter the name of the refchan you want to use')
                end
                if ~isfield(S.(S.method).reference.refchan,'shuffle')
                    S.(S.method).reference.refchan.shuffle = 0;
                end
            end
        end
        
        % other boilerplate cases 
        if ~isfield(S.(S.method),'whatconditions')
            S.(S.method).whatconditions.all = 1;
        end
        
         if ~isfield(S.(S.method),'foi')
            error('plese specify a frequency band of interest!');
        end
        
        if ~isfield(S.(S.method),'woi')
            error('plese specify a frequency band of interest!');
        end
            
        if ~isfield(S.(S.method),'contrast')
            error('please specify a contrast!')
        end
           
    case 'image_mv'
        
        isdesign = struct();
        S.iscustom = 0;
        % Conditions can be specified in S.image_mv.conditions, but thats
        % not where the rest off code will anticipate it, so lets check.
        if isfield(S.(S.method),'whatconditions')
            S.iscustom = 1;
            isdesign.custom.whatconditions = S.(S.method).whatconditions;
        end
        
        if S.iscustom
            
            if isfield(S.(S.method),'woi')
                isdesign.custom.woi = S.(S.method).woi;
            else
                error('plese specify a frequency band of interest!');
            end
            
            if isfield(S.(S.method),'contrast')
                isdesign.custom.contrast = S.(S.method).contrast;
            else
                error('please specify a contrast!')
            end
            
        else % Its a location to a preloaded matrix
            if isfield(S.(S.method),'design')
                isdesign.design = S.(S.method).design;
            end
        end
        
        S.(S.method).isdesign = isdesign;
        
        % catch when foi isnt' specified
        if ~isfield(S.(S.method),'foi')
            error('plese specify a frequency band of interest!');
        end
        
    case 'image_power'
        
        if ~isfield(S.(S.method),'whatconditions')
            S.(S.method).whatconditions.all = 1;
        end
        
        if ~isfield(S.(S.method),'foi')
            error('plese specify a frequency band of interest!');
        end
        
        if ~isfield(S.(S.method),'woi')
            error('plese specify a time window of interest!');
        end
            
        if ~isfield(S.(S.method),'contrast')
            error('please specify a contrast!')
        end
        
end % switch S.method

for ii = 1:numel(opts.val)
    
    tag = opts.val{ii}.tag;
    val = opts.val{ii}.val;
    
    if ~isfield(S.(S.method),tag)
        output.plugin.(S.method).(tag) = val{1};
    else
        output.plugin.(S.method).(tag) = S.(S.method).(tag);
    end
    
end

% determine number of jobs in list then iterate by 1;
jobID = numel(matlabbatch) + 1;
% generate matlabbatch
matlabbatch{jobID}.spm.tools.beamforming.output = output;

if S.run
    out = spm_jobman('run',matlabbatch);
    BF = out{1,1}.BF{:};
else
    BF = [];
end