function spm_save(f,var,varargin)
% Save text and numeric data to file
% FORMAT spm_save(f,var,opts,...)
% f     - filename (can be gzipped) {csv,tsv,json}
% var   - data array or structure
% opts  - optional inputs to be passed on to lower level function
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_save.m 7354 2018-06-22 10:44:22Z guillaume $


ext = lower(spm_file(f,'ext'));
switch ext
    
    case 'gz'
        fgz = spm_file(f,'basename');
        spm_save(fgz,var,varargin{:});
        gzip(fgz);
        delete(fgz);
        
    case {'csv','tsv'}
        if strcmp(ext,'csv')
            delim = ',';
        else
            delim = sprintf('\t');
        end
        if isstruct(var) || iscell(var) || isnumeric(var) || islogical(var)
            if isstruct(var)
                fn = fieldnames(var);
                var = struct2cell(var)';
                for i=1:numel(var)
                    var{i} = var{i}(:);
                    if ~iscell(var{i}), var{i} = cellstr(num2str(var{i},16)); end
                end
                var = [fn'; var{:}];
            elseif iscell(var)
                var = cellfun(@(x) num2str(x,16), var, 'UniformOutput',false);
            elseif isnumeric(var) || islogical(var)
                var = num2cell(var);
                var = cellfun(@(x) num2str(x,16), var, 'UniformOutput',false);
            end
            try, var = strtrim(var); end
            
            fid = fopen(f,'Wt');
            if fid == -1
                error('Unble to write file %s.', f);
            end
            for i=1:size(var,1)
                for j=1:size(var,2)
                    if isempty(var{i,j})
                        var{i,j} = 'n/a';
                    elseif any(var{i,j} == delim)
                        var{i,j} = ['"' var{i,j} '"'];
                    end
                    fprintf(fid,'%s',var{i,j});
                    if j < size(var,2), fprintf(fid,delim); end
                end
                fprintf(fid,'\n');
            end
            fclose(fid);
        elseif isa(var,'table')
            writetable(var,f,'FileType','text','Delimiter',delim);
        else
            error('Unknown data type.');
        end
        
    case 'json'
        if nargin < 3, varargin = {struct([])}; end
        spm_jsonwrite(f,var,varargin{:});
        
    otherwise
        error('Unknown file format.');
        
end
