function spm_save(f,var,varargin)
% Save text and numeric data to file
% FORMAT spm_save(f,var,opts,...)
% f     - filename (can be gzipped) {csv,tsv,json,txt,mat,npy}
% var   - data array or structure
% opts  - optional inputs to be passed on to lower level function
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2018-2023 Wellcome Centre for Human Neuroimaging


ext = lower(spm_file(f,'ext'));
if isempty(ext), ext = ''; end
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
                    if ~iscell(var{i})
                        var{i} = cellstr(num2str(var{i},16));
                        var{i} = strtrim(var{i});
                        var{i}(cellfun(@(x) strcmp(x,'NaN'),var{i})) = {'n/a'};
                    end
                end
                var = [fn'; var{:}];
            elseif iscell(var) || isnumeric(var) || islogical(var)
                if isnumeric(var) || islogical(var)
                    var = num2cell(var);
                end
                var = cellfun(@(x) num2str(x,16), var, 'UniformOutput',false);
                var = strtrim(var);
                var(cellfun(@(x) strcmp(x,'NaN'),var)) = {'n/a'};
            end
            
            fid = fopen(f,'Wt');
            if fid == -1
                error('Unble to write file %s.', f);
            end
            for i=1:size(var,1)
                for j=1:size(var,2)
                    var{i,j} = strrep(var{i,j},'"','""');
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
        
    case {'','txt','md','html'}
        var = cellstr(var);
        fid = fopen(f,'Wt');
        if fid == -1
            error('Unable to write file %s.', f);
        end
        for i=1:numel(var)
            fprintf(fid,'%s\n',var{i});
        end
        fclose(fid);
        
    case 'mat'
        if nargin < 3, varargin = {spm_get_defaults('mat.format')}; end
        if ~isstruct(var)
%             if ~iscellstr([var, varargin])
%                 error('Variable name must be a text scalar.');
%             end
%             var = struct(var,evalin('caller',var));
%             for i=1:numel(varargin)
%                 if varargin{i}(1) ~= '-'
%                     var.(varargin{i}) = evalin('caller',varargin{i});
%                     varargin{i} = '';
%                 end
%             end
            
            var = struct(inputname(2),var);
            for i=1:numel(varargin)
                if varargin{i}(1) ~= '-'
                    name = inputname(i+2);
                    if isempty(name)
                        error('Variable name must be a text scalar.');
                    end
                    var.(name) = varargin{i};
                    varargin{i} = '';
                end
            end
            
            varargin(cellfun('isempty',varargin)) = [];
        end
        save(f,'-struct','var',varargin{:});
        
    case 'npy'
        fid = fopen(f,'W');
        fwrite(fid,[147 'NUMPY'],'uint8');
        fwrite(fid,[2 0],'uint8');
        dt  = containers.Map(...
            {'logical','uint8','uint16','uint32','uint64','int8','int16','int32','int64','single','double'},...
            {'b1','u1','u2','u4','u8','i1','i2','i4','i8','f4','f8'});
        hdr = ['{''descr'': ''<' dt(class(var)) ''', ' ...
            '''fortran_order'': True, ' ...
            '''shape'': (' sprintf('%d, ',size(var)) '), }' sprintf('\n')];
        len = 6+2+4+numel(hdr);
        hdr = [hdr blanks(ceil(len/8)*8 - len)];
        fwrite(fid,numel(hdr),'uint32');
        fwrite(fid,hdr,'uint8');
        fwrite(fid,var,class(var));
        fclose(fid);
        
    otherwise
        error('Unknown file format.');
        
end
