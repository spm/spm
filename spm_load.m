function x = spm_load(f,v,hdr)
% Load text and numeric data from file
% FORMAT x = spm_load(f,v,hdr)
% f   - filename (can be gzipped) {txt,mat,csv,tsv,json,npy}
% v   - name of field to return if data stored in a structure [default: '']
%       or index of column if data stored as an array
% hdr - detect the presence of a header row for csv/tsv [default: true]
%
% x   - corresponding data array or structure
%__________________________________________________________________________
% Copyright (C) 1995-2019 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_load.m 7572 2019-04-12 16:16:32Z guillaume $


%-Get a filename if none was passed
%--------------------------------------------------------------------------
if ~nargin
    [f,sts] = spm_select(1,{...
        'mat',...                        % *.txt, *.mat
        '^.*\.csv$','^.*\.csv.gz$',...   % *.csv, *.csv.gz
        '^.*\.tsv$','^.*\.tsv.gz$',...   % *.tsv, *.tsv.gz
        '^.*\.json$','^.*\.json.gz$',... % *.json, *.json.gz
        '^.*\.npy$','^.*\.npz$',...      % *.npy, *.npz
        });
    if ~sts, x = []; return; end
end

if ~exist(f,'file')
    error('Unable to read file ''%s''',f);
end

if nargin < 2, v = ''; end
if nargin < 3, hdr = true; end % Detect

%-Load the data file
%--------------------------------------------------------------------------
switch spm_file(f,'ext')
    case 'txt'
        x = load(f,'-ascii');
    case 'mat'
        x  = load(f,'-mat');
    case 'csv'
        % x = csvread(f); % numeric data only
        x = dsvread(f,',',hdr);
    case 'tsv'
        % x = dlmread(f,'\t'); % numeric data only
        x = dsvread(f,'\t',hdr);
    case 'json'
        x = spm_jsonread(f);
    case 'gz'
        fz  = gunzip(f,tempname);
        sts = true;
        try
            x   = spm_load(fz{1});
        catch
            sts = false;
        end
        delete(fz{1});
        rmdir(spm_file(fz{1},'path'));
        if ~sts, error('Cannot load ''%s''.',f); end
    case 'npy'
        x = npyread(f);
    case 'npz'
        fz  = unzip(f,tempname);
        sts = true;
        try
            x = spm_load(fz{1});
            if numel(fz) > 1
                warning('Multiple NumPy arrays found and ignored.');
            end
        catch
            sts = false;
        end
        delete(fz{:});
        rmdir(spm_file(fz{1},'path'));
        if ~sts, error('Cannot load ''%s''.',f); end
        
    otherwise
        try
            x = load(f);
        catch
            error('Unknown file format.');
        end
end

%-Return relevant subset of the data if required
%--------------------------------------------------------------------------
if isstruct(x)
    if isempty(v)
        fn = fieldnames(x);
        if numel(fn) == 1 && isnumeric(x.(fn{1}))
            x = x.(fn{1});
        end
    else
        if ischar(v)
            try
                x = x.(v);
            catch
                error('Data do not contain array ''%s''.',v);
            end
        else
            fn = fieldnames(x);
            try
                x = x.(fn{v});
            catch
                error('Invalid data index.');
            end
        end
    end
elseif isnumeric(x)
    if isnumeric(v)
        try
            x = x(:,v);
        catch
            error('Invalid data index.');
        end
    elseif ~isempty(v)
        error('Invalid data index.');
    end
end


%==========================================================================
% function x = dsvread(f,delim)
%==========================================================================
function x = dsvread(f,delim,header)
% Read delimiter-separated values file into a structure array
%  * header line of column names will be used if detected
%  * 'n/a' fields are replaced with NaN

%-Input arguments
%--------------------------------------------------------------------------
if nargin < 2, delim = '\t'; end
if nargin < 3, header = true; end % true: detect, false: no
delim = sprintf(delim);
eol   = sprintf('\n');

%-Read file
%--------------------------------------------------------------------------
S   = fileread(f); % spm_file(f,'local','content');
if isempty(S), x = []; return; end
if S(end) ~= eol, S = [S eol]; end
S   = regexprep(S,{'\r\n','\r','(\n)\1+'},{'\n','\n','$1'});

%-Get column names from header line (non-numeric first line)
%--------------------------------------------------------------------------
h   = find(S == eol,1);
hdr = S(1:h-1);
var = regexp(hdr,delim,'split');
N   = numel(var);
n1  = isnan(cellfun(@str2double,var));
n2  = cellfun(@(x) strcmpi(x,'NaN'),var);
if header && any(n1 & ~n2)
    hdr     = true;
    try
        var = genvarname(var);
    catch
        var = matlab.lang.makeValidName(var,'ReplacementStyle','hex');
        var = matlab.lang.makeUniqueStrings(var);
    end
    S       = S(h+1:end);
else
    hdr     = false;
    fmt     = ['Var%0' num2str(floor(log10(N))+1) 'd'];
    var     = arrayfun(@(x) sprintf(fmt,x),(1:N)','UniformOutput',false);
end

%-Parse file
%--------------------------------------------------------------------------
if strcmpi(spm_check_version,'octave') % bug #51093
    S = strrep(S,delim,'#');
    delim = '#';
end
if ~isempty(S)
    d = textscan(S,'%s','Delimiter',delim);
else
    d = {[]};
end
if rem(numel(d{1}),N), error('Varying number of delimiters per line.'); end
d = reshape(d{1},N,[])';
allnum = true;
for i=1:numel(var)
    sts = true;
    dd = zeros(size(d,1),1);
    for j=1:size(d,1)
        if strcmp(d{j,i},'n/a')
            dd(j) = NaN;
        else
            dd(j) = str2double(d{j,i}); % i,j considered as complex
            if isnan(dd(j)), sts = false; break; end
        end
    end
    if sts
        x.(var{i}) = dd;
    else
        x.(var{i}) = d(:,i);
        allnum     = false;
    end
end

if ~hdr && allnum
    x = struct2cell(x);
    x = [x{:}];
end


%==========================================================================
% function x = npyread(f)
%==========================================================================
function x = npyread(f)
% Read data stored in NumPy NPY format
% https://www.numpy.org/devdocs/reference/generated/numpy.lib.format.html

[fid, msg] = fopen(f,'r','ieee-le');
if fid == -1
    error(msg);
end
%-magic string: \x93NUMPY
magic = fread(fid,[1,6],'*uint8');
if ~isequal(magic,[147 'NUMPY'])
    fclose(fid);
    error('Not an NPY formatted file.');
end
%-major and minor version numbers
ver = fread(fid,2,'*uint8');
if ~ismember(ver(1),[1 2]) || ~ismember(ver(2),0)
    warning('Unsupported version %d.%d of the NPY format.',ver);
end
%-length of the header data
if ver(1) == 1
    len = fread(fid,1,'*uint16');
else
    len = fread(fid,1,'*uint32');
end
%-header data
hdr = deblank(fread(fid,[1,len],'*char'));
hdr = regexp(hdr,'''(\w*)''\s*:\s*(''[<>|=\w]*''|[\(]{1}[\d\s,]*[\)]{1}|\w*)','tokens');
dt  = containers.Map(...
    {'b1','u1','u2','u4','u8','i1','i2','i4','i8','f4','f8'},...
    {'logical','uint8','uint16','uint32','uint64','int8','int16','int32','int64','single','double'});
fmt = struct;
for i=1:numel(hdr)
    switch hdr{i}{1}
        case 'descr'
            fmt.descr = strrep(hdr{i}{2},'''','');
            sb = fmt.descr(1) == '>';
            if ~isKey(dt,fmt.descr(2:end))
                error('Unknown or unsupported data type "%s".',fmt.descr(2:end));
            end
            dt = dt(fmt.descr(2:end));
            fmt.descr = struct('sb',sb,'dt',dt);
        case 'fortran_order'
            fmt.fortran_order = strcmp(hdr{i}{2},'True');
        case 'shape'
            fmt.shape = str2num(strrep(strrep(hdr{i}{2},'(','['),')',']'));
            if numel(fmt.shape) == 1, fmt.shape(2) = 1; end
        otherwise
            warning('Ignoring nknown disctionary key "%s".',hdr{i}{1});
    end
end
if numel(fieldnames(fmt)) ~= 3
    fclose(fid);
    error('Incomplete header data.');
end
%-data
x = fread(fid,prod(fmt.shape),['*' fmt.descr.dt]);
if fmt.descr.sb, x = swapbytes(x); end
if fmt.fortran_order
    x = reshape(x,fmt.shape);
else
    x = permute(reshape(x,fliplr(fmt.shape)),numel(fmt.shape):-1:1);
end
fclose(fid);
