function x = spm_load(f)
% Load text and numeric data from file
% FORMAT x = spm_load(f)
% f  - filename {txt,mat,csv,tsv,json}
% x  - corresponding data array or structure
%__________________________________________________________________________
% Copyright (C) 1995-2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_load.m 6624 2015-12-03 18:58:05Z guillaume $


%-Get a filename if none was passed
%--------------------------------------------------------------------------
if ~nargin
    [f,sts] = spm_select(1,{...
        'mat',...                        % *.txt, *.mat
        '^.*\.csv$','^.*\.csv.gz$',...   % *.csv, *.csv.gz
        '^.*\.tsv$','^.*\.tsv.gz$',...   % *.tsv, *.tsv.gz
        '^.*\.json$','^.*\.json.gz$',... % *.json, *.json.gz
        });
    if ~sts, x = []; return; end
end

if ~exist(f,'file')
    error('Unable to read file ''%s''',f);
end

%-Load the data file
%--------------------------------------------------------------------------
switch spm_file(f,'ext')
    case 'txt'
        x = load(f,'-ascii');
    case 'mat'
        x  = load(f,'-mat');
        fn = fieldnames(x);
        if numel(fn) == 1 && isnumeric(x.(fn{1}))
            x = x.(fn{1});
        end
    case 'csv'
        try
            x = csvread(f);
        catch
            x = dsvread(f,',');
        end
    case 'tsv'
        try
            x = dlmread(f,'\t');
        catch
            x = dsvread(f,'\t');
        end
    case 'json'
        x = spm_jsonread(f);
        % check if x only contains numeric and, if so, return array
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
    otherwise
        try
            x = load(f);
        catch
            error('Unknown file format.');
        end
end


%==========================================================================
function x = dsvread(f,delim)
% Read delimiter-separated values file containing a header line
% 'n/a' fields are replaced with NaN

if nargin < 2, delim = '\t'; end
delim = sprintf(delim);
eol = sprintf('\n');

S   = fileread(f);
if isempty(S), x = []; return; end
if S(end) ~= eol, S = [S eol]; end
S   = regexprep(S,{'\r\n','(\n)\1+'},{'\n','$1'});
h   = find(S == eol,1);
hdr = S(1:h-1);
var = regexp(hdr,delim,'split');
try
    var = genvarname(var);
catch
    var = matlab.lang.makeValidName(var);
    var = matlab.lang.makeUniqueStrings(var);
end
N   = numel(var);
S   = S(h+1:end);

d = textscan(S,'%s','Delimiter',delim);
if rem(numel(d{1}),N), error('Varying number of delimiters per line.'); end
d = reshape(d{1},N,[])';
for i=1:numel(var)
    sts = true;
    dd = zeros(size(d,1),1);
    for j=1:size(d,1)
        if strcmp(d{j,i},'n/a')
            dd(j) = NaN;
        else
            dd(j) = str2double(d{j,i});
            if isnan(dd(j)), sts = false; break; end
        end
    end
    if sts
        x.(var{i}) = dd;
    else
        x.(var{i}) = d(:,i);
    end
end
