function [cmap, warnstr] = pr_getcmap(acmapname)
% Get colormap of name acmapname
% FORMAT [cmap, warnstr] = pr_getcmap(acmapname)
%
% Inputs
% acmapname   - string.  Can be (in order of precedence)
%               - matrix name in base workspace
%               - colour name; one of 'red','green','blue','cyan',
%                 'magenta', 'yellow', 'black', 'white'
%               - filename of .mat or .lut file.  If filename has no
%                 extension, assumes '.mat' extension
%
% Outputs
% cmap        - Nx3 colormap matrix
%               or empty if fails
% warnstr     - warning message if fails
%__________________________________________________________________________

% Matthew Brett
% $Id: pr_getcmap.m 6623 2015-12-03 18:38:08Z guillaume $

cmap = []; warnstr = [];
if nargin < 1
    acmapname = '';
end
if isempty(acmapname)
    warnstr = 'No colormap name passed';
    return
end
% try a matrix first
cmap = evalin('base',acmapname,'[]');
if ~isempty(cmap)
    if size(cmap, 2)~=3
        warnstr = ['Colormap matrix ' acmapname ' was not N by 3'];
        cmap = [];
    end
    return
end
% not a matrix, is it...
% a colour name?
tmp = strcmpi(acmapname, {'red','green','blue','cyan', 'magenta', ...
    'yellow', 'black', 'white'});
coldefs = [1 0 0;
    0 1 0;
    0 0 1;
    0 1 1;
    1 0 1;
    1 1 0;
    0 0 0;
    1 1 1];

if any(tmp)
    coldef = coldefs(tmp, :);
    if ~any(diff(coldef))
        cmap = coldef;
    else
        cmap = (0:63)' * coldef / 63;
    end
    return
end
% is it a file?
oname = acmapname;
[p,f,e] = fileparts(acmapname);
% if no extension, add .mat
if isempty(e)
    e = '.mat';
    acmapname = fullfile(p, [f e]);
end
ef = exist(acmapname, 'file');
% file doesn't exist? Try home directory of this mfile
if ~ef
    mfp = fileparts(which(mfilename));
    acmapname = fullfile(mfp, [f e]);
    ef = exist(acmapname, 'file');
end
if ~ef
    warnstr = ['No matrix or file ''' oname ''''];
    return
end
% found file, get cmap
switch lower(e)
    case '.mat'
        % try for matrix of same name
        s = load(acmapname);
        if isfield(s, f)
            cmap = getfield(s, f);
        else % get first matrix in mat file
            s = struct2cell(s);
            cmap = s{1};
        end
        if size(cmap, 2)~=3
            warnstr = ['Colormap from ' acmapname ' was not an N by 3 matrix'];
            cmap = [];
        end
    case '.lut'
        fid = fopen(acmapname, 'rb');
        if fid~=-1
            cmap = fread(fid, Inf);
            l = length(cmap);
            if ~rem(l,3)
                cmap = reshape(cmap, l/3, 3) / 255;
            else
                warnstr = ['LUT map ' acmapname ' was wrong size'];
            end
            fclose(fid);
        else
            warnstr = ['Cannot open LUT colormap file ' acmapname];
        end
    otherwise
        warnstr = ['Unrecognized file extension ' e ' for ' acmapname];
end
