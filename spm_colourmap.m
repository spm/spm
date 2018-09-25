function varargout = spm_colourmap(varargin)
% Colourmap multi-function
% FORMAT map = spm_colourmap
% Return the colourmap of the current figure as a three-column matrix of
% RGB triplets (between 0.0 and 1.0).
%
% FORMAT [map =] spm_colourmap(map)
% Define a colourmap or set it to the current figure.
% map         - gray, hot, pink, jet, ...: built-in colourmaps {64 x 3}
%             - gray-hot, ...: creates a 'split' colourmap {128 x 3 matrix}
%               The lower half is a gray scale and the upper half is
%               selected colourmap. This colourmap is used for viewing
%               'rendered' SPMs on a PET, MRI or other background images.
% 
% FORMAT [map = ] spm_colourmap(effect[,map])
% Apply an effect to a colourmap then return it or apply it to the current
% figure.
% effect      - 'Invert'   - invert (flip) the colourmap
%               'Brighten' - call MATLAB's brighten with a beta of +0.2
%               'Darken'   - call MATLAB's brighten with a beta of -0.2
%
% FORMAT maps = spm_colourmap('list')
% Return the list of all colourmaps' name (see graph3d).
%
% FORMAT [map =] spm_colourmap('load',fname)
% Load a colourmap from file (*.lut, *.cmap, *.mat) then return it or apply
% it to the current figure.
%
% FORMAT spm_colourmap('save',fname[,map])
% Save a colourmap to file (format according to file extension).
%__________________________________________________________________________
%
% A repository of colourmaps with linearised luminance is available at:
%   https://github.com/CPernet/brain_colours
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_colourmap.m 7428 2018-09-25 13:25:13Z guillaume $


%-FORMAT map = spm_colourmap([map])
%--------------------------------------------------------------------------
if ~nargin
    varargout = { colormap };
    return;
elseif isnumeric(varargin{1}) && size(varargin{1},2) == 3
    map = varargin{1};
    action = 'NOP';
else
    action = varargin{1};
end

%-Colourmap operation
%--------------------------------------------------------------------------
N = 64; % or size(get(gcf,'colormap'),1)

switch lower(action)
    
    case 'load'
        if nargin > 1
            fname = varargin{2};
        else
            fname = spm_select(1,'any','Select colourmap file...');
        end
        map = load_colourmap(fname);
    
    case 'save'
        if nargin > 1
            fname = varargin{2};
        else
            [fname, pth] = uiputfile;
            if isequal(fname,0) || isequal(pth,0)
                return;
            end
            fname = fullfile(pth,fname);
        end
        map = spm_colourmap(varargin{3:end});
        save_colourmap(fname,map);
        return;
                
    case 'default'
        map = 'default';
    
    case 'list'
        varargout = { list_colourmaps };
        return;
        
    case list_colourmaps
        map = feval(lower(action),N);
        
    case 'gray-hot'
        tmp = hot(N + N/4);  tmp = tmp((1:N) + N/4,:);
        map = [gray(N); tmp];
        
    case 'gray-cool'
        if N == 64
            tmp = [zeros(10,1) zeros(10,1) linspace(0.5,1,10)';
                 zeros(31,1) linspace(0,1,31)' ones(31,1);
                 linspace(0,1,23)' ones(23,1) ones(23,1) ];
        else
            tmp = cool(N);
        end
        map = [gray(64); tmp];
        
    case 'gray-pink'
        tmp = pink(N + N/4); tmp = tmp((1:N) + N/4,:);
        map = [gray(N); tmp];
        
    case {'invert','brighten','darken','equalise'}
        map = spm_colourmap(varargin{2:end});
        switch lower(action)
            case 'invert'
                map = flipud(map);
            case 'brighten'
                map = brighten(map, 0.2);
            case 'darken'
                map = brighten(map, -0.2);
            case 'equalise'
                map = equalise_colourmap(map);
        end
        
    case 'nop'
        
    otherwise
        s = regexp(action,'-','split');
        if numel(s) > 1
            map = [];
            for i=1:numel(s)
                map = [map; spm_colourmap(s{i})];
            end
        else
            error('Illegal action specification');
        end
end

%-Update current figure's colourmap or return it
%--------------------------------------------------------------------------
if nargout
    varargout = { map };
else
    colormap(map);
end


%==========================================================================
function maps = list_colourmaps
%==========================================================================
maps = {...
    'parula', 'hsv', 'hot', 'gray', 'bone', 'copper', 'pink', 'white', ...
    'flag', 'lines', 'colorcube', 'jet', 'prism', 'cool', 'autumn', ...
    'spring', 'winter', 'summer'};


%==========================================================================
function map = load_colourmap(fname)
%==========================================================================
ext = spm_file(fname,'ext');
if isempty(ext), ext = 'mat'; fname = spm_file(fname,'ext',ext); end

switch ext
    case 'lut'
        fid = fopen(fname,'rb');
        if fid == -1, error('Cannot open colourmap file.'); end
        map = fread(fid,Inf,'uchar=>double');
        fclose(fid);
        map = reshape(map,[],3);
    case 'cmap'
        map = load(fname,'-ascii');
    case 'mat'
        map = struct2cell(load(fname));
        map = map{1};
    otherwise
        error('Unknown colour map format.');
end
if any(map(:)>1)
    map = map / 255;
end
if size(map, 2) ~= 3
    warning('Colourmap was not an N by 3 matrix.');
    map = [];
end


%==========================================================================
function save_colourmap(fname,map)
%==========================================================================
ext = spm_file(fname,'ext');
if isempty(ext), ext = 'mat'; fname = spm_file(fname,'ext',ext); end

switch ext
    case 'lut'
        fid = fopen(fname,'wb');
        if fid == -1, error('Cannot open colourmap file.'); end
        fwrite(fid,255*map,'uchar');
        fclose(fid);
    case 'cmap'
        save(fname,'map', '-ascii');
    case 'mat'
        save(fname,'map', spm_get_defaults('mat.format'));
    otherwise
        error('Unknown colourmap format.');
end


%==========================================================================
function map = equalise_colourmap(map)
%==========================================================================
% Require equalisecolourmap.m from Peter Kovesi:
%   https://www.peterkovesi.com/matlabfns/Colourmaps/equalisecolourmap.m
% and the Image Processing Toolbox (whitepoint, makecform, applycform)
% map = equalisecolourmap('RGB',map,'CIE76',[1 1 1],1/25*length(map),0,0);
error('Colour contrast equalisation over a colourmap not yet supported.');
