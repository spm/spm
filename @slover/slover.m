function [o, others] = slover(params, others, varargin)
% class constructor for slice overlay (slover) object
% FORMAT [o, others] = slover(params, others, varargin)
%
% Inputs
% params    - either:
%             - action string implementing class methods (see below)
%             - array of image names / vol structs to display
%             - structure with some fields for object (see below)
% others    - structure, containing extra fields for object (or children)
% varargin  - maybe some other parameters for action calls (see below)
%
% Outputs
% o       - slover object
% others  - any unrecognized fields from params, others
%
% Object fields are:
%  - img - array of structs with information for images to display
%        - img structs contain fields
%             type - one of {'truecolour' 'split', 'contour'};
%                   truecolour - displays transparent (see prop) image
%                      overlaid with any previous
%                   split - in defined area, replaces image present (SPM
%                      type activation display)
%                   contour - contour map of image overlaid.  See help
%                      for contours function in matlab
%             vol - vol struct info (see spm_vol)
%                   can also be vol containing image as 3d matrix
%                   set with add_blobs method
%             cmap - colormap for this image
%             nancol - color for NaN. If scalar, this is an index into
%                    the image cmap.  If 1x3 vector, it's a colour
%             prop - proportion of intensity for this cmap/img
%             func - function to apply to image before scaling to cmap
%                    (and therefore before min/max thresholding. E.g. a func of
%                    'i1(i1==0)=NaN' would convert zeros to NaNs
%             range - 2x1 vector of values for image to distribute colormap across
%                    the first row of the colormap applies to the first
%                    value in 'range', and the last value to the second
%                    value in 'range'
%             outofrange - behavior for image values to the left and
%                    right of image limits in 'range'.  Left means
%                    colormap values < 1, i.e for image values <
%                    range(1), if (range(1)<range(2)), and image values >
%                    range(1) where (range(1)>range(2)). If missing,
%                    display min (for Left) and max (for Right) value from colormap.
%                    Otherwise should be a 2 element cell array, where
%                    the first element is the colour value for image values
%                    left of 'range', and the second is for image values
%                    right of 'range'.  Scalar values for
%                    colour index the colormap, 3x1 vectors are colour
%                    values.  An empty array attracts default settings
%                    appropriate to the mode - i.e. transparent colour (where
%                    img(n).type is truecolour), or split colour.  Empty cells
%                    default to 0. 0 specifies that voxels with this
%                    colour do not influence the image (split =
%                    background, true = black)
%            hold  - resampling order for image (see spm_sample_vol) -
%                    default 1
%            background - value when resampling outside image - default
%                    NaN
%            linespec - string, applies only to contour map,
%                    e.g. 'w-' for white continuous lines
%            contours - vector, applies to contour map only, defines
%                    values in image for which to show contours
%                    (see help contours)
%            linewidth - scalar, width in points of contour lines
%
% - transform - either - 4x4 transformation to apply to image slice position,
%             relative to mm given by slicedef, before display
%               or     - text string, one of axial, coronal, sagittal
%                        These orientations assume the image is currently
%                        (after its mat file has been applied) axially
%                        oriented
% - slicedef - 2x3 array specifying dimensions for slice images in mm
%             where rows are x,and y of slice image, and cols are neg max dim,
%             slice separation and pos max dim
% - slices   - vector of slice positions in mm in z (of transformed image)
% - figure    - figure handle for slice display figure
%               The object used for the display is attached as 'UserData'
%               to this figure
% - figure_struct - stored figure parameters (in case figure dies and
%               needs to be recreated)
% - refreshf  - flag - if set or empty, refresh axis info for figure
%             else assume this is OK
% - clf       - flag, non zero -> clear figure before display.  Redundant
%               if refreshf == 0
% - resurrectf - if not zero, and figure (above) does not exist, will
%               attempt to recreate figure with same area properties.
%               Otherwise painting will give an error.
% - userdata  - flag, non zero -> attaches object to figure when ploting,
%               for use by callbacks (default is 1)
% - area      - struct with fields
%                  position - bottom left, x size y size 1x4 vector of
%                      area in which to display slices
%                  units    - one of
%                    inches,centimeters,normalized,points,{pixels}
%                  halign - one of left,{center},right
%                  valign - one of top,{middle},bottom
% - xslices  - no of slices to display across figure (defaults to an optimum)
% - cbar      - if empty, missing, no colourbar.  If an array of integers, then
%             indexes img array, and makes colourbar for each cmap for
%             that img.  Cbars specified in order of appearance L->R
% - labels - struct can be:
%                  - empty (-> default numerical labels)
%                  - 'none' (string) (no labels)
%                  - or contain fields:
%                  colour - colour for label text
%                  size - font size in units normalized to slice axes
%                  format - if = cell array of strings =
%                  labels for each slice in Z.  If is string, specifies
%                  sprintf format string for labelling in distance of the
%                  origin (Xmm=0, Ymm=0) of each slice from plane containing
%                  the AC, in mm, in the space of the transformed image
% - callback - callback string for button down on image panels. THe
%              following examples assume that you have the 'userdata'
%              field set to 1, giving you access to underlying object
%              To print to the matlab window the equivalent position in
%              mm of the position of a mouse click on one of the image
%              slices, set callback to:
%                   'get_pos(get(gcf, ''UserData''))'
%              To print the intensity values of the images at the clicked point:
%                   ['so_obj = get(gcf, ''UserData''); ' ...
%                    'point_vals(so_obj, get_pos(so_obj))']
% - printstr - string for printing slice overlay figure window, e.g.
%              'print -dpsc -painters -noui' (the default)
% - printfile - name of file to print output to; default 'slices.ps'
%
% Action string formats:
% FORMAT [cmap warnstr] = slover('getcmap', cmapname)
% Gets colormap named in cmapname string
%
% FORMAT [mx mn] = slover('volmaxmin', vol)
% Returns maximum and minimum finite values from vol struct 'vol'
%
% FORMAT vol = slover('blobs2vol', XYZ, vals, mat)
% returns (pseudo) vol struct for 3d blob volume specified
% in matrices as above
%
% FORMAT vol = slover('matrix2vol', mat3d, mat)
% returns (pseudo) vol struct for 3d matrix
% input matrices as above
%
% FORMAT obj = slover('basic_ui' [,dispf])
% Runs basic UI to fetch some parameters, does display, returns object
% If optional dispf parameter = 0, supresses display
%__________________________________________________________________________

% Matthew Brett
% $Id: slover.m 6623 2015-12-03 18:38:08Z guillaume $

myclass = 'slover';

% Default object structure
defstruct = struct('img', [], ...
    'transform', 'axial', ...
    'slicedef', [], ...
    'slices', [], ...
    'figure', [], ...
    'figure_struct', [], ...
    'refreshf', 1, ...
    'clf', 1, ...
    'resurrectf', 1, ...
    'userdata', 1, ...
    'area', [], ...
    'xslices', [], ...
    'cbar', [], ...
    'labels', [], ...
    'callback', ';', ...
    'printstr', 'print -dpsc -painters -noui', ...
    'printfile', 'slices.ps');

if nargin < 1
    o = class(defstruct, myclass);
    others = [];
    return
end
if nargin < 2
    others = [];
end

% parse out string action calls (class functions)
if ischar(params)
    switch params
        case 'getcmap'
            if nargin < 2
                error('Need colormap name');
            end
            o = pr_getcmap(others);
            return
        case 'volmaxmin'
            if nargin < 2
                error('Need volume to calculate max/min');
            end
            [o,others] = pr_volmaxmin(others);
            return
        case 'blobs2vol'
            if nargin < 4
                error('Need XYZ, vals, mat');
            end
            o = pr_blobs2vol(others, varargin{:});
            return
        case 'matrix2vol'
            if nargin < 3
                error('Need matrix and mat');
            end
            o = pr_matrix2vol(others, varargin{:});
            return
        case 'basic_ui'
            o = pr_basic_ui(others, varargin{:});
            if ~isempty(o), o = paint(o); end
            return
    end
    
    % if not action string, must be filename(s)
    params = spm_vol(params);
end

% Could these just be image vol structs?
if isfield(params, 'fname')
    for i = 1:numel(params)
        obj.img(i).vol = params(i);
    end
    params = obj;
end

% Deal with passed objects of this (or child) class
if isa(params, myclass)
    o = params;
    % Check for simple form of call
    if isempty(others), return, end
    
    % Otherwise, we are being asked to set fields of object
    [p,others] = mars_struct('split', others, defstruct);
    o = mars_struct('ffillmerge', o, p);
    return
end

% fill params with defaults, parse into fields for this object, children
params = mars_struct('fillafromb', params, others);
[params, others] = mars_struct('ffillsplit', defstruct, params);

% set the slover object
o  = class(params, myclass);

% refill with defaults
o = fill_defaults(o);
