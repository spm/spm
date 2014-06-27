function [varargout] = spm_uitable(varargin)
% WARNING: This feature is not supported in MATLAB
% and the API and functionality may change in a future release.

% UITABLE creates a two dimensional graphic uitable component in a figure window.
%     UITABLE creates a 1x1 uitable object using default property values in
%     a figure window.
%
%     UITABLE(numrows,numcolumns) creates a uitable object with specified
%     number of rows and columns.
%
%     UITABLE(data,columnNames) creates a uitable object with the specified
%     data and columnNames. Data can be a cell array or a vector and
%     columnNames should be cell arrays.
%
%     UITABLE('PropertyName1',value1,'PropertyName2',value2,...) creates a
%     uitable object with specified property values. MATLAB uses default
%     property values for any property not explicitly set. The properties
%     that user can set are: ColumnNames, Data, GridColor, NumColumns,
%     NumRows, Position, ColumnWidth and RowHeight.
%
%     UITABLE(figurehandle, ...) creates a uitable object in the figure
%     window specified by the figure handle.
%
%     HANDLE = UITABLE(...) creates a uitable object and returns its handle.
%
%     Properties:
%
%     ColumnNames:  Cell array of strings for column names.
%     Data:         Cell array of values to be displayed in the table.
%     GridColor:    string, RGB vector.
%     NumColumns:   int specifying number of columns.
%     NumRows:      int specifying number of rows.
%     Parent:       Handle to figure or uipanel. If not specified, it is gcf.
%     Position:     4 element vector specifying the position.
%     ColumnWidth:  int specifying the width of columns.
%     RowHeight:    int specifying the height of columns.
%
%     Enabled:      Boolean specifying if a column is enabled.
%     Editable:     Boolean specifying if a column is editable.
%     Units:        String - pixels/normalized/inches/points/centimeters.
%     Visible:      Boolean specifying if table is visible.
%     DataChangedCallback - Callback function name or handle.
%
%
%     Examples:
%
%     t = uitable(3, 2);
%
%     Creates a 3x2 empty uitable object in a figure window.
%
%     f = figure;
%     t = uitable(f, rand(5), {'A', 'B', 'C', 'D', 'E'});
%
%     Creates a 5x5 uitable object in a figure window with the specified
%     data and the column names.
%
%     data = rand(3);
%     colnames = {'X-Data', 'Y-Data', 'Z-Data'};
%     t = uitable(data, colnames,'Position', [20 20 250 100]);
%
%     Creates a uitable object with the specified data and column names and
%     the specified Position.
%
%     See also AWTCREATE, AWTINVOKE, JAVACOMPONENT, UITREE, UITREENODE

%   Copyright 2002-2006 The MathWorks, Inc.
%   $Revision: 6072 $  $Date: 2006/11/29 21:53:13 $

%   Release: R14. This feature will not work in previous versions of MATLAB.

% $Id: spm_uitable.m 6072 2014-06-27 16:35:30Z guillaume $

% Setup and P-V parsing

if isempty(varargin)
    if ~isempty(javachk('awt')) || spm_check_version('matlab','7.3') <= 0
        varargout{1} = 'off';
    else
        varargout{1} = 'on';
    end
    return
end

if ~isempty(javachk('awt')) || spm_check_version('matlab','7.3') <= 0
    varargout{1} = [];
    varargout{2} = [];
    return;
end

if ischar(varargin{1})
    switch varargin{1}
        case 'set'
            data = varargin{2};
            columnNames = varargin{3};
            [htable,hcontainer] = UiTable(data,columnNames);
            varargout{1} = htable;
            varargout{2} = hcontainer;
        case 'get'
            htable = varargin{2};
            columnNames = get(htable,'columnNames');
            nc = get(htable,'NumColumns');
            nr = get(htable,'NumRows');
            data = get(htable,'data');
            data2 = cell(nr,nc);
            for i=1:nc
                for j=1:nr
                    data2{j,i} = data(j,i);
                end
            end
            varargout{1} = data2;
            varargout{2} = columnNames;
    end
else
    [htable,hcontainer] = UiTable(varargin{:});
    varargout{1} = htable;
    varargout{2} = hcontainer;
end


function [table,container] = UiTable(varargin)

error(nargoutchk(0,2,nargout));

parent = [];
numargs = nargin;

datastatus=false; columnstatus=false;
rownum = 1; colnum = 1; % Default to a 1x1 table.
position = [20 20 200 200];
combo_box_found = false;
check_box_found = false;

import com.mathworks.hg.peer.UitablePeer;

if (numargs > 0 && isscalar(varargin{1}) && ishandle(varargin{1}) && ...
        isa(handle(varargin{1}), 'figure'))
    parent = varargin{1};
    varargin = varargin(2:end);
    numargs = numargs - 1;
end

if (numargs > 0 && isscalar(varargin{1}) &&  ishandle(varargin{1}))
    if ~isa(varargin{1}, 'javax.swing.table.DefaultTableModel')
        error('MATLAB:uitable:UnrecognizedParameter', ['Unrecognized parameter: ', varargin{1}]);
    end
    data_model = varargin{1};
    varargin = varargin(2:end);
    numargs = numargs - 1;

elseif ((numargs > 1) && isscalar(varargin{1}) && isscalar(varargin{2}))
    if(isnumeric(varargin{1}) && isnumeric(varargin{2}))
        rownum = varargin{1};
        colnum = varargin{2};

        varargin = varargin(3:end);
        numargs = numargs-2;
    else
        error('MATLAB:uitable:InputMustBeScalar', 'When using UITABLE numrows and numcols have to be numeric scalars.')
    end

elseif ((numargs > 1) && isequal(size(varargin{2},1), 1) && iscell(varargin{2}))
    if (size(varargin{1},2) == size(varargin{2},2))
        if (isnumeric(varargin{1}))
            varargin{1} = num2cell(varargin{1});
        end
    else
        error('MATLAB:uitable:MustMatchInfo', 'Number of column names must match number of columns in data');
    end
    data = varargin{1};     datastatus        = true;
    coln = varargin{1+1};   columnstatus      = true;

    varargin = varargin(3:end);
    numargs = numargs-2;
end

for i = 1:2:numargs-1
    if (~ischar(varargin{i}))
        error('MATLAB:uitable:UnrecognizedParameter', ['Unrecognized parameter: ', varargin{i}]);
    end
    switch lower(varargin{i})
        case 'data'
            if (isnumeric(varargin{i+1}))
                varargin{i+1} = num2cell(varargin{i+1});
            end
            data        = varargin{i+1};
            datastatus  = true;

        case 'columnnames'
            if(iscell(varargin{i+1}))
                coln            = varargin{i+1};
                columnstatus    = true;
            else
                error('MATLAB:uitable:InvalidCellArray', 'When using UITABLE Column data should be 1xn cell array')
            end

        case 'numrows'
            if (isnumeric(varargin{i+1}))
                rownum = varargin{i+1};
            else
                error('MATLAB:uitable:NumrowsMustBeScalar', 'numrows has to be a scalar')
            end

        case 'numcolumns'
            if (isnumeric(varargin{i+1}))
                colnum = varargin{i+1};
            else
                error('MATLAB:uitable:NumcolumnsMustBeScalar', 'numcolumns has to be a scalar')
            end

        case 'gridcolor'
            if (ischar(varargin{i+1}))
                gridcolor = varargin{i+1};
            else if (isnumeric(varargin{i+1}) && (numel(varargin{i+1}) == 3))
                    gridcolor = varargin{i+1};
                else
                    error('MATLAB:uitable:InvalidString', 'gridcolor has to be a valid string')
                end
            end

        case 'rowheight'
            if (isnumeric(varargin{i+1}))
                rowheight = varargin{i+1};
            else
                error('MATLAB:uitable:RowheightMustBeScalar', 'rowheight has to be a scalar')
            end

        case 'parent'
            if ishandle(varargin{i+1})
                parent = varargin{i+1};
            else
                error('MATLAB:uitable:InvalidParent', 'parent must be a valid handle')
            end

        case 'position'
            if (isnumeric(varargin{i+1}))
                position = varargin{i+1};
            else
                error('MATLAB:uitable:InvalidPosition', 'position has to be a 1x4 numeric array')
            end

        case 'columnwidth'
            if (isnumeric(varargin{i+1}))
                columnwidth = varargin{i+1};
            else
                error('MATLAB:uitable:ColumnwidthMustBeScalar', 'columnwidth has to be a scalar')
            end
        otherwise
            error('MATLAB:uitable:UnrecognizedParameter', ['Unrecognized parameter: ', varargin{i}]);
    end
end

% ---combo/check box detection--- %
if (datastatus)
    if (iscell(data))
        rownum = size(data,1);
        colnum = size(data,2);
        combo_count =0;
        check_count = 0;
        combo_box_data   = num2cell(zeros(1, colnum));
        combo_box_column = zeros(1, colnum);
        check_box_column = zeros(1, colnum);
        for j = 1:rownum
            for k = 1:colnum
                if (iscell(data{j,k}))
                    combo_box_found = true;
                    combo_count = combo_count + 1;
                    combo_box_data{combo_count} = data{j,k};
                    combo_box_column(combo_count ) = k;
                    dc = data{j,k};
                    data{j,k} = dc{1};
                else
                    if(islogical(data{j,k}))
                        check_box_found = true;
                        check_count = check_count + 1;
                        check_box_column(check_count) = k;
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the validity of the parent and/or create a figure.
if isempty(parent)
    parent = gcf; % Get the current figure. Create one if not available
end

if ( columnstatus && datastatus )
    if(size(data,2) ~= size(coln,2))
        error('MATLAB:NeedSameNumberColumns', 'Number of columns in both Data and ColumnNames should match');
    end
elseif ( ~columnstatus && datastatus )
    for i=1:size(data,2)
        coln{i} = num2str(i);
    end
    columnstatus = true;
elseif ( columnstatus && ~datastatus)
    error('MATLAB:uitable:NoDataProvided', 'No Data provided along with ColumnNames');
end

if (~exist('data_model', 'var'))
    data_model = javax.swing.table.DefaultTableModel;
end
if exist('rownum', 'var')
    data_model.setRowCount(rownum);
end
if exist('colnum', 'var')
    data_model.setColumnCount(colnum);
end

table_h= UitablePeer(data_model);

% We should have valid data and column names here.
if (datastatus), table_h.setData(data); end;
if (columnstatus), table_h.setColumnNames(coln); end;

if (combo_box_found),
    for i=1:combo_count
        table_h.setComboBoxEditor(combo_box_data(i), combo_box_column(i));
    end
end
if (check_box_found),
    for i = 1: check_count
        table_h.setCheckBoxEditor(check_box_column(i));
    end
end

% pass the specified parent and let javacomponent decide its validity.
[obj, container] = javacomponent(table_h, position, parent);
% javacomponent returns a UDD handle for the java component passed in.
table = obj;

% Have to do a drawnow here to make the properties stick. Try to restrict
% the drawnow call to only when it is absolutely required.
flushed = false;
if exist('gridcolor', 'var')
    pause(.1); drawnow;
    flushed = true;
    table_h.setGridColor(gridcolor);
end
if exist('rowheight', 'var')
    if (~flushed)
        drawnow;
    end
    table_h.setRowHeight(rowheight);
end
if exist('columnwidth', 'var')
    table_h.setColumnWidth(columnwidth);
end;

% % Add a predestroy listener so we can call cleanup on the table.
% addlistener(table, 'ObjectBeingDestroyed', {@componentDelete});

varargout{1} = table;
varargout{2} = container;

function componentDelete(src, evd)                                     %#ok
% Clean up the table here so it disengages all its internal listeners.
src.cleanup;
