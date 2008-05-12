function [table, container] = uitable(varargin)
% WARNING: This feature is not supported in MATLAB
% and the API and functionality may change in a future release.
% (This code has been hacked from a MATLAB trial function.)

error(javachk('awt'));
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

function componentDelete(src, evd)                                     %#ok
% Clean up the table here so it disengages all its internal listeners.
src.cleanup;
