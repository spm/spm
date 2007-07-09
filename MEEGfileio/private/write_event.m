function write_event(filename, event, varargin);

% WRITE_EVENT writes an event structure to a file, a message daemon
% listening on a network socked, or to another computer connected through
% the serial port.
%
% Use as
%   write_event(filename, event, ...)
%
% The first argument is a string containing the filename. The second
% argument is a structure with the event. Multiple events can be
% represented as a structure array.
%
% Events are represented as
%   event.type     = string
%   event.sample   = expressed in samples, first sample of recording is 1
%   event.value    = number or string
%   event.offset   = expressed in samples
%   event.duration = expressed in samples
%
% Events can also be written to special communication streams
% by specifying the target as URI instead of a filename. Supported are
%   fifo://<filename>
%   tcpsocket://<host>:<port>
%   mysql://<user>:<password>@<host>:<port>
%   serial:<port>?key1=value1&key2=value2&...
%
% See also READ_HEADER, READ_DATA, READ_EVENT, WRITE_DATA


% Copyright (C) 2007, Robert Oostenveld
%
% $Log: write_event.m,v $
% Revision 1.15  2007/06/19 11:11:28  chrhes
% changed the implementation of how the target serial port is located if it
% exists
%
% Revision 1.14  2007/06/19 10:11:37  chrhes
% restricted the functionality of serial port writing such that only the
% content of the field event.value is written as a character array
%
% Revision 1.13  2007/06/13 14:46:35  roboos
% removed type/subtype, added type/value/sample for mysql
%
% Revision 1.12  2007/06/13 08:06:45  roboos
% updated help
%
% Revision 1.11  2007/06/12 19:37:28  roboos
% added support for mysql port specification (default is 3306)
% added support for writing multiple events to mysql, each requires a single query
%
% Revision 1.10  2007/06/12 16:34:49  roboos
% first implementation of writing to a mysql database: approx. 30ms per event
%
% Revision 1.9  2007/06/07 12:43:38  chrhes
% added an option for specifying the maximum length of the event queue for
% the case where the events are being written (appended) to a .mat file
%
% Revision 1.8  2007/06/06 20:14:49  chrhes
% fixed a small bug to do with sting comparison
%
% Revision 1.7  2007/06/06 16:00:31  chrhes
% updated some documentation
%
% Revision 1.6  2007/06/06 15:55:30  chrhes
% extended functionality for serial port writing so that the code recognises
% when to write the whole event structure or only a single control character
%
% Revision 1.5  2007/06/06 15:45:53  chrhes
% implemented option for writing events to the serial port
%
% Revision 1.4  2007/06/06 12:38:57  roboos
% write events to uncompressed matlab v6 file
%
% Revision 1.3  2007/06/06 07:14:29  roboos
% switched to using filetype_check_uri for detection and parsing of filename
% switched to using external struct2msg function
% implemented fcdc_fifo
% added optinoal input arguments (key-value pairs)
%
% Revision 1.2  2007/05/31 09:54:05  roboos
% implemented writing events to a plain matlab file, the default is to append them
%
% Revision 1.1  2007/05/31 09:14:34  roboos
% initial implementation, sofar only tcpsocket
%

% set the defaults
eventformat = keyval('eventformat', varargin); if isempty(eventformat), eventformat = filetype(filename); end
swapping    = keyval('swapping',    varargin); if isempty(swapping),    swapping = 'native';              end
append      = keyval('append',      varargin); if isempty(append),      append = 'yes';                   end
maxqlength  = keyval('maxqlength',  varargin); if isempty(maxqlength),  maxqlength = Inf;                 end

switch eventformat
  case 'disp'
    % display it on screen, this is only for debugging
    disp(event);

  case 'fcdc_serial'
    % serial port on windows or linux platform
    s = [];
    [port, opt] = filetype_check_uri(filename);
    % determine whether any serial port objects are already associated with the
    % target serial port
    temp = instrfind;
    if isa(temp,'instrument')
      % find all serial ports
      i1 = strcmp('serial',{temp(:).Type});
      if any(i1)
        % find all serial ports whose name matches that of the specified port
        i2 = strmatch(lower(['Serial-',port]),lower({temp(find(i1)==1).Name}));
        % set s to the (first) matching port if present (and open if necessary)
        if ~isempty(i2)
          s = temp(i2(1));
          if ~strcmp(s.Status,'open'), fopen(s); end;
        end
      end
    end
    % create, configure a serial port object if necessary and open the port
    if ~isa(s,'serial')
      s = serial(port);
      if ~isempty(opt) && iscell(opt), s = set(s,opt); end;
      fopen(s);
    end

    %     % convert the event structure into an appropriate message
    %     if isfield(event,'type') && strcmp(event.type,'ctrlchar')
    %       % use only a single control character
    %       msg = char(event.value(1));
    %     else
    %       % convert the entire event structure into a message
    %       msg = struct2msg(event);
    %     end

    % write the contents of the field event.value to the serial port as a string
    if isfield(event,'value') && ~isempty(event.value)
      msg = char(event.value);
    else
      msg = [];
    end
    % write the message to the serial port
    if ~isempty(msg) && isa(s,'serial') && strcmp(s.Status,'open')
      fprintf(s,msg,'async');
    else
      error('could not write event to serial port');
    end
 
  case 'fcdc_mysql'
    persistent prev_filename
    % write to a MySQL server listening somewhere else on the network
    % there should be a database named 'fieldtrip' containing the table named 'event'
    % FIXME the port is currently not yet used
    if ~strcmp(filename, prev_filename)
      % close the database
      mysql('close');
      prev_filename = [];
    end
    if ~strcmp(filename, prev_filename)
      % open the database
      [user, password, server, port] = filetype_check_uri(filename);
      if ~isempty(port)
        server = sprintf('%s:%d', server, port);
      end
      try
        mysql('open', server, user, password);
        prev_filename = filename;
      catch
        prev_filename = [];
      end
    end
    for i=1:length(event)
      % make a structure with the same elements as the fields in the database table
      s = struct;
      % the type, value and sample field also exist as elements in the
      % table and as such can be used for filtering
      if isa(event(i).type, 'char')
        s.type = event(i).type;
      end
      if isa(event(i).value, 'numeric') && max(size(event(i).value))==1
        s.value = event(i).value;
      end
      if isa(event(i).sample, 'numeric') && max(size(event(i).sample))==1
        s.sample = event(i).sample;
      end
      % convert the event into a network message
      s.data    = struct2msg(event(i));
      s.length  = length(s.data);
      cmd = insert_query('fieldtrip.event', s);
      % insert the structure into the database table
      mysql(cmd);
    end

  case 'fcdc_fifo'
    % these are opened in blocking mode, i.e. reading/writing will block until boths sides are connected
    fifo = filetype_check_uri(filename);
    if ~exist(fifo)
      error(sprintf('the FIFO %s does not exist', fifo));
    end
    fid = fopen(fifo, 'w');
    for i=1:length(event)
      % convert the event into a network message
      msg = struct2msg(event(i));
      num = fwrite(fid, msg);
      if num~=length(msg)
        error(sprintf('problem writing to FIFO %s', fifo));
      end
    end
    fclose(fid);

  case 'fcdc_tcpsocket'
    % TCP network socket
    [host, port] = filetype_check_uri(filename);
    con = pnet('tcpconnect', host, port);
    if con<0
      error('problem opening network connection');
    end
    for i=1:length(event)
      % convert the event into a network message
      msg = struct2msg(event(i));
      % tell the message daemon that a message will be sent, and send it
      pnet(con,'write', 'c', swapping);
      pnet(con,'write', uint8(0), swapping);                    % type
      pnet(con,'write', uint8(0), swapping);                    % subtype
      pnet(con,'write', uint16(length(msg)), swapping);         % size
      pnet(con,'write', char(msg), swapping);                   % message data
    end
    pnet(con,'close')

  otherwise
    % assume that it is a file. Since the file probably does not yet
    % exist, determine its type by only looking at the extension
    if filetype_check_extension(filename, '.mat')
      % write the events to a matlab file
      if exist(filename) && strcmp(append, 'yes')
        try
          tmp = load(filename, 'event');
          event = cat(1, tmp.event(:), event(:));
        catch
          event = event(:);
        end
        % optionally restric the length of the event queue to flush old events
        if isfinite(maxqlength) && isreal(maxqlength) && (maxqlength>0) && (length(event)>maxqlength)
          event = event(end-maxqlength+1:end);
          % NOTE: this could be done using the filter event function, but
          % then this is just a temporary solution that will probably be
          % removed in a future versions of the code
        end
        save(filename, 'event', '-append', '-v6');
        % NOTE: the -append option in this call to the save function does
        % not actually do anything useful w.r.t. the event variable since the
        % events are being appended in the code above and the the save function
        % will just overwrite the existing event variable in the file.
        % However, if there are other variables in the file (whatever) then the
        % append option preservs them
      else
        save(filename, 'event', '-v6');
      end
    else
      error('unsupported file type')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction, this is copied from testing/db_insert.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = insert_query(tablename, s);
name  = fieldnames(s);
value = {};
numel = length(name);
for i=1:numel
  value{i} = getfield(s, name{i});
end
% these are the table rows
str = ['insert into ' tablename ' ( '];
for i=1:numel
  str = [str name{i} ', '];
end
str = str(1:(end-2));  % remove the last ', '
str = [str ' ) values ( '];
% these are the corresponding values
for i=1:numel
  if isempty(value{i})
    str = [str '""' ];
  elseif isnumeric(value{i})
    str = [str num2str(value{i})];
  elseif isstr(value{i})
    str = [str '"' value{i} '"' ];
  else
    error('unsuported data type in structure');
  end
  str = [str ', '];
end
str = str(1:(end-2));  % remove the last ', '
str = [str ' ) ;'];

