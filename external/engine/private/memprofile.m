function varargout = memprofile(varargin)

% MEMPROFILE ON starts the profiler and clears previously recorded
% profile statistics.
%
% MEMPROFILE OFF stops the profiler.
%
% MEMPROFILE RESUME restarts the profiler without clearing
% previously recorded memory statistics.
%
% MEMPROFILE CLEAR clears all recorded profile statistics.
%
% MEMPROFILE REPORT displays a summary of the recorded statistics.
%
% STATS = MEMPROFILE('INFO') returns a structure containing the
% current profiler statistics.
%
% Example use:
%   memprofile on
%   x = {};
%   for i=1:100
%     x{i} = zeros(1000,1000); % 8kB per item
%     disp(i);
%     pause(0.1);
%   end
%   memprofile report
%
%   stat = memprofile('info');
%   plot([stat.time], [stat.mem])
%   xlabel('time(s)'); ylabel('memory (bytes)')
%
% See also PROFILE

% -----------------------------------------------------------------------
% Copyright (C) 2011, Robert Oostenveld
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/
% -----------------------------------------------------------------------

varargout = {struct('mem',NaN,'time',NaN,'bytes',NaN)};
