% read_eeglabdata() - import EEGLAB dataset files
%
% Usage:
%   >> dat = read_eeglabdata(filename);
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'begsample' - [integer] first sample to read
%   'endsample' - [integer] last sample to read
%   'begtrial'  - [integer] first trial to read, mutually exclusive with begsample+endsample
%   'endtrial'  - [integer] last trial to read, mutually exclusive with begsample+endsample
%   'chanindx'  - [integer] list with channel indices to read
%   'header'    - FILEIO structure header
%
% Outputs:
%   dat    - data over the specified range
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2008-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: read_eeglabdata.m,v $
% Revision 1.1  2008/04/18 14:04:48  roboos
% new implementation by Arno, shoudl be tested
%

function dat = read_eeglabdata(filename, varargin);

if nargin < 1
  help read_eeglabdata;
  return;
end;

header    = keyval('header',     varargin);
begsample = keyval('begsample',  varargin);
endsample = keyval('endsample',  varargin);
begtrial  = keyval('begtrial',   varargin);
endtrial  = keyval('endtrial',   varargin);
chanindx  = keyval('chanindx',   varargin);

if isempty(header)
  header = read_eeglabheader(filename);
end;

if isempty(begsample), begsample = 1; end;
if isempty(endsample), endsample = header.nSamples; end;
if isempty(begtrial), begtrial = 1; end;
if isempty(endtrial), endtrial = header.nTrials; end;
if isempty(chanindx), chanindx = [1:header.nChans]; end;

if ischar(header.orig.data)
  if strcmpi(header.orig.data(end-2:end), 'set'),
    header.ori = load('-mat', filename);
  else
    fid = fopen(header.orig.data);
    if fid == -1, error('Cannot not find data file'); end;
    dat = fread(fid,prod([ header.nChans header.nSamples*header.nTrials ]),'float');
    fclose(fid);
    dat = reshape(dat, header.nChans, header.nSamples, header.nTrials);
  end;
else
  dat = header.orig.data;
end;

dat = dat(chanindx, begsample:endsample, begtrial:endtrial);
