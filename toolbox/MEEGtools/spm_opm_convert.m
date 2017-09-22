function Dout= spm_opm_convert(array,fnamedat,fs,scale)
% Converts array into SPM MEEG object
% FORMAT Dout = spm_opm_convert(array,fnamedat,fs,scale)
%
% array        - A numeric Array of 2,3,4 dimensions(channels,time,trials)
% fnamedat     - String specifing output path of object(include extension .dat)
% fs           - Specify sampling frequency
% scale        - Scale factor to convert to fT (defaults to 1)
% _________________________________________________________________________
%     Copyright (C) <2017>  <Tim Tierney>
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License s published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%__________________________________________________________________________

% determine output filename
[a, b] = fileparts(fnamedat);
outMat = fullfile(a,[b,'.mat']);

% if scale is not supplied set a default  of 1 
if nargin < 4
    scale = 1;
end
    
array = array.*scale;

% find number of dimensions and decide what to do based on result
dim = size(array);
L = length(dim);

if L>4
    % throw error for having unsupported number of dimensions
    ME = MException('array should not have more than 4 dimensions');
elseif L ==4
    % create MEG object with appropriate size
    Dout = meeg(dim(1),dim(2),dim(3),dim(4));
    
    % make it a blank object
    Dout= blank(Dout,fnamedat);
    
    % set the smaple rate 
    Dout = Dout.fsample(fs);
    
    % fill in data with supplied array
    Dout(1:dim(1),1:dim(2),1:dim(3),1:dim(4)) = array;
elseif L ==3
    % same with less dimensions
    Dout = meeg(dim(1),dim(2),dim(3));
    Dout= blank(Dout,fnamedat);
    Dout = Dout.fsample(fs);
    Dout(1:dim(1),1:dim(2),1:dim(3)) = array;
    
elseif L ==2 
    %same with less dimensions
    Dout = meeg(dim(1),dim(2),1);
    Dout= blank(Dout,fnamedat);
    Dout = Dout.fsample(fs);
    Dout(1:dim(1),1:dim(2),1) = array;
else
    % throw exception if  I really don't know what to do
    ME = MException('array must have between 2 and  4 dimensions');
end

% set the filename and save
Dout = fname(Dout,outMat);
Dout.save
end