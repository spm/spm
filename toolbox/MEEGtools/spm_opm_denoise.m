function Dout= spm_opm_denoise(D,refD,derivative,gs,update,prefix)
% Denoises OPM data with specified regressors, optionally with derivatives
% and global signal
% FORMAT Dout - spm_opm_denoise(D,refD,derivative,gs,update,prefix)
%
% D            - MEEG object
% refD         - array or MEEG object containing regresssors for denoising.
% derivative   - boolean indicating whether to use derivatives or not.
% gs           - boolean indicating whether to use Global Signal or not.
% update       - boolean indicating whether to create MEEG object.
% prefix       - string to prefix file path with if update is TRUE.
% _________________________________________________________________________
% Copyright (C) <2017>  <Tim Tierney>
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
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



% check to make sure dimensions match between reference and MEEG object
ndim = size(D);
ldim = size(refD);
if (~all(ndim(2:3)==ldim(2:3)))
    ME = MException('Dimensions of D and refD should be equal');
    throw(ME)
end


% ref, sensors and residual objects 
ref = refD;
sensors = D;
res = zeros(size(sensors));

% need to loop over Ntrials of size winSize
Ntrials = size(sensors,3);
winSize = size(sensors,2);

% loop (inefficient due to continued memory reallocation but ...)
for i =1:Ntrials
     % add a mean column to the reference regressors
     intercept = ones(winSize,1);
     reference = ref(:,:,i)';
     reference=[reference intercept];
    
     % optionally add derivatives
     if(derivative)
        drefer=diff(reference);
        drefer=[drefer(1,:); drefer];
        reference=[drefer reference];
     end
     
    % optionally add global signal 
    if(gs)
        trial = D(:,:,i)';
        gsrefer =  mean(trial,2);
        reference=[gsrefer reference];
    end
    % reference is column major so transpose sensors
    beta = pinv(reference)*sensors(:,:,i)';
    res(:,:,i) = (sensors(:,:,i)'- reference*beta)';
end

% make sure output has the singleton dimension if matrix supplied
if ((length(size(res)))==2)
    outsize = [size(res) 1];
else
    outsize = size(res);
end


if(update)
    % name, clone and fill with residual data 
    inFile = fnamedat(D);
    [a b]=fileparts(inFile);
    outfile = fullfile(a,[prefix b '.dat']);
    Dout = clone(D,outfile,outsize);
    Dout(:,:,:) = res;
else
    Dout = res;
end

end