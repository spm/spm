function res = fmove(this, fname)
% Method for moving or changing name of data file
% FORMAT res = fmove(this, fname)
%
% fname can be
% - path\filename -> data moved and renamed
% - path          -> data moved only
% - filename      -> data renamed only
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id: fmove.m 4484 2011-09-13 11:30:19Z guillaume $

D = struct(this);

[pth,fn,c] = fileparts(fname);
if isempty(fn)
    disp('No filename provided, using same fname.')
    [a,fn,c] = fileparts(D.fname);
end

% check data before hand, maybe not needed
[res,D]=checkmeeg(D);
orig_pfname = fullfile(D.path,D.fname);
orig_fnamedat = D.data.y.fname;

% change filename(s)
D.fname = [fn,'.mat'];
[a,b,extdat] = fileparts(D.data.fnamedat);
D.data.fnamedat = [fn,extdat];

% change directory & create it, if needed
if ~isempty(pth)
    if ~exist(pth,'dir'), mkdir(pth); end;
    D.path = pth;
end
D.data.y.fname = fullfile(D.path,D.data.fnamedat);

% move/rename data file
movefile(orig_fnamedat,D.data.y.fname);
delete(orig_pfname)
save(fullfile(D.path,D.fname),'D','-V6')
