function udat = spm_loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_loaduint8.m 1131 2008-02-06 11:17:09Z spm $

if size(V.pinfo,2)==1 && V.pinfo(1) == 2,
    mx = 255*V.pinfo(1) + V.pinfo(2);
    mn = V.pinfo(2);
else,
    spm_progress_bar('Init',V.dim(3),...
        ['Computing max/min of ' spm_str_manip(V.fname,'t')],...
        'Planes complete');
    mx = -Inf; mn =  Inf;
    for p=1:V.dim(3),
        img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
        mx  = max([max(img(:))+paccuracy(V,p) mx]);
        mn  = min([min(img(:)) mn]);
        spm_progress_bar('Set',p);
    end;
end;
spm_progress_bar('Init',V.dim(3),...
    ['Loading ' spm_str_manip(V.fname,'t')],...
    'Planes loaded');

udat = uint8(0);
udat(V.dim(1),V.dim(2),V.dim(3))=0;
rand('state',100);
for p=1:V.dim(3),
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    acc = paccuracy(V,p);
    if acc==0,
        udat(:,:,p) = uint8(round((img-mn)*(255/(mx-mn))));
    else,
        % Add random numbers before rounding to reduce aliasing artifact
        r = rand(size(img))*acc;
        udat(:,:,p) = uint8(round((img+r-mn)*(255/(mx-mn))));
    end;
    spm_progress_bar('Set',p);
end;
spm_progress_bar('Clear');
return;

function acc = paccuracy(V,p)
% if ~spm_type(V.dim(4),'intt'),
if ~spm_type(V.dt(1),'intt'),
        acc = 0;
else,
        if size(V.pinfo,2)==1,
                acc = abs(V.pinfo(1,1));
        else,
                acc = abs(V.pinfo(1,p));
        end;
end;

