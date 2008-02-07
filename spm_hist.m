function h = spm_hist(ind,val)
% Generate a weighted histogram
% FORMAT h = spm_hist(ind,wt)
% ind - indices (unsigned byte)
% val - weights
%
% spm_hist.c should be compiled for optimal speed.  If it is not
% compiled, then the following code is run:
%     h = full(sparse(double(ind)+1,ones(size(ind)),wt,256,1));
% For Matlab 7, the accumarray function could be used 
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_hist.m 1143 2008-02-07 19:33:33Z spm $


persistent flg

%if isempty(flg),
%    flg = 0;
%    try,
%        [pth,nam,ext] = fileparts(which('spm_hist.c'));
%        mex(fullfile(pth,[nam,ext]),'-O','-outdir',pth);
%    catch,
%    end;
%end;

if isempty(flg),
    try
        h   = accumarray(double(ind(:))+1,double(val(:)),[256 1]);
    catch
        flg = true;
        h   = full(sparse(double(ind)+1,ones(size(ind)),val,256,1));
    end;
else
    h = full(sparse(double(ind)+1,ones(size(ind)),val,256,1));
end;
