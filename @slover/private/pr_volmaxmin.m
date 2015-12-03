function [mx,mn] = pr_volmaxmin(vol)
% Return max and min value in image volume
% FORMAT [mx,mn] = pr_volmaxmin(vol)
%
% Input
% vol      - image name or vol struct
%
% Outputs
% mx       - maximum
% mn       - minimum
%__________________________________________________________________________

% Matthew Brett
% $Id: pr_volmaxmin.m 6623 2015-12-03 18:38:08Z guillaume $

if nargin < 1
    error('Need volume to process');
end
if ischar(vol)
    vol = spm_vol(vol);
end
if ~isstruct(vol)
    error('vol did not result in vol struct');
end
if mars_struct('isthere', vol, 'imgdata')
    tmp = vol.imgdata(isfinite(vol.imgdata));
    mx = max(tmp);
    mn = min(tmp);
else
    mx = -Inf;mn=Inf;
    for i=1:vol.dim(3),
        tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),[0 NaN]);
        tmp = tmp(find(isfinite(tmp(:))));
        if ~isempty(tmp)
            mx = max([mx; tmp]);
            mn = min([mn; tmp]);
        end
    end
end
