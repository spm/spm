function ret = spm_ov_goto_max(varargin)
% Goto maximum intensity tool - plugin for spm_orthviews
%
% This tool provides capabilities similar to the "Goto ... maximum"
% functionality in spm_mip_ui.m. When the tool is called for the first
% time, it has to read the whole image data file. This might result in a
% slow response depending on the image dimensions.
%
% This routine is a plugin to spm_orthviews. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the MATLAB prompt.
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% $Id: spm_ov_goto_max.m 4996 2012-10-11 18:28:37Z guillaume $

global st;
if isempty(st)
    error('goto_max: This routine can only be called as a plugin for spm_orthviews!');
end

if nargin < 2
    error('goto_max: Wrong number of arguments. Usage: spm_orthviews(''goto_max'', cmd, volhandle, varargin)');
end

cmd = lower(varargin{1});
volhandle = varargin{2};

switch cmd

    % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, 'Label', 'Goto maximum');
        item1 = uimenu(item0, 'Label', 'Goto global maximum',...
            'Callback', ['spm_orthviews(''goto_max'', ''global'',' num2str(volhandle) ');']);
        item2 = uimenu(item0, 'Label', 'Goto nearest local maximum',...
            'Callback', ['spm_orthviews(''goto_max'', ''local'',' num2str(volhandle) ');']);
        ret = item0;

    case 'global'
        if ~isfield(st.vols{volhandle}, 'goto_max')
            [dat, xyz] = spm_read_vols(st.vols{volhandle});
            [unused, mxind] = max(dat(:));
            st.vols{volhandle}.goto_max.globalmm = xyz(:, mxind);
        end
        posmm = st.vols{volhandle}.premul*[st.vols{volhandle}.goto_max.globalmm;1];
        spm_orthviews('reposition', posmm(1:3));
        
    case 'local'
        % Read volume, goto local maxima higher or equal current voxel
        % intensity
        dat       = spm_read_vols(st.vols{volhandle});
        [x, y, z] = ndgrid(1:st.vols{volhandle}.dim(1), 1:st.vols{volhandle}.dim(2), 1:st.vols{volhandle}.dim(3));
        xyz       = [x(:) y(:) z(:)]';
        posvx     = round((st.vols{volhandle}.premul*st.vols{volhandle}.mat)\[spm_orthviews('pos');1]);
        try
            if isfinite(dat(posvx(1), posvx(2), posvx(3)))
                sel = isfinite(dat(:)) & dat(:) >= dat(posvx(1), posvx(2), posvx(3)) - eps;
            else
                sel = isfinite(dat(:));
            end
        catch
            sel = isfinite(dat(:));
        end
        [unused, unused, XYZ]  = spm_max(dat(sel), xyz(:,sel));
        XYZdist      = bsxfun(@minus,XYZ,posvx(1:3));
        [unused, nmaxind] = min(sum(XYZdist.^2));
        posmm        = st.vols{volhandle}.premul*st.vols{volhandle}.mat*[XYZ(:, nmaxind); 1];
        spm_orthviews('reposition', posmm(1:3));
end
