function ret = spm_ov_set_ac(varargin)

% set_ac tool - plugin for spm_orthviews
%
% This tool sets the origin to the current crosshair position.
%
% This routine is a plugin to spm_orthviews for SPM2. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the matlab prompt.
%_____________________________________________________________________________
% spm_ov_set_ac.m,v 1.2 2006/11/17 Christian Gaser
% $Id: spm_ov_set_ac.m 1895 2008-07-09 08:13:57Z volkmar $

rev = '$Rev: 1895 $'; %#ok<NASGU>

global st;
if isempty(st)
    error('set_ac: This routine can only be called as a plugin for spm_orthviews!');
end;

if nargin < 2
    error('set_ac: Wrong number of arguments. Usage: spm_orthviews(''set_ac'', cmd, volhandle, varargin)');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};
switch cmd
    %-------------------------------------------------------------------------
    % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, 'Label', 'Set AC', 'Callback', ...
            ['spm_orthviews(''set_ac'',''context_init'', ', ...
            num2str(volhandle), ');'], 'Tag', ['set_ac_', num2str(volhandle)]);
        ret = item0;

    case 'context_init'

        pos = spm_orthviews('pos');
        P = st.vols{volhandle}.fname;
        fprintf('Shift origin of %s by %.1f %.1f %.1f mm\n',P,pos);
        st.vols{volhandle}.premul = spm_matrix(-pos');
        M = spm_get_space(P);
        spm_get_space(P,st.vols{volhandle}.premul*M);

    otherwise
        fprintf('spm_orthviews(''set_ac'', ...): Unknown action %s', cmd);
end;