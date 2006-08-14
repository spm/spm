function ret = spm_ov_rgb(varargin)
% RGB overlays
% A shorthand to overlaying the absolute value of three different images
% onto a displayed image in colours red, green and blue. The overlay images
% are optionally masked and multiplied with a scaling image. The displayed
% overlay images are the absolute value of the given overlays.
%
% This routine is a plugin to spm_orthviews for SPM5. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the matlab prompt.
%_______________________________________________________________________
%
% @(#) $Id: spm_ov_rgb.m 591 2006-08-14 11:06:49Z volkmar $

rev = '$Revision: 591 $';

global st;
if isempty(st)
  error('rgb: This routine can only be called as a plugin for spm_orthviews!');
end;

if nargin < 2
  error('rgb: Wrong number of arguments. Usage: spm_orthviews(''rgb'', cmd, volhandle, varargin)');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};

switch cmd
  
  %-------------------------------------------------------------------------
  % Context menu and callbacks
  case 'context_menu'  
    item0 = uimenu(varargin{3}, 'Label', 'RGB overlays');
    item1 = uimenu(item0, 'Label', 'Add', 'Callback', ...
	['feval(''spm_ov_rgb'',''context_init'', ', ...
	  num2str(volhandle), ');'], 'Tag', ['RGB_0_', num2str(volhandle)]);
    item1 = uimenu(item0, 'Label', 'Help', 'Callback', ...
	  ['feval(''spm_help'',''' mfilename ''');']);
    
  case 'context_init'
    Finter = spm_figure('FindWin', 'Interactive');
    spm_input('!DeleteInputObj',Finter);
    [Vqfnames sts] = spm_select(3, 'image',...
                                'Components of 1st eigenvector', ...
                                [], pwd, 'evec1.*');
    if ~sts return; end;
    Vq = spm_vol(Vqfnames);
    [Vfafname sts] = spm_select([0 1],'image','FA image (optional)', [], ...
                                pwd, 'fa.*'); 
    Vfa = spm_vol(Vfafname);
    [Vmaskfname sts] = spm_select([0 1],'image','Mask image (optional)');
    Vmask = spm_vol(Vmaskfname);
    spm('pointer','watch');
    Vamq = rmfield(Vq,'private');
    for k=1:3
      [p n e v]=fileparts(Vq(k).fname);
      sel = 2*isempty(Vmask)+isempty(Vfa);
      switch(sel)
	case 0, %both Vmask and Vfa set
	  Vamq(k).fname=fullfile(p,['abs_msk_fa_' n e v]);
	  spm_imcalc([Vq(k) Vfa Vmask],Vamq(k),'abs(i1.*i2.*i3)',{[],1,[]});
	case 1, %only Vmask set
	  Vamq(k).fname=fullfile(p,['abs_msk_' n e v]);
	  spm_imcalc([Vq(k) Vmask],Vamq(k),'abs(i1.*i2)',{[],1,[]});
	case 2, %only Vfa set
	  Vamq(k).fname=fullfile(p,['abs_fa_' n e v]);
	  spm_imcalc([Vq(k) Vfa],Vamq(k),'abs(i1.*i2)',{[],1,[]});
	case 3, %nothing set
	  Vamq(k).fname=fullfile(p,['abs_' n e v]);
	  spm_imcalc(Vq(k),Vamq(k),'abs(i1)',{[],1,[]});
      end;
    end;
    spm_orthviews('addcolouredimage',volhandle,Vamq(1).fname,[1 0 0]);
    spm_orthviews('addcolouredimage',volhandle,Vamq(2).fname,[0 1 0]);
    spm_orthviews('addcolouredimage',volhandle,Vamq(3).fname,[0 0 1]);
    spm_orthviews('redraw');
    
    spm_input('!DeleteInputObj',Finter);
 case 'redraw'
  % Do nothing
  otherwise    
    fprintf('spm_orthviews(''rgb'', ...): Unknown action %s', cmd);
end;

spm('pointer','arrow');
