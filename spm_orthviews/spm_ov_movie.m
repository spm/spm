function ret = spm_ov_movie(varargin)

% Movie tool - plugin for spm_orthviews
%
% This plugin allows an automatic "fly-through" through all displayed
% volumes. Apart from pre-defined trajectories along the x-, y- and z-axis,
% resp., it is possible to define custom start and end points (in mm) for
% oblique trajectories.
%
% @(#) spm_ov_movie.m,v 1.3 2003/06/04 14:37:49 glauche Exp

global st;
if isempty(st)
  error('movie: This routine can only be called as a plugin for spm_orthviews!');
end;

if nargin < 2
  error('movie: Wrong number of arguments. Usage: spm_orthviews(''roi'', cmd, volhandle, varargin)');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};

switch cmd
  
  %-------------------------------------------------------------------------
  % Context menu and callbacks
  case 'context_menu'  
    item0 = uimenu(varargin{3}, 'Label', 'Movie tool', 'Callback', ...
	['feval(''spm_ov_movie'',''context_init'', ', ...
	  num2str(volhandle), ');'], 'Tag', ['MOVIE_0_', num2str(volhandle)]);
    
  case 'context_init'
    Finter = spm_figure('FindWin', 'Interactive');
    spm_input('!DeleteInputObj',Finter);
    dir=spm_input('Select movie direction', '!+1', 'b', 'x|y|z|custom', ...
	[1 2 3 0], 1);
    if dir==0
      mstart=spm_input('First point (mm)', '!+1', 'e', '', [1 3]);
      mend  =spm_input('Final point (mm)', '!+1', 'e', '', [1 3]);
    else
      mstart=[0 0 0];
      mend=[0 0 0];
      dirs='XYZ';
      tmp=spm_input([dirs(dir) ' intervall (mm)'], '!+1', 'e', ...
	  num2str(st.bb(:,dir)', '%.1f %.1f'), 2);
      mstart(dir)=tmp(1);
      mend(dir)=tmp(2);
    end;
    ds=spm_input('Step size (mm)', '!+1', 'e', '1', 1);
    opos=spm_orthviews('pos');
    d=mend-mstart;
    l=sqrt(d*d');
    d=d./l;
    for k=0:ds:l
      spm_orthviews('reposition', mstart+k*d);
    end;
    spm_orthviews('reposition', opos);
    spm_input('!DeleteInputObj',Finter);
  otherwise    
    fprintf('spm_orthviews(''movie'', ...): Unknown action %s', cmd);
end;


spm('pointer','arrow');
