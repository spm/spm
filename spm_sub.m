function spm_sub(s)
% displays the son subroutines called by s
% FORMAT spm_sub(s)
% s   - filename {e.g. spm_routine.m}
%____________________________________________________________________________
% 
% spm_sub first invokes spm_help_disp and then finds all occurences
% of spm_* (the spm routines called or referenced) in the specified
% file.  These strings are displayed in Figure(2) with appropriate
% callback strings, that again call spm_sub.  This recursive arrangement
% allows successive access to lower and lower levels of subroutines.
%
%__________________________________________________________________________
% %W% %E%

% display the current help comments
%----------------------------------------------------------------------------
spm_help_disp(s); figure(2); cla; axis off

% find subroutines and create text objects and their CallBacks
%----------------------------------------------------------------------------
fid   = fopen(s,'r');
S     = setstr(fread(fid))'; fclose(fid);
q     = findstr(S,'spm_');
H     = 0.7;
y     = H;
x     = 0;
U     = s;

for i = 1:length(q)
    d = [0:32] + q(i);
    Q = S(d(d <= length(S)));
    d = find((Q == ';') | (Q == '(') | (Q == 10) | (Q == '.') | (Q == ' '));
	if length(d)
	  Q   = [Q(1:(min(d) - 1)) '.m'];
	  if exist(Q) == 2 & ~length(findstr(U,Q))
	    c = ['spm_sub(''' Q ''');'];
	    text(x,y,Q,'ButtonDownFcn',c,'FontWeight','Bold','Color',[1 0 0])
	    y = y - 0.06;
	    U = [U ';' Q];
	  end
    end
    if y < -0.1; x = x + 0.5; y = H; end
end

% heading
%----------------------------------------------------------------------------
if y < H | x > 0
	text(0.32,0.92,'SPM routines','FontWeight','Bold','FontSize',16);
end
