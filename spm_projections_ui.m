function [R1,R2,R3,R4,R5,R6]=spm_projections_ui(Action,P2)
% used to review results of statistical analysis (SPM{Z})
% FORMAT spm_projections_ui
%_______________________________________________________________________
%
% spm_projections_ui allows the SPM{Z} created by spm_spm.m to be re-displayed
% and characterized in terms of regionally significant effects.  Note that
% the threshold does not have to be the same as in the original analysis
%
% Multiple [orthogonal] contrasts can be specified to produce a SPM{Z} that
% reflects the significance of two or more effects:-
%
% specifying a vector for the contrasts causes the second (and ensuing)
% contrasts to mask the first.  Non-orthogonal compounds are allowed.
% The ensuing voxels reach criteria (at the uncorrected threshold
% specified) for all the contrasts.  The statistic constituting the
% SPM{Z} is the mean of the n component Z values divided by sqrt(n).
% Given the orthogonality constraint above, this statistic is itself
% Gaussian.  This use of spm_projections_ui.m is useful for testing
% multiple hypothese simultaneously.  For example the conjunction
% of activations in several task-pairs or (in multifactorial designs)
% the analysis of interaction effects in (and only in) areas subject to a
% main effect.  The resulting p values are generally so small that one
% can forgo a correction for multiple comparisons.
%
% see spm_projections.m for further details
%
%_______________________________________________________________________
% %W% Karl Friston %E%

%-Format arguments
if nargin==0, Action='Display'; end

if strcmp(lower(Action),lower('Display'))
%=======================================================================

%-Get SPMt.mat for analysis, and load results
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
spm_clf(Finter)
set(Finter,'Name','SPM{Z} projections')

tmp = spm_get(1,'.mat','select SPMt.mat for analysis','SPMt');
global CWD
CWD = strrep(tmp,'/SPMt.mat','');
K   = [];

load([CWD,'/SPM'])
load([CWD,'/XYZ'])
load([CWD,'/SPMt'])

%-Get contrast[s]
%-----------------------------------------------------------------------
i    = 0;
while any(i < 1 | i > size(CONTRAST,1))
	i = spm_input(sprintf('contrast[s] ? 1 - %i',size(CONTRAST,1)),1);
	c = CONTRAST(i,:);
end

%-Get height threshold [default = 3.2]
%-----------------------------------------------------------------------
u    = spm_input('height threshold {Z value}',2,'e',3.2);

%-Get extent threshold [default = E{n} - expected voxels per cluster]
% Omit spatial extent threshold for multiple contrasts.
%-----------------------------------------------------------------------
if size(c,1) == 1
	[P,EN,Em,En,Pk] = spm_P(1,W,u,0,S);
	k    = spm_input('extent threshold {voxels}',3,'e',round(En));
else
	k    = 0;
end


%-Pass arguments to spm_projections
%-----------------------------------------------------------------------
set(Finter,'Name','Thankyou','Pointer','watch')

[t,XYZ] = spm_projections(SPMt(i,:),XYZ,u,k,V,W,S,[K H C B G],c,df);

%-If no suprathreshold voxels then return
%-----------------------------------------------------------------------
if isempty(t), spm_clf(Finter), return, end

%-Store essential variables in UserData of various text objects,
% Create button to write filtered SPM{Z}
%=======================================================================

set(Finter,'Name','SPM{Z} projections','Pointer','Arrow')

h = findobj(Fgraph,'Tag','Empty');
set(h(1),'Tag','XYZ','UserData',XYZ)
set(h(2),'Tag','t','UserData',t)
set(h(3),'Tag','V','UserData',V)
set(h(4),'Tag','u','UserData',u)
set(h(5),'Tag','k','UserData',k)
if exist('FLIP')~=1, FLIP=0; end
set(h(6),'Tag','FLIP','UserData',FLIP)

%-Determine positioning for Buttons (as in spm_input)
%-------------------------------------------------------------------
set(Finter,'Units','pixels')
FigPos = get(Finter,'Position');
Xdim = FigPos(3); Ydim = FigPos(4);
a = 5.5/10;
YPos = 5;
y = Ydim - 30*YPos;
PPos = [10,     y, (a*Xdim -20),     20];
MPos = [PPos(1), PPos(2), Xdim-50, 20];
RPos = MPos + [0,0,30,0];

%-Delete any previous inputs using position YPos
%-------------------------------------------------------------------
delete(findobj(Finter,'Tag',['GUIinput_',int2str(YPos)]))


uicontrol(Finter,'Style','PushButton',...
	'String','Write filtered SPM{Z} to Analyze file',...
	'Tag','WriteButton',...
	'Callback',['set(gco,''Visible'',''off''),',...
		'spm_projections_ui(''WriteFiltered'')'],...
	'Interruptible','yes',...
	'Tag',['GUIinput_',int2str(YPos)],...
	'Position',RPos)
%uicontrol(Finter,'Style','PushButton',...
%	'String','Get details of filtered image',...
%	'Callback',[...
%		'[XYZ,t,V,u,k,FLIP]=',...
%			'spm_projections_ui(''GetFilteredPrams'');',...
%		'disp(''XYZ,t,V,u,k,FLIP now exist in base workspace'')'],...
%	'Interruptible','yes',...
%	'Tag',['GUIinput_',int2str(YPos+1)],...
%	'Position',RPos-[0 30 0 0])

%-Finished
%-----------------------------------------------------------------------
return


elseif strcmp(lower(Action),lower('WriteFiltered'))
%=======================================================================
% spm_projections_ui('WriteFiltered',FName)
[XYZ,t,V,u,k,FLIP] = spm_projections_ui('GetFilteredPrams');
if isempty(XYZ), disp('No Voxels'), return, end

%-Get filename
%-----------------------------------------------------------------------
if nargin<2
	FName=spm_input('Filename ?',5,'s','SPM_filtered');
else
	FName = P2;
end
str=sprintf('spm{Z}-filtered: u= %5.3f, k=%d',u,k);

%-Reconstruct filtered image from XYZ & t
%-----------------------------------------------------------------------
n       = size(XYZ,2);
%-Unflip flipped images - Negate X values if FLIP==1
rcp     = round(XYZ./meshgrid([1-2*FLIP;1;1].*V(4:6),1:n)' + ...
	meshgrid(V(7:9),1:n)');
DimMult = cumprod([1,V(1:2)']);
OffSets = meshgrid([0,1,1],1:n)';
e       = ((rcp-OffSets)'*DimMult')';
t       = t.*(t>0); %-Ensure positivity of z-values
T       = zeros(1,prod(V(1:3)));
T(e)    = t;

%-Write out to analyze file
%-----------------------------------------------------------------------
spm_hwrite([FName,'.hdr'],V(1:3),V(4:6),1/16,2,0,V(7:9),str);
fid = fopen([FName,'.img'],'w');
fwrite(fid,T*16,spm_type(2));
fclose(fid);

spm_clf(Finter)

return

elseif strcmp(lower(Action),lower('GetFilteredPrams'))
%=======================================================================
% [XYZ,t,V,u,k,FLIP] = spm_projections_ui('GetFilteredPrams')
Fgraph = spm_figure('FindWin','Graphics');
R1 = get(findobj(Fgraph,'Tag','XYZ'),'UserData');
R2 = get(findobj(Fgraph,'Tag','t'),'UserData');
R3 = get(findobj(Fgraph,'Tag','V'),'UserData');
R4 = get(findobj(Fgraph,'Tag','u'),'UserData');
R5 = get(findobj(Fgraph,'Tag','k'),'UserData');
R6 = get(findobj(Fgraph,'Tag','FLIP'),'UserData');
return

else
%=======================================================================
error('Unknown action string')

%=======================================================================
end
