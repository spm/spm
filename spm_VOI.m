% Tabular display of adjusted data
% FORMAT spm_VOI
%_______________________________________________________________________
%
% spm_VOI is called by spm_results and takes variables in working memory to
% compute a p value corrected for a specified volume of interest for,   and
% centred on, the current voxel.
%_______________________________________________________________________
% %E% Karl Friston %W%


% Find nearest voxel [in a Euclidean sense] in point list & update GUI
%-----------------------------------------------------------------------
[L,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),XYZ);
spm_XYZreg('SetCoords',L,hReg);
A     = spm_clusters(XYZ,V([4 5 6]));
ML    = XYZ(:,find(A == A(i)));
Z     = t(i);
N     = length(ML);


% Specify search volume
%--------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
set(Finter,'Name','small volume correction');
SPACE  = spm_input('search volume',1,'b','Sphere|Box|Cluster',['S' 'B' 'V']);
if     SPACE == 'S'
	r    = spm_input('radius of spherical VOI {mm}',2);
	str  = sprintf('%0.1f mm sphere',r);
	r    = r/exp(mean(log(V(4:6))));

elseif SPACE == 'B'
	r    = spm_input('box dimensions [k l m] {mm}',2);
	str  = sprintf('%0.1f x %0.1f x %0.1f mm box',r(1),r(2),r(3));
	r    = r(:)./V(4:6);

elseif SPACE == 'V'
	r    = ML;
	str  = sprintf('%0.0f voxel cluster',N);

end

% p values
%-----------------------------------------------------------------------
[Pu EC] = spm_Pec(SPMdist,SPACE,Z,r,W);		% [un]corrected p value (u)
Pn      = 1 - spm_kcdf(N,U,W);			%   uncorrected p value (k)



% Display results
%=======================================================================

%-Delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
spm_results_ui('ClearPane',Fgraph)

% Table headings
%-----------------------------------------------------------------------
axes('Position',[0.1 0.06 0.8 0.46]); axis off
lab     = ['Small volume correction:  ',spm_str_manip(CWD,'a50')];
text(0,23,['Search volume:  ' str],'FontSize',12,'FontWeight','Bold');
text(0,24,lab			  ,'FontSize',12,'FontWeight','Bold');

line([0 1],[22 22],'LineWidth',3,'Color',[0 0 0])

text(0.18,21,'corrected p value'  		,'FontSize',10);
text(0.42,21,'cluster-level {k}'		,'FontSize',10);
text(0.62,21,['voxel-level {' SPMdist '}'] 	,'FontSize',10);
text(0.86,21,'x,y,z {mm}'         		,'FontSize',10);

line([0 1],[20 20],'LineWidth',3,'Color',[0 0 0])


% Print cluster and maximum voxel-level p values {Z}
%-----------------------------------------------------------------------
text(0.18,19,sprintf('%-0.3f',Pu)		,'FontSize',10);
text(0.42,19,sprintf('%-0.3f (%i)',Pn,N)	,'FontSize',10);
text(0.62,19,sprintf('%-0.3f (%0.2f)',EC(1),Z)	,'FontSize',10);
text(0.84,19,sprintf('%-6.0f',L)		,'Fontsize',10);

line([0 1],[1 1],'LineWidth',3,'Color',[0 0 0])
