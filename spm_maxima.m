
% Tabular display of adjusted data
% FORMAT spm_maxima
%____________________________________________________________________________
%
% spm_maxima is called by spm_sections_ui and takes variables in working
% memory to produce a table of maxima within the selected region.  The region
% and its associated maxima are characterized in terms of its cluster,
% and voxel-level p values (corrected and uncorrected).
%
% Upto 16 maxima are displayed that are all at least 8mm apart.  If the maximum
% selected corresponds to one of these maxima, its location is displayed in red
%
%__________________________________________________________________________
% %W% %E%


% characterize point list in terms of maxima and regions
%----------------------------------------------------------------------------
[N Z M A] = spm_max(t,XYZ,V([4 5 6]));


% find nearest maximum [in a Euclidean sense] in the point list
%----------------------------------------------------------------------------
[d i] = min(sum(([(M(1,:) - L(1));(M(2,:) - L(2));(M(3,:) - L(3))]).^2));
L     = M(:,i);


% reset the pointer and position strings created by spm_sections_ui.m
%----------------------------------------------------------------------------
if V(3) == 1
	set(h1,'String',sprintf('%0.0f',L(1)));
	set(h2,'String',sprintf('%0.0f',L(2)));
	set(X1,'Position',[L(1)  L(2) 1]);
else
	set(h1,'String',sprintf('%0.0f',L(1)));
	set(h2,'String',sprintf('%0.0f',L(2)));
	set(h3,'String',sprintf('%0.0f',L(3)));
	set(X1,'Position',[(124 + L(2))  (248 + L(1)) 1]);
	set(X2,'Position',[(124 + L(2))  (112 - L(3)) 1]);
	set(X3,'Position',[(276 + L(1))  (112 - L(3)) 1]);
end

% select region and compute p values for all its maxima
%----------------------------------------------------------------------------
d         = A == A(i);
N         = N(d);
Z         = Z(d);
M         = M(:,d);


% delete previous axis
%----------------------------------------------------------------------------
subplot(2,1,2); delete(gca)
subplot(2,1,2); axis off


% display (sorted on Z values)
%===========================================================================


% table headings
%---------------------------------------------------------------------------
y   = 24;
axes('Position',[0.1 0.06 0.8 0.46]); axis off
text(0,y,['P values & statistics:   ' CWD],'FontSize',12,'FontWeight','Bold');
y   = y - 1;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y   = y - 1;

% Construct tables
%===========================================================================
text(0.18,y,'cluster-level {k,Z}','FontSize',10);
text(0.42,y,'voxel-level {Z}'    ,'FontSize',10);
text(0.64,y,'uncorrected'        ,'FontSize',10);
text(0.84,y,'location {mm}'      ,'FontSize',10);
y  = y - 1;

line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y  = y - 1;

% cluster-level p values {k}
%---------------------------------------------------------------------------
[u i] = max(Z);				% largest Z value
j     = find(A == A(i));		% maxima in cluster
Pk    = spm_P(1,W,U,N(i),S);		% cluster-level p value
Pu    = spm_P(1,W,u,0,S);		% voxel-level p value
Pz    = 1 - spm_Ncdf(u);		% uncorrected p value

% print cluster and maximum voxel-level p values {Z}
%---------------------------------------------------------------------------
str   = sprintf('%-0.3f   (%i, %0.2f)',Pk,N(i),u);
text(0.18,y,str,'FontSize',8,'FontWeight','Bold')
str   = sprintf('%-0.3f   (%-0.2f)',Pu,u);
text(0.44,y,str,'FontSize',8,'FontWeight','Bold')
str   = sprintf('%-0.3f',Pz);
text(0.68,y,str,'FontSize',8,'FontWeight','Bold')
 
d = text(0.82,y,sprintf('%-6.0f',M(:,i)),'Fontsize',10,'FontWeight','Bold');
if all(~(M(:,i) - L))
	set(d,'Color',[1 0 0]); end

y     = y - 1;

% print region and 4 secondary maxima (8mm apart)
%---------------------------------------------------------------------------
[l q] = sort(-Z);			% sort on Z value
D     = i;
for i = 2:length(q)
	d     =  min(sum((M(:,D) - M(:,q(i))*ones(1,size(D,2))).^2));
	if (d > 64 ) & (length(D) < 16);

		% voxel-level p values {Z}
		%-----------------------------------------------------------
		Pu  = spm_P(1,W,Z(q(i)),0,S);	% voxel-level p value
		Pz  = 1 - spm_Ncdf(Z(q(i)));	% uncorrected p value

		str = sprintf('%-0.3f   (%-0.2f)',Pu,Z(q(i)));
		text(0.44,y,str,'FontSize',8)
		str = sprintf('%-0.3f',Pz);
		text(0.68,y,str,'FontSize',8)
		str = sprintf('%-6.0f',M(:,q(i)));
		text(0.84,y,str,'FontSize',8)
		D   = [D q(i)];
		y   = y - 1;
	end
end

line([0 1],[y y],'LineWidth',1,'Color',[0 0 0])
y     = y - 1;

% footnote with SPM parameters
%---------------------------------------------------------------------------
str = sprintf('Height threshold {u} = %0.2f, p = %0.3f',U,1 - spm_Ncdf(U));
text(0,y,str,'FontSize',8);
y  = y - 1;
str = sprintf('Extent threshold {k} = %i voxels',k);
text(0,y,str,'FontSize',8);


set(gca,'Ylim',[y 24]);
