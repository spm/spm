
% Tabular display of adjusted data
% FORMAT spm_maxima
%_______________________________________________________________________
%
% spm_maxima is called by spm_sections_ui and takes variables in working
% memory to produce a table of maxima within the selected region.  The region
% and its associated maxima are characterized in terms of its cluster,
% and voxel-level p values (corrected and uncorrected).
%
% The maxima displayed are all at least 8mm apart.  If the maximum
% selected corresponds to one of these maxima, its location is
% displayed in red.  Selecting the voxel coordinates of other 
% maxima causes a green pointer to appear in the appropriate place on
% maximum intensity projection.
%
%_______________________________________________________________________
% %W% %E%


% characterize point list in terms of maxima and regions
%-----------------------------------------------------------------------
[N Z M A] = spm_max(t,XYZ,V([4 5 6]));


% find nearest maximum [in a Euclidean sense] in the point list
%-----------------------------------------------------------------------
[d i] = min(sum(([(M(1,:) - L(1));(M(2,:) - L(2));(M(3,:) - L(3))]).^2));
L     = M(:,i);


% reset the pointer and position strings created by spm_sections_ui.m
%-----------------------------------------------------------------------
if V(3) == 1
	set(h1,'String',sprintf('%0.0f',L(1)));
	set(h2,'String',sprintf('%0.0f',L(2)));
	set(hXstr,'String',sprintf('x = %0.0f',L(1)));
	set(hYstr,'String',sprintf('y = %0.0f',L(2)));
	set(X1,'Position',[L(1)  L(2) 1]);
else
	set(h1,'String',sprintf('%0.0f',L(1)));
	set(h2,'String',sprintf('%0.0f',L(2)));
	set(h3,'String',sprintf('%0.0f',L(3)));
	set(hXstr,'String',sprintf('x = %0.0f',L(1)));
	set(hYstr,'String',sprintf('y = %0.0f',L(2)));
	set(hZstr,'String',sprintf('z = %0.0f',L(3)));

	set(X1,'Position',[(P1 + L(2))  (P2 + L(1)) 1]);
	set(X2,'Position',[(P1 + L(2))  (P3 - L(3)) 1]);
	set(X3,'Position',[(P4 + L(1))  (P3 - L(3)) 1]);
end

% select region and compute p values for all its maxima
%-----------------------------------------------------------------------
d         = A == A(i);
N         = N(d);
Z         = Z(d);
M         = M(:,d);


% delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
subplot(2,1,2); delete(gca), spm_figure('DeletePageControls')
subplot(2,1,2); axis off


% display (sorted on Z)
%=======================================================================


% table headings
%-----------------------------------------------------------------------
axes('Position',[0.1 0.06 0.8 0.46]); axis off
text(0,24,['Statistics:  ',spm('DirTrunc',CWD)],...
	'FontSize',12,'FontWeight','Bold');
line([0 1],[23 23],'LineWidth',3,'Color',[0 0 0])

% Construct tables
%=======================================================================
if SPMZ

	text(0.18,22,'cluster-level {k,Z}','FontSize',10);
	text(0.42,22,'voxel-level {Z}'    ,'FontSize',10);
	text(0.62,22,'uncorrected k & Z'  ,'FontSize',10);
	text(0.86,22,'x,y,z {mm}'         ,'FontSize',10);

elseif SPMF

	text(0.18,22,'cluster-level {k,F}','FontSize',10);
	text(0.42,22,'voxel-level {F}'    ,'FontSize',10);
	text(0.62,22,'uncorrected k & F'  ,'FontSize',10);
	text(0.86,22,'x,y,z {mm}'         ,'FontSize',10);


end

line([0 1],[21 21],'LineWidth',3,'Color',[0 0 0])

% Pagination variables
%-----------------------------------------------------------------------
y     = 20;
hPage = [];
nPage = 1;
Vis   = 'on';

% cluster-level p values {k}
%-----------------------------------------------------------------------
[u i] = max(Z);				% largest Z value

if SPMZ

Pk    = spm_P(1,W,U,N(i),S);		% cluster-level p value
Pkn   = spm_Pkn(N(i),(u - U),W,u,S);	% cluster-level (Bivariate)
Pu    = spm_P(1,W,u,0,S);		% voxel-level p value
Pn    = 1 - spm_kcdf(N(i),U,W);		% uncorrected p value (k)
Pz    = 1 - spm_Ncdf(u);		% uncorrected p value (Z)

elseif SPMF

Pk    = NaN;				% cluster-level p value
Pkn   = NaN;				% cluster-level (Bivariate)
Pu    = spm_pF(S,W,df,u);		% voxel-level p value
Pn    = NaN;				% uncorrected p value (k)
Pz    = 1 - spm_Fcdf(u,df);		% uncorrected p value (F)

end

% print cluster and maximum voxel-level p values {Z}
%-----------------------------------------------------------------------
str   = sprintf('%-0.3f   (%i, %0.2f)',Pkn,N(i),u);
h     = text(0.18,y,str,'FontSize',8,'FontWeight','Bold');
hPage = [hPage, h];
str   = sprintf('%-0.3f   (%-0.2f)',Pu,u);
h     = text(0.44,y,str,'FontSize',8,'FontWeight','Bold');
hPage = [hPage, h];
str   = sprintf('%-0.3f',Pn);
h     = text(0.64,y,str,'FontSize',8,'FontWeight','Bold');
hPage = [hPage, h];
str   = sprintf('%-0.3f',Pz);
h     = text(0.74,y,str,'FontSize',8,'FontWeight','Bold');
hPage = [hPage, h];
h     = text(0.84,y,sprintf('%-6.0f',M(:,i)),...
	'Fontsize',8,'FontWeight','Bold','UserData',M(:,i));
hPage = [hPage, h];
if all(~(M(:,i) - L)), set(h,'Color','r','FontAngle','Italic'); end

% Callback for voxel co-ordinates
%-----------------------------------------------------------------------
tmp = [	'l = get(gco,''UserData'');',...
	'if V(3)==1,',...
		'set(X1g,''Position'',[l(1) l(2) 1],''Visible'',''on''),',...
	'else,',...
		'set(X1g,''Position'',[(P1 + l(2))  (P2 + l(1)) 1]),',...
		'set(X2g,''Position'',[(P1 + l(2))  (P3 - l(3)) 1]),',...
		'set(X3g,''Position'',[(P4 + l(1))  (P3 - l(3)) 1]),',...
		'set([X1g X2g X3g],''Visible'',''on''),',...
	'end,',...
	'set(gcf,''WindowButtonUpFcn'',[',...
		'''set([X1g,X2g,X3g],''''Visible'''',''''off''''),',...
		'set(gcf,''''WindowButtonUpFcn'''','''''''')''])'];

set(h,'ButtonDownFcn',tmp,'Interruptible','no')

y     = y - 1;

% print region and all secondary maxima (8mm apart)
%-----------------------------------------------------------------------
[l q] = sort(-Z);			% sort on Z value
D     = i;
for i = 2:length(q)
	d     =  min(sum((M(:,D) - M(:,q(i))*ones(1,size(D,2))).^2));
	if (d > 64 )

		% Paginate if necessary
		%-------------------------------------------------------
		if y < 2
			h     = text(0.5,-2,['Page ',num2str(nPage)],...
				'FontSize',8,'FontAngle','Italic',...
				'Visible',Vis);
			hPage = [hPage, h];
			spm_figure('NewPage',hPage)
			hPage = [];
			nPage = nPage + 1;
			Vis   = 'off';
			y     = 20;
		end

		% voxel-level p values {Z}
		%-------------------------------------------------------
		if SPMZ

		Pu    = spm_P(1,W,Z(q(i)),0,S);		% voxel-level p value
		Pz    = 1 - spm_Ncdf(Z(q(i)));		% uncorrected p value

		elseif SPMF

		Pu    = spm_pF(S,W,df,Z(q(i)));		% voxel-level p value
		Pz    = 1 - spm_Fcdf(Z(q(i)),df);	% uncorrected p value

		end

		str   = sprintf('%-0.3f   (%-0.2f)',Pu,Z(q(i)));
		h     = text(0.44,y,str,'FontSize',8,'Visible',Vis);
		hPage = [hPage, h];
		str   = sprintf('%-0.3f',Pz);
		h     = text(0.74,y,str,'FontSize',8,'Visible',Vis);
		hPage = [hPage, h];
		str   = sprintf('%-6.0f',M(:,q(i)));
		h     = text(0.84,y,str,'FontSize',8,'Visible',Vis,...
			'UserData',M(:,q(i)));
		set(h,'ButtonDownFcn',tmp,'Interruptible','no')
		if all(~(M(:,q(i)) - L))
			set(h,'Color','r','FontAngle','Italic'); end
		hPage = [hPage, h];
		D     = [D q(i)];
		y     = y - 1;
	end
end

% number and paginate last page (if pagination was required)
%-----------------------------------------------------------------------
if strcmp(Vis,'off')
	h     = text(0.5,-2,['Page ',num2str(nPage)],...
		'FontSize',8,'FontAngle','Italic',...
		'Visible',Vis);
	hPage = [hPage, h];
	spm_figure('NewPage',hPage);
end
line([0 1],[1 1],'LineWidth',1,'Color',[0 0 0])

