
% Tabular display of adjusted data
% FORMAT spm_maxima
%_______________________________________________________________________
%
% spm_maxima is called by spm_results and takes variables in
% working memory to produce a table of maxima within the selected
% region.  The region and its associated maxima are characterized in
% terms of its cluster, and voxel-level p values (corrected and
% uncorrected).
%
% The maxima displayed are all at least 8mm apart.  If the maximum
% selected corresponds to one of these maxima, its location is
% displayed in red italics.  Selecting the voxel coordinates of a
% maxima causes a green pointer to appear in the appropriate place on
% maximum intensity projection.
%
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%


% Characterize point list in terms of maxima and regions
%-----------------------------------------------------------------------
[N Z ML A] = spm_max(t,XYZ,V([4 5 6]));

% Find nearest maximum [in a Euclidean sense] in point list & update GUI
%-----------------------------------------------------------------------
[L,i,d] = spm_XYZreg('NearestXYZ',L,ML);
spm_XYZreg('SetCoords',L,hReg);

% Select region and compute p values for all its maxima
%-----------------------------------------------------------------------
d     = A == A(i);
N     = N(d);
Z     = Z(d);
ML    = ML(:,d);


% Display (sorted on Z/F)
%=======================================================================

%-Delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
spm_results_ui('ClearPane',Fgraph)

% Table headings
%-----------------------------------------------------------------------
axes('Position',[0.1 0.06 0.8 0.46]); axis off
text(0,24,['Statistics:  ',spm_str_manip(CWD,'a50')],...
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

% Cluster-level p values {k}
%-----------------------------------------------------------------------
[u i] = max(Z);				% largest Z value

if SPMZ

Pk    = spm_P(1,W,U,N(i),S);		% cluster-level p value
Pkn   = spm_Pkn(N(i),(u - U),W,U,S);	% cluster-level (Bivariate)
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

% Print cluster and maximum voxel-level p values {Z}
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
h     = text(0.84,y,sprintf('%-6.0f',ML(:,i)),...
	'Fontsize',8,'FontWeight','Bold',...
	'ButtonDownFcn','spm_mip_ui(''ShowGreens'')',...
	'Interruptible','off',...
	'UserData',ML(:,i));
hPage = [hPage, h];
if all(~(ML(:,i) - L)), set(h,'Color','r','FontAngle','Italic'); end

y     = y - 1;

% Print region and all secondary maxima (8mm apart)
%-----------------------------------------------------------------------
[l q] = sort(-Z);			% sort on Z value
D     = i;
for i = 2:length(q)
	d     =  min(sum((ML(:,D) - ML(:,q(i))*ones(1,size(D,2))).^2));
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
		str   = sprintf('%-6.0f',ML(:,q(i)));
		h     = text(0.84,y,str,'FontSize',8,...
				'Visible',Vis,...
				'ButtonDownFcn',...
					'spm_mip_ui(''ShowGreens'')',...
				'Interruptible','off',...
				'UserData',ML(:,q(i)));
		if all(~(ML(:,q(i)) - L))
			set(h,'Color','r','FontAngle','Italic'); end
		hPage = [hPage, h];
		D     = [D q(i)];
		y     = y - 1;
	end
end

% Number and paginate last page (if pagination was required)
%-----------------------------------------------------------------------
if strcmp(Vis,'off')
	h     = text(0.5,-2,['Page ',num2str(nPage)],...
		'FontSize',8,'FontAngle','Italic',...
		'Visible',Vis);
	hPage = [hPage, h];
	spm_figure('NewPage',hPage);
end
line([0 1],[1 1],'LineWidth',1,'Color',[0 0 0])

