
% Tabular display of adjusted data
% FORMAT spm_maxima
%____________________________________________________________________________
%
% spm_maxima is called by spm_sections_ui and takes variables in working
% memory to produce a table of maxima within the selected region.  The region
% and its associated maxima are characterized in terms of their corrected
% p values based on spatial extent [region] and height [maxima] and their
% uncorrected p values [maxima].
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
Pn        = spm_Pn(N,W,U,S);
Pz        = spm_Pz(W,Z,S);
Pu        = 1 - spm_Ncdf(Z);


% delete previous axis
%----------------------------------------------------------------------------
subplot(2,1,2); delete(gca)
subplot(2,1,2); axis off


% display (sorted on Z values)
%===========================================================================

% table headings
%---------------------------------------------------------------------------
y  = 24;
text(0,y,['Regional maxima:    ' CWD],'Fontsize',16,'FontWeight','Bold');
y  = y - 1.2;
line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y  = y - 1;

%---------------------------------------------------------------------------
text(0.00,y,'region');
text(0.10,y,'size {k}');
text(0.24,y,'P(n      > k)');
text(0.29,y,'max','Fontsize',8);
text(0.42,y,'Z');
text(0.48,y,'P(Z      > u)');
text(0.53,y,'max','Fontsize',8);
text(0.64,y,'(Uncorrected)','Fontsize',10);
text(0.84,y,'{x,y,z mm}');
y  = y - 1;

line([0 1],[y y],'LineWidth',3,'Color',[0 0 0])
y  = y - 1;

% list of maxima
%---------------------------------------------------------------------------
[j i] = max(Z);				% largest Z value

% print region and largest maximum
%-------------------------------------------------------------------
text(0.00,y,sprintf('%-0.0f',1),     'Fontsize',10,'FontWeight','Bold')
text(0.10,y,sprintf('%-0.0f',N(i))  ,'Fontsize',10,'FontWeight','Bold')
text(0.24,y,sprintf('%-0.3f',Pn(i)) ,'Fontsize',10,'FontWeight','Bold')
text(0.42,y,sprintf('%-0.2f',Z(i))  ,'Fontsize',10,'FontWeight','Bold')
text(0.54,y,sprintf('%-0.3f',Pz(i)) ,'Fontsize',10,'FontWeight','Bold')
text(0.64,y,sprintf('(%-0.3f)',Pu(i)) ,'Fontsize',10)
d = text(0.82,y,sprintf('%-6.0f',M(:,i)),'Fontsize',10,'FontWeight','Bold');
if all(~(M(:,i) - L))
	set(d,'Color',[1 0 0]); end

y     = y - 1;

% print region and 4 secondary maxima (8mm apart)
%-------------------------------------------------------------------
[l k] = sort(-Z);			% sort on Z value
D     = i;
for i = 2:length(k)
	d     =  min(sum((M(:,D) - M(:,k(i))*ones(1,size(D,2))).^2));
	if (d > 64 ) & (length(D) < 16);
		text(0.42,y,sprintf('%-0.2f',Z(k(i)))       ,'Fontsize',10);
		text(0.54,y,sprintf('%-0.3f',Pz(k(i)))      ,'Fontsize',10);
		text(0.64,y,sprintf('(%-0.3f)',Pu(k(i)))    ,'Fontsize',10);
		d = text(0.82,y,sprintf('%-6.0f',M(:,k(i))) ,'Fontsize',10);
		if all(~(M(:,k(i)) - L))
			set(d,'Color',[1 0 0]); end
		D = [D k(i)];
		y = y - 1;
  end
end

line([0 1],[y y],'LineWidth',1,'Color',[0 0 0])
y     = y - 1;

% volume, resels and smoothness 
%---------------------------------------------------------------------------
D     = length(W);					% dimension of SPM
RESEL = S*prod(V([1:D] + 3))/prod(FWHM);		% RESELS

% footnote with SPM parameters
%---------------------------------------------------------------------------
d  = sprintf('Threshold = %0.2f; Volume [S] = %0.0f voxels',U,S);
text(0,y,d,'Fontsize',10);
y  = y - 1;
if D == 2
    d = sprintf('FWHM = [%0.1f %0.1f] mm (i.e. %0.0f RESELS) ',FWHM,RESEL);
end
if D == 3
    d = sprintf('FWHM = [%0.1f %0.1f %0.1f] mm (i.e. %0.0f RESELS)',FWHM,RESEL);
end
text(0,y,d,'Fontsize',10);


set(gca,'Ylim',[y 24]);
