
% Tabular display of data pertaining the SPM{F} and parameter estimates
% FORMAT spm_spmF_table
%____________________________________________________________________________
%
% spm_spmF_table is called by spm_spmF_ui and takes variables in working
% memory to produce a table of parameter estimates for s selected voxel
%
%__________________________________________________________________________
% %W% %E%

% find nearest voxel [in a Euclidean sense] in the point list of locations XYZ
%----------------------------------------------------------------------------
[d i] = min(sum(([(XYZ(1,:) - L(1));(XYZ(2,:) - L(2));(XYZ(3,:) - L(3))]).^2));
L     = XYZ(:,i);


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

% Z score and P(Zmax > u)
%----------------------------------------------------------------------------
Y      	= BETA(:,i);				% parameter estimates
PF     	= 1 - spm_Fcdf(SPMF(i),df);		% uncorrected p value
CPF 	= spm_pF(S, Wresid, Fdf, SPMF(i));	% Corrected p value

% delete previous axis
%----------------------------------------------------------------------------
subplot(2,1,2); delete(gca)
subplot(2,2,3);
imagesc([H C G B]*inv(diag(max(abs([H C G B])))))
title('Design Matrix')
xlabel('parameter')
ylabel('observation')

PWD    = pwd;
PWD    = PWD([max(find(pwd == '/')):length(pwd)]);
subplot(2,2,4); axis off
text(0,1.1,PWD,'FontSize',16,'FontWeight','bold');
text(0,1.0,sprintf('Location {x,y,z}  =  %0.0f %0.0f %0.0f mm',L'),...
'FontSize',12,'FontWeight','bold');
str    = sprintf('F = %0.3f df: %0.0f, %0.0f',SPMF(i),df(1),df(2));
text(0,0.9,str,'FontSize',12,'FontWeight','bold');
str    = sprintf('p = %0.3f (uncorrected) ',PF);
text(0,0.8,str,'FontSize',12);
str    = sprintf('p = %0.3f  (corrected) ',CPF);
text(0,0.7,str,'FontSize',12,'FontWeight','bold');

text(0,0.6,sprintf('Estimated parameters {%0.0f - %0.0f}',...
[1 length(Y)]),'FontSize',12,'FontWeight','bold');


% display data
%----------------------------------------------------------------------------
y     = 0.55;
x     = 0;
for i = 1:length(Y)
	text(x,y,sprintf('%0.0f  -  %-12.2f',i,Y(i)),'FontSize',10)
	y = y - 0.06;
	if y < 0; y = 0.55; x = x + 0.4; end
end
