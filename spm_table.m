
% Tabular display of adjusted data
% FORMAT spm_table
%____________________________________________________________________________
%
% spm_table is called by spm_sections_ui and takes variables in working
% memory to produce a table of adjusted mean activities for the voxel
% in question
% These mean activities are the average estimates over subjects.  If these
% estimates derive from an activation study (the effects are factors
% described by the H partition) then an effect size is calculated based
% on the specified contrast (e.g. the activation in rCBF equivalents)
%
%__________________________________________________________________________
% %W% %E%


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
Y      = XA(:,i);
Z      = t(i);					% Z value
Pz     = spm_Pz(W,Z,S);				% corrected p value
Pu     = 1 - spm_Ncdf(Z);			% uncorrected p value


% condition effects or covariate ?
%----------------------------------------------------------------------------
d      = min(find(CONTRAST(con,:)));

% delete previous axis
%----------------------------------------------------------------------------
subplot(2,1,2); delete(gca)
subplot(2,1,2); axis off
text(0,1.1,CWD,'FontSize',16,'FontWeight','bold');
text(0,1.0,sprintf('Location {x,y,z}  =  %0.0f %0.0f %0.0f mm',L'),...
'FontSize',12,'FontWeight','bold');
str    = sprintf('Z = %0.3f;   P(Zmax > Z) = %0.2f (corrected) and %0.2e (uncorrected)',Z,Pz,Pu);
text(0,0.9,str,'FontSize',12,'FontWeight','bold');


if d > (size(H,2) + size(K,2));	% covariate

	if size(B,2)
		y     = 0;
		for i = 1:size(B,2)
			y = y + Y(find(B(:,i) > 0));
		end
		Y     = y/size(B,2);
	end
	text(0,0.8,sprintf('[Adjusted] activity (scans %0.0f - %0.0f )',...
	[1 length(Y)]),'FontSize',12,'FontWeight','bold');

else				% condition means

	% get voxel-specific adjusted data 
	%-------------------------------------------------------------------
	y     = [];
	for j = 1:size(H,2); y(j) = mean(Y(H(:,j))); end
	Y     = y;
	d     = CONTRAST(con,:);
	d     = d/sum(d(d > 0));
	ES    = Y*d([1:size(H,2) + size(K,2)])';
	str   = ...
	sprintf('Adjusted activity (conditions %0.0f - %0.0f)',1,length(Y));
	str   = ...
	[str sprintf('  Effect size: %0.2f (%0.2f percent)',ES,ES*100/mean(Y))];
	text(0,0.8,str,'FontSize',12,'FontWeight','bold');
end

% display data
%----------------------------------------------------------------------------
y     = 0.7;
x     = 0;
for i = 1:length(Y)
	text(x,y,sprintf('%0.0f  %-12.2f',i,Y(i)),'FontSize',10)
	y = y - 0.05;
	if y < -.24; y = 0.7; x = x + 0.16; end
end
