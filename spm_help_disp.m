function spm_help_disp(s)
% displays commented text in a graphics window
% FORMAT spm_help_disp(s)
% s   - filename {e.g. spm_routine.m}
%____________________________________________________________________________
%
% spm_help_disp reads an ASCII file until an empty line is encountered
% If the line begins with '%' it is displayed.  The second line is displayed
% in bold slightly above the remainder (the first line is ignored)
%
% This facility is specifically for displaying the commented headers
% of MATLAB scipts (and suitably configured manual pages) in a figure
% that has 'Tag' = 'Graphics' (if it exists).
%
%__________________________________________________________________________
% %W% %E%

%----------------------------------------------------------------------------
h   = findobj(get(0,'Children'),'flat','Tag','Graphics');
if isempty(h); h = gcf; end
figure(h)
spm_clf; axis off

fid = fopen(s,'r');
S   = setstr(fread(fid))';
q   = min([length(S) findstr(S,setstr([10 10]))]);	% find empty lines
q   = find(S(1:q(1)) == 10);				% find line breaks
y   = 0.84;

text(-0.1,0.94,[s ' - help'],'FontSize',16,'FontWeight','bold');

% loop over lines of text
%----------------------------------------------------------------------------
for i = 1:(length(q) - 1)
	d = S((q(i) + 1):(q(i + 1) - 1));
	if d(1) == abs('%');
		h = text(-0.1,y,d(2:length(d)),'FontSize',10);
		if y == 0.84	% if the first line displayed
			set(h,'FontWeight','bold');
			y = 0.82;
		end
		y = y - 0.02;
	end
end


% create control object if scrolling is required
%----------------------------------------------------------------------------
if y < -0.1;
	uicontrol(gcf,'Style','Pushbutton','String','next',...
	'Callback',['spm_position(gca,1,[0 0.8]); refresh'],...
        'Position',[500 050 60 020]);
	uicontrol(gcf,'Style','Pushbutton','String','last',...
	'Callback',['spm_position(gca,1,[0 -.8]); refresh'],...
        'Position',[500 020 60 020]);
end

fclose(fid);
