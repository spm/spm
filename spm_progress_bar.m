% Display a 'Progress Bar'
% FORMAT spm_progress_bar('Init',height,xlabel,ylabel)
% Initialises the bar in the 'Interactive' window.
%
% FORMAT spm_progress_bar('Set',value)
% Sets the height of the bar itsself.
%
% FORMAT spm_progress_bar('Clear')
% Clears the 'Interactive' window.
%

% %W% John Ashburner %E%

function spm_progress_bar(action,arg1,arg2,arg3,arg4)
global pb_pointer pb_name ax
if (nargin == 0)
	spm_progress_bar('Init');
else
	if (strcmp(action,'Init'))
		if (nargin<4)
			arg3 = '';
			if (nargin<3)
				arg2 = 'Computing';
				if (nargin<2)
					arg1 = 1;
				end
			end
		end
		fg = spm_figure('FindWin','Interactive');
		pb_pointer = get(fg,'Pointer');
		pb_name    = get(fg,'Name');
		figure(fg);
		spm_progress_bar('Clear');
		set(fg,'Pointer','watch');
		set(fg,'Name',pb_name);
		ax = axes('Position', [0.45 0.2 0.1 0.6],...
			'XTick',[],...
			'Xlim', [0 1],...
			'Ylim', [0 arg1],...
			'Box', 'on');
		xlabel(arg2);
		ylabel(arg3);
		title('0% Complete');
		tim = clock;
		tim = tim(4:6);
		t1=text(1.1,arg1/2,0,sprintf('Began %2.0f:%2.0f:%2.0f',tim(1),tim(2),tim(3)));
		set(t1,'Tag','StartTime');
		line('Xdata',[0.5 0.5], 'Ydata',[0 0],...
			'LineWidth',16, 'Color', [1 0 0],'Tag','ProgressBar');
		drawnow;
	elseif (strcmp(action,'Set'))
		if (nargin<2)
			arg1 = 0;
		end
		F = spm_figure('FindWin','Interactive');
		br = findobj(F,'Tag','ProgressBar');
		if (~isempty(br))
			set(br,'Ydata',[0 arg1]);
			lim = get(get(br,'Parent'),'Ylim');lim=lim(2);
			title(sprintf('%.0f%% Complete',100*arg1/lim));
			drawnow;
			% fprintf(['%.1f%%' 8 8 8 8 8 8 8 ],100*arg1/lim);
		end
	elseif (strcmp(action,'Clear'))
		fg = spm_figure('FindWin','Interactive');
		spm_figure('Clear',fg);
		set(fg,'Pointer',pb_pointer);
		set(fg,'Name',pb_name);
		drawnow;
	end
end
