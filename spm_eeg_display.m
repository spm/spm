function D = spm_eeg_display(S)
% (internal) function to display epoched EEG/MEG channel data
% FORMAT D = spm_eeg_display(S)
%
% struct S must contain:
% Hfig		- handle of figure
% D			- EEG data struct
%
% The following variables are stored in D.gfx:
% 		Hflathead	- axes handle for flathead
% 		Heegaxes	- struct of axes handles for EEG plots
%_______________________________________________________________________
%
% spm_eeg_display is an internally used function that plots EEG/MEG traces.
%_______________________________________________________________________
% %W% Stefan Kiebel %E%


try
    Hfig = S.Hfig;
catch
	error('Argument ''Hfig'' is missing');
end
% D = get(Hfig, 'UserData');

try
    D = S.D;
catch
	error('Argument ''D'' is missing');
end

% display flathead with plots in it
if D.gfx.first == 1
    D.gfx.first = 0;
    % first call to spm_eeg_display
    figure(Hfig);
    Pos = [-0.05 0.34 0.95 0.62];

	% Compute width of display boxes
	Csetup = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));
	p = Csetup.Cpos(:, D.channels.order);
	Rxy = Csetup.Rxy;
	lp = size(p, 2);
	
	if lp > 1
		% more than 1 channel
		for i = 1:lp
			for j = 1:lp
				d(i,j) = sqrt(sum((p(:,j)-p(:,i)).^2));
				alpha(i,j) = acos((p(1,j)-p(1,i))/(d(i,j)+eps));
			end
		end
		d = d/2;
	
		alpha(alpha > pi/2) = pi-alpha(alpha > pi/2);
		Talpha = asin(1/(sqrt(1+Rxy^2)));

		for i = 1:lp
			for j = 1:lp
				if alpha(i,j) <= Talpha
					x(i,j) = d(i,j)*cos(alpha(i,j));
				else
					x(i,j) = Rxy*d(i,j)*cos(pi/2-alpha(i,j));
				end
			end
		end
        
        
		% half length of axes in x-direction
		Lxrec = min(x(find(x~=0)));
		
	else
		% 1 channel
		Lxrec = 0.3;
	end
	
	% l is used for putting the channel names into the axes
	l = Lxrec/10;

	% plot the graphs
	for i = 1:lp
		% small textboxes with channel names in lower left corner of each plot
		uicontrol('Style', 'text', 'Units', 'normalized', 'Position', ...
			[(p(1,i)-Lxrec+l)*Pos(3)+Pos(1) (p(2,i)-Lxrec/Rxy+l)*Pos(4)+Pos(2) 4*l l],...
			'String', D.channels.name{i}, 'BackgroundColor', get(gca, 'Color'),...
			'FontWeight', 'bold');
		
		% uicontextmenus for axes
		Heegmenus(i) = uicontextmenu;
		Heegmenus_add(i) = uimenu(Heegmenus(i), 'Label', sprintf('Add %s to graph', D.channels.name{i}),...
			'CallBack', 'spm_eeg_update_gfx(''addchannel'', gcbo)',...
			'UserData', i);
		
		Heegmenus_rm(i) = uimenu(Heegmenus(i), 'Label',...
			sprintf('Remove %s from graph', D.channels.name{i}),...
			'CallBack', 'spm_eeg_update_gfx(''rmchannel'', gcbo)', 'UserData', i);
		
		if ~ismember(i, D.gfx.Cdisplay)
			set(Heegmenus_rm(i), 'Enable', 'off');
		end
		
		Heegaxes(i) = axes('Box', 'on', 'Position',...
			[(p(1,i)-Lxrec)*Pos(3)+Pos(1) (p(2,i)-Lxrec/Rxy)*Pos(4)+Pos(2) 2*Lxrec*Pos(3) 2*Lxrec/Rxy*Pos(4)],...
			'Parent', Hfig, 'UIContextMenu', Heegmenus(i));
		
		
		for j = 1:length(D.gfx.Tdisplay)
			plot(D.data(i, :, D.gfx.Tdisplay(j)), D.gfx.linestyle{j});
		end
		
		set(gca, 'YLim', [-D.gfx.scale/2 D.gfx.scale/2],...
			'XLim', [1 D.Nsamples], 'XTick', [], 'YTick', []);
	end
	
	% construct graph with scalings
	% axes('Position', [0.05 0.85 2*Lxrec*Pos(3) 2*Lxrec/Rxy*Pos(4)]);
	axes('Position', [0.05 0.87 0.2 0.07]);
	set(gca, 'YLim', [-D.gfx.scale/2 D.gfx.scale/2], 'XLim', [1 D.Nsamples],...
		'XTick', [], 'YTick', [], 'LineWidth', 2, 'XAxisLocation', 'top');
	text(0, -D.gfx.scale/2, sprintf(' %d \\muV', D.gfx.scale), 'Interpreter', 'Tex',...
		'FontSize', 16, 'VerticalAlignment', 'bottom',...
		'Tag', 'scaletext2');
	text(D.Nsamples, D.gfx.scale/2, sprintf('%d ms', D.Nsamples), 'Interpreter', 'Tex',...
		'FontSize', 16, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
	
	D.gfx.Heegaxes = Heegaxes;
	D.gfx.Heegmenus = Heegmenus;
	D.gfx.Heegmenus_add = Heegmenus_add;
	D.gfx.Heegmenus_rm = Heegmenus_rm;
	D.gfx.Lxrec = Lxrec;
	
else
	
	% plot the graphs
	for i = 1:D.Nchannels
		axes(D.gfx.Heegaxes(i));
		cla
		set(D.gfx.Heegaxes(i), 'NextPlot', 'add');
		
		for j = 1:length(D.gfx.Tdisplay)
			plot(D.data(i, :, D.gfx.Tdisplay(j)), D.gfx.linestyle{j});
		end
		
		set(gca, 'YLim', [-D.gfx.scale/2 D.gfx.scale/2],...
			'XLim', [1 size(D.data, 2)]);
	end

	
end
set(Hfig, 'UserData', D);
