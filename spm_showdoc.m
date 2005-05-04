function spm_showdoc(c)
% Show SPM documentation
%
% This function extracts and displays SPM documentation.
%
% Example usage:
%   diary('spmdoc.txt'); spm_showdoc; diary off
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_showdoc.m 112 2005-05-04 18:20:52Z john $


if nargin==0,
	c = spm_config;
end;
fprintf('\n\nSPM DOCUMENTATION\n\nTable of Contents\n----- -- --------\n');
contents(c);
fprintf('%s\n%s\n\n\n',repmat('_',1,80),repmat('_',1,80));
maintext(c);
return;

function contents(c,lev)
if nargin==1,
	lev = '';
end;
if isfield(c,'name'),
	str = [lev repmat(' ',1,25-length(lev)) c.name];
	fprintf('%s\n',str);
        i = 0;
	if isfield(c,'values'),
		for ii=1:length(c.values),
			if isstruct(c.values{ii}) && isfield(c.values{ii},'name'),
				i    = i+1;
				lev1 = sprintf('%s%d.', lev, i);
				contents(c.values{ii},lev1);
			end;
		end;
	end;
	if isfield(c,'val'),
		for ii=1:length(c.val),
			if isstruct(c.val{ii}) && isfield(c.val{ii},'name'),
				i    = i+1;
				lev1 = sprintf('%s%d.', lev, i);
				contents(c.val{ii},lev1);
			end;
		end;
	end;
end;
return;

function maintext(c, lev)
if nargin==1,
	lev = '';
end;
if ~isempty(lev) && sum(lev=='.')==1,
	disp(repmat('_',1,80));
	fprintf('\n');
end;
if isfield(c,'name'),
	str   = sprintf('%s %s', lev, c.name);
	under = repmat('-',1,length(str));
	fprintf('%s\n%s\n', str, under);
	if isfield(c,'modality'),
		fprintf('Only for');
		for i=1:numel(c.modality),
			fprintf(' %s',c.modality{i});
		end;
		fprintf('\n\n');
	end;
	if isfield(c,'help');
		hlp = spm_justify(80,c.help);
		disp(strvcat(hlp{:}));
	end;

	switch (c.type),
	case {'repeat'},
		if length(c.values)==1,
			fprintf('\nRepeat "%s", any number of times.\n',c.values{1}.name);
        else
			fprintf('\nAny of the following options can be chosen, any number of times\n');
			i = 0;
			for ii=1:length(c.values),
				if isstruct(c.values{ii}) && isfield(c.values{ii},'name'),
					i    = i+1;
					fprintf('    %2d) %s\n', i,c.values{ii}.name);
				end;
			end;
		end;
		fprintf('\n');

	case {'choice'},
		fprintf('\nAny one of these options can be selected:\n');
		i = 0;
		for ii=1:length(c.values),
			if isstruct(c.values{ii}) && isfield(c.values{ii},'name'),
				i    = i+1;
				fprintf('    %2d) %s\n', i,c.values{ii}.name);
			end;
		end;
		fprintf('\n');
	case {'branch'},
		fprintf('\nThis item contains %d fields:\n', length(c.val));
		i = 0;
		for ii=1:length(c.val),
			if isstruct(c.val{ii}) && isfield(c.val{ii},'name'),
				i    = i+1;
				fprintf('    %2d) %s\n', i,c.val{ii}.name);
			end;
		end;
		fprintf('\n');
	case {'menu'},
		fprintf('\nOne of these values is chosen:\n');
		for k=1:length(c.labels),
			fprintf('    %2d) %s\n', k, c.labels{k});
		end;
		fprintf('\n');
	case {'files'},
		if length(c.num)==1 && isfinite(c.num(1)) && c.num(1)>=0,
			fprintf('\nA "%s" file is selected by the user.\n',c.filter);
        else
			fprintf('\n"%s" files are selected by the user.\n',c.filter);
		end;

	case {'entry'},
		switch c.strtype,
		case {'e'},
			d = 'Evaluated statements';
		case {'n'},
			d = 'Natural numbers';
		case {'r'},
			d = 'Real numbers';
		case {'w'},
			d = 'Whole numbers';
		otherwise,
			d = 'Values';
		end;
		fprintf('\n%s are typed in by the user.\n',d);
	end;

	i = 0;
	fprintf('\n');
	if isfield(c,'values'),
		for ii=1:length(c.values),
			if isstruct(c.values{ii}) && isfield(c.values{ii},'name'),
				i    = i+1;
				lev1 = sprintf('%s%d.', lev, i);
				maintext(c.values{ii},lev1);
			end;
		end;
	end;
        if isfield(c,'val'),
                for ii=1:length(c.val),
                        if isstruct(c.val{ii}) && isfield(c.val{ii},'name'),
                                i    = i+1;
                                lev1 = sprintf('%s%d.', lev, i);
                                maintext(c.val{ii},lev1);
                        end;
                end;
        end;
	fprintf('\n');
end;
