function res = repl(obj, ind)
% Method for getting all replication counts, over trials
% FORMAT res = repl(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: repl.m 1373 2008-04-11 14:24:03Z spm $

if isfield(obj.trials(1), 'repl')
    res = cat(1, obj.trials.repl);
else
    res = zeros(length(obj.trials), 1);
end

if nargin > 1
    res = res(ind);
end
