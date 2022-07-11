function res = repl(this, varargin)
% Method for getting replication counts, over trials
% FORMAT res = repl(this, index, nrepl)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = getset(this, 'trials', 'repl', varargin{:});
