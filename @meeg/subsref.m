function varargout=subsref(this,subs)
% SUBSREF Subscripted reference
% An overloaded function...
% _________________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Stefan Kiebel
% $Id: subsref.m 6600 2015-11-12 13:07:41Z christophe $

if isempty(subs)
    return;
end

if this.Nsamples == 0
    error('Attempt to reference a field of an empty meeg object.');
end

switch subs(1).type
    case '()'
        if ~islinked(this), error('The object is not linked to data file'); end
        if numel(subs)~= 1, error('Expression too complicated');            end

        if this.montage.Mind==0
            varargout = {double(subsref(this.data, subs))};
        else
            Mem_max = 200*2^20; % Limit memory usage to about 200Mb
            vect_fl = 0; dat3D = 0;
            dim = size(this);
            if numel(dim)>2 && dim(3)>1
                % assume at most 3D data
                dat3D = 1;
            end
            if ischar(subs.subs{1})
                % need to handle the case of a ':' argument
                if ~strcmp(subs.subs{1},':'), 
                    error('This shouldn''t happen....'); 
                end
                if length(subs.subs) == 1
                    chanidx = 1:dim(1); 
                    subs.subs{2} = ':'; 
                    if dat3D, subs.subs{3} = ':'; end
                    vect_fl = 1;
                else
                    chanidx = 1:dim(1);
                end
            else
                chanidx = subs.subs{1};
            end
            
            % check if correct channel index
            if any(chanidx > nchannels(this))
                error('channel index higher than number of channels in current montage')
            end
            % get corresponding rows of 'tra' matrix
            traidx = this.montage.M(this.montage.Mind).tra(chanidx,:);
            % change subs to use only the necessary channels from data
            lchan_o = find(any(traidx,1));
            subs.subs{1} = lchan_o;
            
            % need to handle the case of ':' arguments
            if ischar(subs.subs{2})
                if ~strcmp(subs.subs{2},':'), error('This shouldn''t happen....'); end
                subs.subs{2} = 1:dim(2);
                Ntb = dim(2);
            else
                Ntb = length(subs.subs{2});
            end
            if dat3D && ischar(subs.subs{3})
                if ~strcmp(subs.subs{3},':'), error('This shouldn''t happen....'); end
                subs.subs{3} = 1:dim(3);
            end
            Mem_load = 8*Ntb*length(lchan_o);
            if Mem_load<=Mem_max % small chunk loaded
                if dat3D
                    subs_c = subs;
                    for ii=1:numel(subs.subs{3})
                        subs_c.subs{3} = subs.subs{3}(ii);
                        varargout{1}(:,:,ii) = ...
                            traidx(:,lchan_o)*double(subsref(this.data, subs_c));
                    end
                else
                    varargout = {traidx(:,lchan_o)*double(subsref(this.data, subs))};
                end
            else % otherwise split data reading into chunks
                Ntb_chunk = round(Mem_max/length(lchan_o)/8);
                Nchunk = ceil(Ntb/Ntb_chunk);
                varargout{1} = zeros(length(chanidx),Ntb);
                for ii=1:Nchunk
                    subs_ch = subs;
                    if ii<Nchunk
                        ll = (1:Ntb_chunk)+(ii-1)*Ntb_chunk;
                    else
                        ll = ((ii-1)*Ntb_chunk+1):Ntb;
                    end
                    subs_ch.subs{2} = subs_ch.subs{2}(ll);
                    if dat3D
                        for jj=1:numel(subs.subs{3})
                            subs_ch.subs{3} = subs.subs{3}(jj);
                            varargout{1}(:,ll,jj) = ...
                                traidx(:,lchan_o)*double(subsref(this.data, subs_ch));
                        end
                    else
                        varargout{1}(:,ll) = ...
                            traidx(:,lchan_o)*double(subsref(this.data, subs_ch));
                    end
                end
            end
            if vect_fl, varargout{1} = varargout{1}(:); end
        end

    case '{}'
    case '.'
        if ismethod(this, subs(1).subs)
            if numel(subs) == 1
                varargout = {feval(subs(1).subs, this)};
            elseif (numel(subs) == 2) && isequal(subs(2).type,  '()')
                varargout = {feval(subs(1).subs, this, subs(2).subs{:})};
            elseif (numel(subs)> 2) && isequal(subs(2).type,  '()')
                varargout{1} = builtin('subsref', ...
                    feval(subs(1).subs, this, subs(2).subs{:}),  subs(3:end));
            else
                varargout{1} = builtin('subsref', feval(subs(1).subs, this),  subs(2:end));
            end
        elseif isfield(this.other, subs(1).subs)
            field = getfield(this.other, subs(1).subs);
            if numel(subs)==1
                varargout = {field};
            else
                varargout{1} = builtin('subsref', field, subs(2:end));
            end
        else
            error('Reference to non-existent or private meeg method or field.');
        end
    otherwise
        error('Unfamiliar referencing type');
end
