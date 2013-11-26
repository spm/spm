function cfg_dbstop(fh)
fhinfo = functions(fh);
try
    switch fhinfo.type
        case 'simple'
            dbstop(fhinfo.function);
        case 'scopedfunction'
            evalc(sprintf('dbstop in %s at %s', fhinfo.parentage{2}, fhinfo.parentage{1}));
        case 'nested'
            fnames = regexp(fhinfo.function, '^(?<func>[^/]+)/(?<nfunc>.*)$', 'names');
            evalc(sprintf('dbstop in %s at %s', fnames.func, fnames.nfunc));
        case 'anonymous'
            fnames = regexp(func2str(fh), '@\(.*\)(?<afunc>[^\(]*)\(.*\)','names');
            if any(exist(fnames.afunc) == [2:6 8])
                dbstop(fnames.afunc);
            else
                [~, func] = fileparts(fhinfo.file);
                evalc(sprintf('dbstop in %s at %s', func, fnames.afunc));
            end
        otherwise
            error('matlabbatch:cfg_dbstop:noBreakPointSet', 'Unknown function type.');
    end
catch
    warning('matlabbatch:cfg_dbstop:noBreakPointSet', 'Don''t know how to set a breakpoint for ''%s''.\n', func2str(fh));
    disp(fhinfo);
end