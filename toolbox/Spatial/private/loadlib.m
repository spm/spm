function loadlib(nam)
if nargin<1, nam = 'pushpull'; end

if ~libisloaded(nam)
    pth = mfilename('fullpath');
    [d,~,~] = fileparts(pth);
    mext = mexext;
    if ispc, soext = 'dll'; else soext = 'so'; end
    libname = fullfile(d,"..","lib",['lib' mext(4:end)],[nam '.' soext]);
    hdrname = fullfile(d,"..","C",  [nam '.h']);
    [notfound,warnings] = loadlibrary(libname,hdrname);
    if ~isempty(notfound)
        fprintf('\n --- Displaying debugging information ---\n');
        disp(notfound)
        disp(warnings)
        fprintf('\n ---                                  ---\n');
    end
end

