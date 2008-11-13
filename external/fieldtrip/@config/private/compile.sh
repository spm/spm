#
# shell helper script to compile all mex files
#

if   [ -a /opt/matlab72/bin/mex ]; then
  MEX=/opt/matlab72/bin/mex
elif [ -a /Applications/MATLAB72/bin/mex ]; then
  MEX=/Applications/MATLAB75/bin/mex
elif [ -a /Applications/MATLAB-2007b/bin/mex ]; then
  MEX=/Applications/MATLAB-2007b/bin/mex
else
  echo could not locate mex compiler
  exit 1
fi

$MEX increment.c
$MEX reset.c
$MEX deepcopy.c

