mex -v CFLAGS="-ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread -Wall" shoot3.c shoot_diffeo3d.c shoot_optim3d.c shoot_multiscale.c shoot_regularisers.c shoot_expm3.c shoot_dartel.c shoot_boundary.c -O -largeArrayDims -DNEUMANN

