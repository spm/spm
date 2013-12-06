The averaged traces are in a structure called allsubj (load DATA.mat)
Subfields are:
FN = Fast Noisy (i.e., 22 deg/sec)
FS = Fast Smooth
SN = Slow Noisy (i.e., 18 deg/sec)
SS = Slow Smooth

Within these subfields: ‘.one’ contains the averaged eye position, measured in pixels (from 600 to -600) over several thousand milliseconds. 

The target paths are in the ‘target’ structures target18.one and target22.one

The occluder was present between target18.plotocc.x(1,1) and (1,2), and between (2,1) and (2,2), in milliseconds.


See spm_SEM_demo.m for an illustration of how these data are inverted




