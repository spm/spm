# Short explanation of principles behind FieldMap toolbox for SPM

## 1. Short explanation about why field-disturbances cause distortions in EPI images

We can think of it in the following way; Each location (voxel) in
object will emit a series of echoes with its own characteristic 
"fingerprint". This "fingerprint" will consist of the frequency
within the echoes, which will uniquely identify its position in
the frequency encode direction (typically left->right). It will
also consist of the phase difference between consecutive echoes
(caused by the phase-encode blips) that identifies its position
in the phase-encode direction (typically ant->post).

So, let us for example say we have a voxel to the right of 
the object, and in the middle in the up-down direction. This 
voxel will then experience a higher field in the readout 
phase (lets ignore positive/negative readouts for the time 
being) and no field changes due to the phase-encode blips. 
To assign this signal to the right voxel the "reconstruction" 
will then search for a signal with high frequency during 
read-out, and whose phase doesn't change from one echo to 
the other. It does so by taking the inner product between 
the echo-train and a vector that has this (expected) behaviour.
Let us now consider what happens if this voxel (in the absence 
of gradients) experiences a slightly higher  field than we "think" 
(and we'll soon discuss why). This means that during read-out it 
will emit a signal with slightly higher frequency than we would 
expect, and if the "reconstruction" is unaware of this, it will
"misplace" the intensity to the right of its proper place. The
distance by which intensity is misplaced depends on the relative
magnitude of the read-out gradient and the "field disturbance".
Typically the gradient is very strong relative to the "disturbance"
and left->right misplacement is negligible.

However, the slightly higher than tentative frequency also means 
that there is a continuous phase-development (remember that the
"disturbance" is "on" all the time) which will lead to a phase-
difference between consecutive echoes. Remember that position in
the up->down direction was encoded in this phase-difference, which
means that the "reconstruction" will assign the intensity also to
the wrong position in that direction. This time the encoding 
gradient has an "unfair" disadvantage since it is on only for a
very short time (the "blip") whereas the "disturbance" will affect
the phase throughout the whole readout as well. Therefore the 
displacements in the phase-encode direction can be quite
substantial, and it may be a good idea to do something about it.
The problem can be alleviated by faster gradients, which will
shorten the time for read-out and hence the time during which
the "disturbance" is allowed to do its dirty deed.
It can also be alleviated by using bigger blips (i.e.
changing the balance between the blips and the "disturbance").
This would then mean lower resolution for a given FOV, UNLESS
additional means of spatial encoding is used (i.e. SMASH/SENSE).
Hence, for the latest generation scanners this is less of a 
problem than it used to be, unless of course one wants to start
pushing resolution.

## 2. Short explanation of why there are field disturbances

Different types off materials are more or less susceptible to
magnetisation by an external field. Air for example will not 
"resist" magnetisation much at all, and the resulting field
M will be very close to the external field B. On the other 
hand water (for example) will put up a bit of a fight (i.e. is
less susceptible) and the resulting field will be slightly lower.
In the vicinity of junctions between different materials (i.e.
around air-cavities in the head) the field will vary in a quite
non-intuitive manner that depends on the size and shape of the
cavity, and its orientation relative to the magnetic flux (which
goes in the z-direction in an MR-scanner).

This means that the problems are located mainly in the frontal
lobe (caused by the sinuses), the orbito-frontal cortex 
(caused by the roof of the palate) and the temporal lobes
(the ear canals). 

Another source of susceptibility distortions (though more
rarely discussed) is the blood in large vessels, as can be
observed around the sagittal sinus.

Since the field (and hence the distortions) depends on the   
orientation of the cavities (and hence the object) the field
and the ensuing distortions will depend on (be a function of)
the exact position of the subject.

## 3. How can we measure the field?

If we know the field it is an easy task to re-sample the 
distorted images such that the misplaced intensity is "replaced".
This can be done either in (i)mage- or k-space, but for reasons
of convenience is mostly done in i-space.

Let us now consider a non-EPI sequence, which could very 
simplistically be described as:
```
for all PE-gradients
{ 
  Excitation - PE-gradient on-off - Simultaneous FE-gradient readout
}
```
In this sequence, the time between excitation and the readout of
each echo is exactly the same for each and every echo, which 
means that any field "disturbance" will affect each echo in
exactly the same way (i.e. it has the same time to influence
the signal). Therefore, in contrast to EPI, the relative phase
between consecutive echos (which is what encodes position in the
PE-direction) is _completely_ unaffected.

However, the overall phase of the data is still affected (i.e. a
spatially varying phase is introduced, but one which is common
for all PE-steps. Hence, we should in principle be able to deduce
local field-strength from the phase of the reconstructed images.
In practise this doesn't work since there are numerous other 
factors (mainly timing errors) that affect the phase. 

The easiest way around this is to perform two acquisitions with
slightly different echo-times. The other factors are not affected
by this, so any difference in phase between the two images can
be attributed to a longer time for (disturbance induced) 
phase to evolve when the echo-time is longer. Hence, the
difference in phase between two images acquired with different 
times will tell us how fast the phase changes. I.e. we can get
a frequency-map from the phase difference (divided by 2pi) 
divided by the echo-time difference. From the frequency map
one can easily get a field-map through the Larmor-equation.
Well then, if it is that easy why don't we all do it?

## 4. Ifs and buts I: Phase-wrapping

Let us say we do two measurements with an echo time difference
of 10ms and that we in a given voxel observe a phase of 20degrees
in the first image and 80degrees in the second image. We will
deduce that the phase has evolved 60degrees in 10ms, and hence
that the field is off by (60/360)/10e-3 ~ 17Hz compared to what
we "think". The only problem is that exactly the same result
would have been obtained if the field had been 
```
((60+360)/360)/10e-3 ~ 117Hz off, or 
((60-360)/360)/10e-3 ~ -83Hz off etc etc.
```
So, we have an ambiguity here and one that isn't easily resolved
because of its non-linear nature.

One solution (the one adopted by the scanner manufacturers, cause
they don't like post-processing) is to use such a short echo-time
difference that it is highly unlikely that one will ever obtain
a |phase-difference| larger than pi (180degrees). The problem with
that approach is that most values that we observe will then be
around a few degrees, i.e. in the noise range, which means we will
not have very accurate values for most parts of the brain.
The other solution is to use "phase-unwrapping". This implies 
deducing on the true phase by considering also the phase of
neighbouring voxels. Say for example that we observe the 
following phase along a row of 10 voxels.
```
 50  95 150 -160 -150 -170 160 130 110 100
```
Clearly most of us would assume that the "true"  phase-
differences were rather 
```
 50  95 150  200  210  190 160 130 110 100
```
Hence, the simplest possible "phase-unwrapping" would simply
proceed along rows (or columns) and whenever it observed a 
jump in phase from one voxel to the next greater than pi
(180degrees) it would add (or subtract) 2pi (360degrees) to
the present voxel.

The problem with that is that once a mistake has been made
(i.e. deducing a wrap when there is none, or vice versa)
all the voxels downstream will be affected by that error.
Mistakes will typically happen when passing through regions
where there is little signal (and hence low SNR for the
phase) e.g. areas of air or bone.

As it turns out in practise, this will almost always happen
if one goes about the phase-unwrapping in a "naive" way.
There are (at least) two possible "non-naive" ways of performing 
the phase unwrapping to avoid this problem.

**4.1** One is to unwrap from some seed point in a high SNR voxel
into areas with gradually decreasing SNR. That way any mistakes
are likely to be made late in the process, and should hence have
less influence. The implementation will use a water-shed like
algorithm which can be thought of as pouring water into a basin
where the depth is inversely proportional to the SNR. This is the
"traditional" approach that has been taken to unwrapping.
The problem with this type of solution is that it is non-trivial
to implement a fast water-shed algorithm, and that the final
result will depend critically on the step-size of our filling
of the basin.

**4.2** Another is to divide the brain/object into a set of (3D)
connected regions, where within each region the phase lies
within some smallish interval (e.g. pi/4). This guarantees
that there are no phase-wraps within a region, and hence that
all "wrapping surfaces" run along surfaces between regions.
One can then start to merge regions according to a cunning
scheme, and for each merge wrap (or not) one of the regions
depending on the average phase-difference along the surface
between the regions. This approach has been suggested by Mark
Jenkinson from the evil fMRIB empire. It is also non-trivial
to implement.

In this toolbox we have implemented both approaches. The
comparisons we have made indicate that the latter approach
is faster and much more robust, and is hence the default choice.

## 5. Ifs and buts II: Choice of phase-mapping sequence

The description of a phase-mapping sequence above was of a 
non-EPI based sequence. The same measurement could be performed
using a dual echo-time EPI sequence, in which case the fieldmap 
is obtained in the distorted space (as opposed to "true" space as
is the case of non-EPI). This is in itself not a problem and all
we need to do is to invert the EPI-fieldmap prior to using it (which
can be thought of as correcting it by itself). Hence, if we measure
two fieldmaps in the same subject, once with an EPI sequence and 
once with a non-EPI sequence we would expect them to be identical
after inversion of the EPI-base one.

They are not.

It is our experience that the difference between a non-EPI and 
an EPI based fieldmap is much greater than the within-sequence
difference, and it is not at all clear to us why this should be
the case.

It is consequently not clear to us which one represents the "truth"
and hence should be used. Arguments for and against using an EPI
based fieldmap would be.

**5.1** For

1. It is based on the same sequence that we want to correct and 
should therefore be affected by the same factors (e.g. Maxwell
effects).
2. The magnitude images underlying the field-map are similar to
those we want to correct, which renders alignment very simple.
3. It is quick.

**5.2** Against

1. It necessitates inversion of the fieldmap which means that
local errors (e.g. due to failed unwrapping) will be propagated
to other parts of the map.

I think on the whole we cautiously recommend using an EPI-based
sequence.

## 6. The field changes when subject moves

As described above the susceptibility induced field changes are 
dependent on the exact position of the object. This means that as
the subject moves the field, and hence also the apparent shape of
the object, changes. The relative shape-changes resulting from
typical movements are small compared to the "overall" distortion.
Hence it is not a big issue in terms of anatomical fidelity.
However, it will result in residual movement-related variance in
the data.

In order to correct this one would in principle have to acquire
one field-map for each point of the time-series. This would double
the TR which would rarely be acceptable. In  addition it has been
shown by Hutton et al. that such a strategy may actually 
_introduce_ variance in time-series where movement are small 
(typical).

One solution to this would be to combine Unwarp (see main SPM2 
documentation for a description of principles behind Unwarp) with
one measured fieldmap per subject. Unwarp utilises the observed 
shape changes in the time-series to estimate the rate of change of
the field w.r.t. subject movements. It CANNOT however estimate the
field per se, since that doesn't introduce any variance across the
series. 

By using the FieldMap toolbox to calculate a "static" field (i.e. a 
field that contributes equally to all images in the timeseries)
and then use this field map with Unwarp to calculate the
movement-induced contributions to each time point it should be
possible to correct for (almost) all adverse effects of
susceptibility induced distortions.

## 7. Who do I credit when using the toolbox?

Primarily Chloe and Jesper by buying drinks whenever you encounter
either at a conference or similar.

Academically the implementations in FieldMap/Unwarp are based on
the following papers:

1. The original idea to use fieldmaps for distortion correction:

>  Jezzard P & Balaban RS. 1995. Correction for geometric distortion in
>  echo planar images from Bo field variations. MRM 34:65-73.

2. There isn't really a paper describing exactly our implementation
     of watershed based unwrapping. An example of a paper about
     the initial Fieldmap implementation and use of a similar algorithm is:

>  Hutton C, Bork A, Josephs O, Deichmann R, Ashburner J,
>  Turner R. 2002. Image distortion correction in fMRI: A
>  quantitative evaluation. NeuroImage 16:217-240.

3. FieldMap when using the Mark3D (default) or Mark2D options

>  Jenkinson M. 2003. Fast, automated, N-dimensional phase-
>  unwrapping algorithm. MRM 49:193-197.

4. Unwarp (i.e. the Realign & Unwarp option in SPM2)

>  Andersson JLR, Hutton C, Ashburner J, Turner R, Friston K.
>  2001. Modelling geometric deformations in EPI time series.
>  NeuroImage 13:903-919.

5. Combined use of FieldMap and Unwarp (fully implemented in SPM5)

>  Hutton C, Deichmann R, Turner R, Andersson JLR. 2004. Combined 
>  correction for geometric distortion and its interaction with 
>  head motion in fMRI. Proceedings of ISMRM 12, Kyoto, Japan.

```
Good Luck! Jesper Andersson & Chloe Hutton
```
