---
title: 'SPM 25: open source neuroimaging analysis software'
tags:
  - neuroimaging
  - neuroscience
  - SPM
  - MRI
  - MEG
  - EEG
  - PET
authors:
  - name: Tim M. Tierney
    affiliation: 1	
  - name: Nicholas Alexander
    affiliation: 1	
  - name: John Ashburner
    affiliation: 1	
  - name: Nicole Labra Avila
    affiliation: 1	
  - name: Yaël Balbastre
    affiliation: 1	
  - name: Gareth Barnes
    affiliation: 1	
  - name: Yulia Bezsudnova
    affiliation: 1	
  - name: Mikael Brudfors
    affiliation: 1	
  - name: Korbinian Eckstein
    affiliation: 2	
  - name: Guillaume Flandin
    affiliation: 1	
  - name: Karl Friston
    affiliation: 1	
  - name: Amirhossein Jafarian
    affiliation: 1		
  - name: Olivia Kowalczyk
    affiliation: 1	
  - name: Vladimir Litvak
    affiliation: 1	
  - name: Johan Medrano
    affiliation: 1	
  - name: Stephanie Mellor
    affiliation: 1	
  - name: George O'Neill
    affiliation: 1	
  - name: Thomas Parr
    affiliation: 1	
  - name: Adeel Razi
    affiliation: 1	
  - name: Ryan Timms
    affiliation: 1	
  - name: Peter Zeidman
    orcid: 0000-0003-3610-6619
    corresponding: true
    affiliation: 1
affiliations:
 - name: University College London, UK
   index: 1
   ror: 02jx3x895
 - name: Independent Researcher, Germany
   index: 2
date: 18 December 2024
bibliography: paper.bib

---

# Summary

Statistical Parametric Mapping (SPM) is an integrated set of methods for testing hypotheses about the brain's structure and function, using data from medical imaging devices. These methods are implemented in an open source software package, `SPM`, which has been in continuous development for more than 30 years by an international community of developers. This paper reports the release of `SPM 25.01`, a major new version of the software that incorporates novel analysis methods, optimisations of existing methods, as well as improved practices for open science and software development.

# Statement of need

`SPM` introduced many of the statistical methods that underpin cognitive and clinical neuroimaging research today, including:

- The voxel-wise application of General Linear Models (GLMs) to neuroimaging data [@friston1994statistical].
- Convolution modelling of functional MRI (fMRI) signals using haemodynamic response functions [@friston1994analysis].
- Correction for multiple comparisons using topological methods (Random Field Theory, RFT) [@worsley1996unified].
- Event-related fMRI [@josephs1997event].
- Voxel-Based morphometry (VBM) for detecting changes in anatomy [@ashburner2000voxel].
- Dynamic Causal Modelling (DCM) for state-space modelling using variational Bayesian methods [@friston2003dynamic].
- Source localisation for M/EEG data using variational Bayesian methods [@phillips2005empirical].

These methods share certain key principles: the use of generative models, the application of well-motivated parametric statistics and a commitment to open science practices. They are included in a major new release of `SPM`, which addresses a series of needs in the neuroimaging community, set out below.

## Open development

`SPM` was previously developed and tested using a private Subversion server within University College London. To enable community engagement in the future development of `SPM` and to increase transparency, development has recently moved to a public [Github repository](https://github.com/spm/spm). `SPM 25.01` is the first release of the software following the move to Github. The key advantages of using Github thus far have been:

- Introducing automated unit and regression tests across platforms.
- Automating the build process to conveniently generate and release source code and compiled versions.
- Issue tracking and distributing tasks among developers.

## Documentation and training

The documentation for `SPM` was previously spread across multiple locations, most of which could not be edited by the community. `SPM 25.01` is accompanied by a new [documentation website](https://www.fil.ion.ucl.ac.uk/spm/docs/), the source code for which is hosted in a public [Github repository](https://github.com/spm/spm-docs). The new website has step-by-step tutorials on all of SPM's main features, as well as freely available video recordings of lectures from previous SPM courses covering the mathematical theory.

## Major new features 

`SPM 25` includes 10 years of new developments since the last major release (`SPM 12`, dated 1st October 2014). This section highlights some of the most significant new features for different neuroimaging modalities.

### MRI

- Multi-Brain Toolbox [@brudfors2020flexible]. Generates population average-shaped brains, enabling more precise spatial normalisation with the option to automatically label brain structures.

- SCOPE Toolbox. Generates MRI fieldmaps using phase-encode-reversed pairs of MRI images (blip-up and blip-down images), similar to the Topup tool in FSL.

### M/EEG

- Methods for spectral decomposition - `SPM 25.01` offers an implementation of an existing approach called FOOOF (Specparam) in the MEEGtools toolbox, based on code from Brainstorm [@donoghue2020parameterizing], as well as a new Bayesian implementation that introduces formal statistical testing, called Bayesian Spectral Decomposition (BSD) [@medrano2024bsd].

- Support for fusion of different MEG sensor types and EEG sensors in beamforming with pre-whitening [@westner2022unified].

- Support for MEG BIDS for specification of events, channels and fiducials [@westner2022unified].

### Bayesian statistics

- Parametric Empirical Bayes (PEB) [@friston2016bayesian] extends the Dynamic Causal Modelling (DCM) framework to include random effects modelling of neural connectivity parameters, enabling people to test hypotheses about the similarities and differences among research participants.

- Bayesian model reduction (BMR) [@friston2018bayesian] enables statistical evidence to be rapidly scored for large numbers of competing models, where models differ only in their priors.

### OPMs

A major recent innovation in neuroimaging is MEG using Optically Pumped Magnetometers (OPMs), which enable free movement of the head and body during neural recordings [@boto2018moving]. This makes MEG available to new experimental paradigms (e.g., experiments involving free movement [@mellor2023real]), new body parts (e.g., the spinal cord [@spedden2024towards]) and new study populations who may not be amenable to traditional MEG (e.g., people with epilepsy [@mellor2024detection]). Developing analysis tools for OPM data is a major focus for SPM, with recently added features including:

- File IO for all major OPM manufacturers (Quspin, Cerca, Mag4Health, Fieldline).
- Methods to simulate arbitrary OPM arrays of differing densities and vector measurements.
- OPM interference cancellation algorithms for low channel systems: Homogeneous Field Correction [@tierney2021modelling].
- OPM interference cancellation algorithms for large channel systems: Adaptive Multipole Models [@tierney2024adaptive].

### SPM without MATLAB

Approximately 90% of the `SPM 25.01` source code is written in MATLAB and the remainder is C++. This code has been highly optimised and thoroughly tested over 30 years of development. We have therefore carefully considered how to capitalise on the stability of the `SPM` software, while making it more  accessible for people who do not have access to a MATLAB license, or who prefer to write their analysis code in other languages. 

Our strategy is as follows:

- `SPM 25` will be the first version of SPM to be fully accessible from the Python programming language, without requiring MATLAB, using a new Python wrapper called [spm-python](https://github.com/spm/spm-python). This is in the final stages of development and will be released in the first quarter of 2025.

- [SPM Standalone](https://www.fil.ion.ucl.ac.uk/spm/docs/installation/standalone/) is the compiled version of `SPM` that can be run from the command line without a MATLAB license. This enables people to run neuroimaging analyses from command line scripts written in any language, or using the GUI. It is now generated automatically with each new release, as part of the Github-based build process. 

- [Docker and Singularity containers](https://www.fil.ion.ucl.ac.uk/spm/docs/installation/containers/) are additionally provided and are now generated automatically as part of SPM's Github build process.

# Software versions

`SPM 25.01` is the first release of `SPM` to use calendar versioning, thus `SPM 25.01` is the version issued in January 2025. All releases are available via [https://github.com/spm/spm/releases](https://github.com/spm/spm/releases). 

# Acknowledgements

Tim M. Tierney is funded by an Epilepsy Research UK fellowship (FY2101). Johan Medrano is supported by the Discovery Research Platform for Naturalistic Neuroimaging funded by Wellcome
[226793/Z/22/Z]. Peter Zeidman is funded by an MRC Career Development Award [MR/X020274/1]. A full list of authors of `SPM 25.01` can be found in the file `AUTHORS.txt` supplied with the software.

# References